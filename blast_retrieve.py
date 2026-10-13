#!/usr/bin/env python3
"""Step 1 of the frozen retrieval protocol: broad homolog retrieval.

One search per structural reference, each against refseq_protein. RefSeq is
chosen over nr deliberately: WP_ accessions resolve to genomic locations through
NCBI's Identical Protein Groups, which is what the synteny stage needs. Searching
nr would return many sequences with no usable genomic context.

Nothing is filtered here beyond BLAST's own E-value. The 30-90% identity window
and every other gate is applied downstream, so the raw retrieval stays auditable.

THE SEARCH PARAMETERS BELOW ARE FROZEN. This rewrite (2026-10-12) changed the
job-management infrastructure only: program, database, E-value, hit-list size and
filtering are byte-identical to the superseded version, which is kept as
`blast_retrieve_superseded_20261012.py`. The state file records them and this
script refuses to resume a state written under different parameters.

WHY IT WAS REWRITTEN
--------------------
Two retrieval attempts (2026-10-10, 2026-10-11) were lost. Both submitted all
five searches at once, held the RIDs only in memory, polled on a fixed 60 s
cycle, and died with the terminal. NCBI RIDs expire after about 24 h, so the
results -- if they ever landed -- were unrecoverable, and neither run ever wrote
a single result file.

DESIGN
------
  * every RID is written to `retrieval/state.json` the moment it is issued,
    with its family, submission time and status;
  * a rerun resumes RIDs that are still valid and resubmits the ones that have
    expired (by NCBI status UNKNOWN, or by age past RID_TTL_H);
  * searches are submitted ONE AT A TIME -- a family is not submitted until the
    previous one has reached a terminal state (--max-concurrent to override);
  * each result is validated (XML parses, blastp, hits present) and written
    atomically the moment it arrives, so an interruption cannot truncate it;
  * polling is bounded: per-invocation wall budget, per-job exponential backoff,
    and explicit terminal states;
  * READY, WAITING, FAILED, UNKNOWN, network failure and malformed response are
    distinguished and recorded separately;
  * state is flushed after every transition and on SIGINT/SIGTERM, so a closed
    terminal or a sleeping machine costs nothing but the resubmission of
    whatever was still queued;
  * "RETRIEVAL COMPLETE" is printed only when all five families hold a validated
    result or an explicitly documented failure. Anything else exits non-zero and
    names what is outstanding.

USAGE
-----
  python3 blast_retrieve.py --status                 # state only, no network
  python3 blast_retrieve.py --families r2pyocin      # single-family test
  python3 blast_retrieve.py --max-minutes 240        # resume / run the rest
  python3 blast_retrieve.py --max-concurrent 1       # default: strictly serial

Exit codes: 0 all families terminal (results or documented failures)
            2 incomplete but resumable -- rerun the same command
            3 terminal with at least one abandoned family
Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import signal
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "retrieval"
STATE_PATH = OUT / "state.json"
LOG_PATH = HERE / "retrieval.log"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

# ---- FROZEN search parameters -- do not edit ---------------------------------
PROGRAM = "blastp"
DATABASE = "refseq_protein"
EXPECT = 1e-5
HITLIST = 1000
FILTER = "F"

# ---- infrastructure parameters -- tunable, no scientific content ------------
SUBMIT_GAP = 12         # s, NCBI asks for >=10 s between submissions
POLL_MIN = 60           # s, NCBI asks for >=60 s between polls of one job
POLL_MAX = 600          # s, backoff ceiling
POLL_GROWTH = 1.5
RID_TTL_H = 20.0        # h, resubmit before NCBI's ~24 h expiry
MAX_ATTEMPTS = 3        # submissions per family before it is abandoned
MAX_CONSECUTIVE_ERRORS = 8
MIN_RESULT_BYTES = 10_000
HTTP_TRIES = 3
HTTP_BACKOFF = (5, 15, 45)
USER_AGENT = "sheath-evo-pilot/2.0"

TERMINAL = {"done", "no_hits", "abandoned"}
SCHEMA = 2

_stop = False


def now() -> datetime:
    return datetime.now(timezone.utc)


def iso(dt: datetime | None = None) -> str:
    return (dt or now()).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_iso(s: str) -> datetime:
    return datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def log(msg: str) -> None:
    line = f"{iso()}  {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a") as fh:
        fh.write(line + "\n")


class NetworkFailure(Exception):
    """Transport-level failure: no answer from NCBI, or an HTTP error."""


class Malformed(Exception):
    """NCBI answered, but not with something this protocol recognises."""


# ---------------------------------------------------------------- state -------

def atomic_write_text(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def atomic_write_bytes(path: Path, data: bytes) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def frozen_params() -> dict:
    return {"program": PROGRAM, "database": DATABASE, "expect": EXPECT,
            "hitlist": HITLIST, "filter": FILTER}


def blank_family() -> dict:
    return {"status": "unsubmitted", "rid": None, "submitted_at": None,
            "last_poll_at": None, "poll_count": 0, "backoff_s": POLL_MIN,
            "consecutive_errors": 0, "attempts": [], "result": None, "note": None}


def load_state(families: list[str]) -> dict:
    if STATE_PATH.exists():
        state = json.loads(STATE_PATH.read_text())
        if state.get("schema") != SCHEMA:
            sys.exit(f"ABORT: {STATE_PATH} has schema {state.get('schema')}, "
                     f"expected {SCHEMA}. Inspect it by hand.")
        if state.get("params") != frozen_params():
            sys.exit(f"ABORT: {STATE_PATH} was written under different search "
                     f"parameters:\n  on disk: {state.get('params')}\n"
                     f"  current: {frozen_params()}\nThe frozen protocol has "
                     "changed or the state belongs to another run. Not resuming.")
    else:
        state = {"schema": SCHEMA, "params": frozen_params(),
                 "created_at": iso(), "families": {}}
    for fam in families:
        state["families"].setdefault(fam, blank_family())
    return state


def save_state(state: dict) -> None:
    state["updated_at"] = iso()
    atomic_write_text(STATE_PATH, json.dumps(state, indent=2, sort_keys=True))


# ---------------------------------------------------------------- http --------

def request(params: dict, method: str) -> str:
    """One NCBI call, with bounded retries. Raises NetworkFailure on exhaustion."""
    last = None
    for attempt in range(HTTP_TRIES):
        try:
            if method == "POST":
                req = urllib.request.Request(
                    BLAST, data=urllib.parse.urlencode(params).encode(),
                    headers={"User-Agent": USER_AGENT})
            else:
                req = urllib.request.Request(
                    BLAST + "?" + urllib.parse.urlencode(params),
                    headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=300) as fh:
                if fh.status != 200:
                    raise NetworkFailure(f"HTTP {fh.status}")
                return fh.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as exc:
            last = f"HTTPError {exc.code}"
        except Exception as exc:                      # URLError, timeout, reset
            last = f"{type(exc).__name__}: {exc}"
        if attempt + 1 < HTTP_TRIES:
            log(f"      transport failure ({last}); retry in {HTTP_BACKOFF[attempt]}s")
            time.sleep(HTTP_BACKOFF[attempt])
    raise NetworkFailure(last or "unknown transport failure")


def submit(name: str, seq: str) -> str:
    txt = request({"CMD": "Put", "PROGRAM": PROGRAM, "DATABASE": DATABASE,
                   "QUERY": seq, "HITLIST_SIZE": HITLIST, "EXPECT": EXPECT,
                   "FILTER": FILTER}, "POST")
    m = re.search(r"^\s*RID = (\S+)", txt, re.M)
    if not m:
        raise Malformed(f"no RID in Put response for {name} "
                        f"({len(txt)} bytes, starts {txt[:120]!r})")
    return m.group(1)


def poll(rid: str) -> tuple[str, str]:
    """-> (state, detail). state in ready|no_hits|waiting|failed|unknown."""
    txt = request({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}, "GET")
    m = re.search(r"Status=(\w+)", txt)
    if not m:
        raise Malformed(f"no Status field for {rid} ({len(txt)} bytes)")
    status = m.group(1).upper()
    if status == "READY":
        hits = re.search(r"ThereAreHits=(\w+)", txt)
        if hits and hits.group(1).lower() == "yes":
            return "ready", "READY, hits present"
        return "no_hits", "READY, ThereAreHits=no"
    if status in ("WAITING", "SEARCHING", "RUNNING"):
        return "waiting", status
    if status == "FAILED":
        return "failed", "NCBI reports the search itself failed"
    if status == "UNKNOWN":
        return "unknown", "RID expired or not known to NCBI"
    raise Malformed(f"unrecognised Status={status} for {rid}")


def fetch(rid: str) -> str:
    return request({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
                    "ALIGNMENTS": HITLIST, "DESCRIPTIONS": HITLIST}, "GET")


def validate(xml: str, rid: str) -> tuple[int, str]:
    """-> (n_hits, query_def). Raises Malformed on anything unusable."""
    if len(xml.encode()) < MIN_RESULT_BYTES:
        raise Malformed(f"{rid}: result is {len(xml.encode())} bytes, "
                        f"below the {MIN_RESULT_BYTES} floor")
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as exc:
        raise Malformed(f"{rid}: XML does not parse ({exc})") from exc
    if root.tag != "BlastOutput":
        raise Malformed(f"{rid}: root element is {root.tag!r}, not BlastOutput")
    prog = (root.findtext("BlastOutput_program") or "").strip()
    if prog != PROGRAM:
        raise Malformed(f"{rid}: program is {prog!r}, expected {PROGRAM!r}")
    db = (root.findtext("BlastOutput_db") or "").strip()
    if DATABASE not in db:
        raise Malformed(f"{rid}: database is {db!r}, expected {DATABASE!r}")
    n_hits = len(root.findall(".//Hit"))
    if n_hits == 0:
        raise Malformed(f"{rid}: status said hits were present but the XML "
                        "contains no Hit elements")
    return n_hits, (root.findtext("BlastOutput_query-def") or "").strip()


# ---------------------------------------------------------------- driver ------

def reconcile_from_disk(state: dict) -> None:
    """A validated XML on disk is the truth, whatever the state file says."""
    for fam, f in state["families"].items():
        path = OUT / f"{fam}.xml"
        if f["status"] == "done":
            if not path.exists():
                log(f"  {fam}: state says done but {path.name} is missing "
                    "-> back to unsubmitted")
                state["families"][fam] = blank_family()
            continue
        if path.exists():
            try:
                n_hits, _ = validate(path.read_text(), f.get("rid") or path.name)
            except Malformed as exc:
                log(f"  {fam}: {path.name} present but unusable ({exc}); ignoring it")
                continue
            data = path.read_bytes()
            f.update(status="done", note="result already on disk at startup",
                     result={"path": path.name, "bytes": len(data),
                             "sha256": hashlib.sha256(data).hexdigest(),
                             "n_hits": n_hits, "fetched_at": f.get("last_poll_at")})
            log(f"  {fam}: already retrieved ({n_hits} hits) -> done")


def expired(f: dict) -> bool:
    if not f.get("submitted_at"):
        return False
    age_h = (now() - parse_iso(f["submitted_at"])).total_seconds() / 3600.0
    return age_h > RID_TTL_H


def close_attempt(f: dict, outcome: str, note: str) -> None:
    if f.get("rid"):
        f["attempts"].append({"rid": f["rid"], "submitted_at": f["submitted_at"],
                              "closed_at": iso(), "polls": f["poll_count"],
                              "outcome": outcome, "note": note})
    f["note"] = note


def do_submit(state: dict, fam: str, seq: str) -> None:
    f = state["families"][fam]
    if len(f["attempts"]) >= MAX_ATTEMPTS:
        f.update(status="abandoned",
                 note=f"{MAX_ATTEMPTS} submissions all expired or failed; "
                      "not resubmitted automatically")
        log(f"  {fam}: ABANDONED after {MAX_ATTEMPTS} attempts")
        save_state(state)
        return
    try:
        rid = submit(fam, seq)
    except (NetworkFailure, Malformed) as exc:
        f["consecutive_errors"] += 1
        f["note"] = f"submission failed: {exc}"
        log(f"  {fam}: submission failed ({exc})")
        save_state(state)
        return
    f.update(status="waiting", rid=rid, submitted_at=iso(), last_poll_at=None,
             poll_count=0, backoff_s=POLL_MIN, consecutive_errors=0,
             note="submitted, awaiting results")
    save_state(state)                      # RID persisted before anything else
    log(f"  {fam}: submitted RID={rid} (attempt {len(f['attempts']) + 1}"
        f"/{MAX_ATTEMPTS}) -> state.json")
    time.sleep(SUBMIT_GAP)


def do_poll(state: dict, fam: str) -> None:
    f = state["families"][fam]
    if expired(f):
        close_attempt(f, "expired_by_age", f"RID {f['rid']} passed {RID_TTL_H} h "
                                           "without returning; resubmitting")
        f.update(status="unsubmitted", rid=None, submitted_at=None)
        log(f"  {fam}: RID aged out past {RID_TTL_H} h -> will resubmit")
        save_state(state)
        return

    try:
        st, detail = poll(f["rid"])
    except NetworkFailure as exc:
        f["consecutive_errors"] += 1
        f["last_poll_at"] = iso()
        f["backoff_s"] = min(int(f["backoff_s"] * POLL_GROWTH), POLL_MAX)
        f["note"] = f"transport failure on poll: {exc}"
        log(f"  {fam}: network failure {f['consecutive_errors']}"
            f"/{MAX_CONSECUTIVE_ERRORS} ({exc})")
        if f["consecutive_errors"] >= MAX_CONSECUTIVE_ERRORS:
            close_attempt(f, "network", f"{MAX_CONSECUTIVE_ERRORS} consecutive "
                                        f"transport failures; last: {exc}")
            f["status"] = "abandoned"
            log(f"  {fam}: ABANDONED -- NCBI unreachable")
        save_state(state)
        return
    except Malformed as exc:
        f["consecutive_errors"] += 1
        f["last_poll_at"] = iso()
        f["backoff_s"] = min(int(f["backoff_s"] * POLL_GROWTH), POLL_MAX)
        f["note"] = f"malformed status response: {exc}"
        log(f"  {fam}: MALFORMED status response ({exc})")
        save_state(state)
        return

    f["poll_count"] += 1
    f["last_poll_at"] = iso()
    f["consecutive_errors"] = 0

    if st == "waiting":
        f["backoff_s"] = min(int(f["backoff_s"] * POLL_GROWTH), POLL_MAX)
        f["note"] = f"waiting ({detail})"
        log(f"    {fam} ({f['rid']}): waiting  [poll {f['poll_count']}, "
            f"next in {f['backoff_s']}s]")
        save_state(state)
        return

    if st == "no_hits":
        close_attempt(f, "no_hits", detail)
        f["status"] = "no_hits"
        log(f"  {fam}: NO HITS -- terminal, documented")
        save_state(state)
        return

    if st == "failed":
        close_attempt(f, "failed", detail)
        f.update(status="unsubmitted", rid=None, submitted_at=None)
        log(f"  {fam}: NCBI reports FAILED -> will resubmit")
        save_state(state)
        return

    if st == "unknown":
        close_attempt(f, "expired", detail)
        f.update(status="unsubmitted", rid=None, submitted_at=None)
        log(f"  {fam}: RID UNKNOWN (expired) -> will resubmit")
        save_state(state)
        return

    # ---- READY with hits: fetch, validate, persist atomically ---------------
    log(f"    {fam} ({f['rid']}): READY -- fetching")
    try:
        xml = fetch(f["rid"])
        n_hits, qdef = validate(xml, f["rid"])
    except (NetworkFailure, Malformed) as exc:
        f["consecutive_errors"] += 1
        f["backoff_s"] = min(int(f["backoff_s"] * POLL_GROWTH), POLL_MAX)
        f["note"] = f"READY but fetch/validate failed: {exc}"
        log(f"  {fam}: fetch or validation failed ({exc}); will retry")
        if f["consecutive_errors"] >= MAX_CONSECUTIVE_ERRORS:
            close_attempt(f, "unfetchable",
                          f"READY but never fetched cleanly; last: {exc}")
            f["status"] = "abandoned"
            log(f"  {fam}: ABANDONED -- READY but result unusable")
        save_state(state)
        return

    data = xml.encode()
    path = OUT / f"{fam}.xml"
    atomic_write_bytes(path, data)
    f.update(status="done", note=f"validated: {n_hits} hits",
             result={"path": path.name, "bytes": len(data),
                     "sha256": hashlib.sha256(data).hexdigest(),
                     "n_hits": n_hits, "query_def": qdef, "fetched_at": iso()})
    close_attempt(f, "done", f"{n_hits} hits")
    f["status"] = "done"
    save_state(state)
    log(f"  {fam}: saved {path.name} ({len(data)} bytes, {n_hits} hits, "
        f"sha256 {f['result']['sha256'][:12]}...)")


def report(state: dict, families: list[str]) -> int:
    rows, outstanding, failures = [], [], []
    for fam in families:
        f = state["families"][fam]
        res = f.get("result") or {}
        rows.append(f"  {fam:12s} {f['status']:12s} {str(f.get('rid') or '-'):14s} "
                    f"{str(res.get('n_hits', '-')):>6s} hits  "
                    f"attempts={len(f['attempts'])}  {f.get('note') or ''}")
        if f["status"] not in TERMINAL:
            outstanding.append(fam)
        elif f["status"] in ("no_hits", "abandoned"):
            failures.append(f"{fam} ({f['status']}: {f.get('note')})")
    print("\n" + "=" * 78)
    print("RETRIEVAL STATE")
    print("=" * 78)
    print("\n".join(rows))
    print(f"\n  state file : {STATE_PATH}")
    print(f"  results    : {OUT}")

    if outstanding:
        log(f"RETRIEVAL INCOMPLETE -- still outstanding: {', '.join(outstanding)}. "
            "Rerun the same command to resume; RIDs are on disk.")
        return 2
    if failures:
        log("RETRIEVAL COMPLETE with documented failures: " + "; ".join(failures))
        return 3
    log(f"RETRIEVAL COMPLETE -- {len(families)}/{len(families)} families hold a "
        "validated result.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--families", default=None,
                    help="comma-separated subset, e.g. r2pyocin (default: all)")
    ap.add_argument("--max-minutes", type=float, default=30.0,
                    help="wall budget for this invocation (default 30)")
    ap.add_argument("--max-concurrent", type=int, default=1,
                    help="searches in flight at once (default 1: strictly serial)")
    ap.add_argument("--status", action="store_true",
                    help="print state and exit; touches no network")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    refs = json.loads((HERE / "reference_sequences.json").read_text())
    families = ([f.strip() for f in args.families.split(",")] if args.families
                else list(refs))
    for fam in families:
        if fam not in refs:
            sys.exit(f"ABORT: unknown family {fam!r}. Known: {', '.join(refs)}")

    state = load_state(families)

    def on_signal(signum, _frame):
        global _stop
        _stop = True
        log(f"received signal {signum} -- finishing the current step and saving state")
    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, on_signal)

    if args.status:
        return report(state, families)

    log(f"=== retrieval run: {len(families)} famil"
        f"{'y' if len(families) == 1 else 'ies'} ({', '.join(families)}), "
        f"budget {args.max_minutes:g} min, max {args.max_concurrent} in flight ===")
    reconcile_from_disk(state)
    save_state(state)

    deadline = time.monotonic() + args.max_minutes * 60
    while not _stop:
        fams = {fam: state["families"][fam] for fam in families}
        if all(f["status"] in TERMINAL for f in fams.values()):
            break
        if time.monotonic() > deadline:
            log(f"wall budget of {args.max_minutes:g} min reached")
            break

        in_flight = [fam for fam, f in fams.items() if f["status"] == "waiting"]
        # submit, strictly serially, only into free slots
        for fam, f in fams.items():
            if _stop or time.monotonic() > deadline:
                break
            if f["status"] == "unsubmitted" and len(in_flight) < args.max_concurrent:
                do_submit(state, fam, refs[fam]["seq"])
                if state["families"][fam]["status"] == "waiting":
                    in_flight.append(fam)

        # poll whatever is due
        due = []
        for fam in in_flight:
            f = state["families"][fam]
            last = f.get("last_poll_at")
            wait = f["backoff_s"] if last else POLL_MIN
            if last is None or (now() - parse_iso(last)).total_seconds() >= wait:
                due.append(fam)
        for fam in due:
            if _stop or time.monotonic() > deadline:
                break
            do_poll(state, fam)
            time.sleep(2)

        if _stop:
            break
        fams = {fam: state["families"][fam] for fam in families}
        if all(f["status"] in TERMINAL for f in fams.values()):
            break
        # sleep until the next job is due, but never past the deadline
        waits = []
        for fam, f in fams.items():
            if f["status"] == "waiting" and f.get("last_poll_at"):
                elapsed = (now() - parse_iso(f["last_poll_at"])).total_seconds()
                waits.append(max(5.0, f["backoff_s"] - elapsed))
        nap = min(waits) if waits else 15.0
        nap = min(nap, max(1.0, deadline - time.monotonic()))
        if nap > 0:
            time.sleep(nap)

    save_state(state)
    if _stop:
        log("interrupted -- state saved; rerun the same command to resume")
    return report(state, families)


if __name__ == "__main__":
    sys.exit(main())
