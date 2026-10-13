import json,sys,urllib.request
def search(v,rows=100):
    body={"query":{"type":"terminal","service":"full_text","parameters":{"value":v}},
          "return_type":"entry","request_options":{"paginate":{"start":0,"rows":rows},"results_verbosity":"compact"}}
    r=urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",data=json.dumps(body).encode(),headers={"Content-Type":"application/json"})
    try: return json.load(urllib.request.urlopen(r)).get("result_set",[])
    except: return []
KW=['sheath','trunk','contract','extend','tail','baseplate','tube']
for term in sys.argv[1:]:
    print("##### %s"%term)
    for i in search(term):
        try: d=json.load(urllib.request.urlopen("https://data.rcsb.org/rest/v1/core/entry/%s"%i))
        except: continue
        ti=d['struct'].get('title','')
        if not any(k in ti.lower() for k in KW): continue
        ci=d.get('rcsb_primary_citation',{})
        print("  %s %s %s | %s | %s %s PM:%s"%(i,d['rcsb_entry_info'].get('resolution_combined'),d['rcsb_accession_info'].get('initial_release_date','')[:10],ti[:115],ci.get('journal_abbrev'),ci.get('year'),ci.get('pdbx_database_id_PubMed')))
