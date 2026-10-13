import json,sys,urllib.request,time
def search(v,rows=100):
    body={"query":{"type":"terminal","service":"full_text","parameters":{"value":v}},
          "return_type":"entry","request_options":{"paginate":{"start":0,"rows":rows},"results_verbosity":"compact"}}
    r=urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",
        data=json.dumps(body).encode(),headers={"Content-Type":"application/json"})
    try:
        d=json.load(urllib.request.urlopen(r))
    except Exception as e:
        return []
    return d.get("result_set",[])
def meta(i):
    try:
        d=json.load(urllib.request.urlopen("https://data.rcsb.org/rest/v1/core/entry/%s"%i))
    except Exception as e:
        return None
    ci=d.get('rcsb_primary_citation',{})
    return dict(id=i,title=d['struct'].get('title'),
        res=d.get('rcsb_entry_info',{}).get('resolution_combined'),
        rel=d.get('rcsb_accession_info',{}).get('initial_release_date','')[:10],
        ne=d['rcsb_entry_info']['polymer_entity_count'],
        cite="%s | %s %s PM:%s"%(ci.get('title'),ci.get('journal_abbrev'),ci.get('year'),ci.get('pdbx_database_id_PubMed')))
if __name__=="__main__":
    seen={}
    for term in sys.argv[1:]:
        ids=search(term)
        print("##### %s -> %d"%(term,len(ids)))
        for i in ids:
            if i in seen: continue
            seen[i]=1
            m=meta(i)
            if m: print("  %s [%s %s] %s\n      %s"%(m['id'],m['res'],m['rel'],m['title'][:150],m['cite'][:140]))
