import glob
#import process_HMM
import os
import mg_pipeline


def generate_dom_tbl(out_fn, hmm_id = "dbCAN", is_bin=False, ):
    #hmm_id = "dbCAN"
    if is_bin:
        fns = glob.glob("*_5000/Markers/bins/" + hmm_id + "/*.dom.tbl")
    else:
        fns = glob.glob("*_5000/Markers/" + hmm_id + "/*.dom.tbl")

    ids = []
    hmm_score_threshold=200
    hmm_tc_fn = "/home/siukinng/db/Markers/" + hmm_id + "/" + hmm_id + ".tc"
    
    if not os.path.isfile(hmm_tc_fn):
            hmm_tc_fn = None
            
    print("Processing " + hmm_id + " and export the results to " + out_fn)
            
    dom_tbls = {}
    for fn in fns:
       #hmm_orf_dict = process_HMM.postprocess_HMMER_search(fn, hmm_score_threshold=hmm_score_threshold, hmm_tc_fn=hmm_tc_fn)
       hmm_orf_dict = mg_pipeline.postprocess_HMMER_search_by_fn(fn, hmm_score_threshold=hmm_score_threshold, hmm_tc_fn=hmm_tc_fn)
       
       id = (os.path.basename(fn)).replace(".dom.tbl", "")
       id = id.replace("."+hmm_id, "")
       ids.append(id)
       #dom_tbl = process_HMM.generate_dom_tbl(hmm_orf_dict)
       dom_tbl = mg_pipeline.generate_dom_tbl(hmm_orf_dict)
       
       dom_tbls[id] = dom_tbl
       
    
    dom_ids = []
    for k in dom_tbls.keys():
        dom_ids.extend(dom_tbls[k].keys())
    dom_ids = list(set(dom_ids))
    
    
    dom_list = {}
    dom_list["header"] = [0 for i in range(0, len(fns))]
    for dom_id in dom_ids:
        dom_list[dom_id] = [0 for i in range(0, len(fns))]
        
    for idx, id in enumerate(ids):
        dom_list["header"][idx] = id
        print("Processing " + id)
        for dom_id in dom_ids:
            if dom_id in dom_tbls[id].keys():
                #print("Dom=" + dom_id)
                seq_ids = [d[0] for d in dom_tbls[id][dom_id]]
                seq_ids = list(set(seq_ids))
                #dom_list[dom_id][idx] = len(dom_tbls[id][dom_id])
                dom_list[dom_id][idx] = len(seq_ids)
    
    with open("/home/siukinng/db/Markers/" + hmm_id + "/" + hmm_id + ".hmm.info") as IN:
        hmm_info = IN.read().splitlines()
    hmm_info = {v.split("\t")[0]:v.split("\t")[1] for v in hmm_info}
        
    # if hmm_id == "specI":
    #     with open("/home/siukinng/db/Markers/specI/specI.hmm.info") as IN:
    #         specI_info = IN.read().splitlines()
    #     specI_info = {v.split("\t")[0]:v.split("\t")[1] for v in specI_info}
    
    with open(hmm_id + ".tbl", "w") as OUT:
        OUT.write("HEADER\tDESC\t" + "\t".join(dom_list["header"]) + "\n")
        for dom_id in dom_ids:
            if dom_id in hmm_info.keys():
                info = hmm_info[dom_id]
            else:
                info = dom_id
            OUT.write(dom_id + "\t" + info + "\t" + "\t".join(str(v) for v in dom_list[dom_id]) + "\n") 



def append_abund_to_dom_tbl(tbl_fn, out_fn, summary_dir="."):
    print("Processing " + tbl_fn)
    
    with open(tbl_fn) as IN:
        lines = IN.read().splitlines()
    header = lines[0]
    header = header.split("\t")[2:]
    del lines[0]
    
    nr_headers = list(set([h.split(".")[0] for h in header]))
    group_summary = {}
    for nr_header in nr_headers:
        summary_fn = "/disk/rdisk08/siukinng/MG/scaffolds_5000/" + nr_header + "_5000/"+ nr_header+ ".summary"
        with open(summary_fn) as IN:
            summary = IN.read().splitlines()
        del summary[0]
        summary = {(s.split("\t")[0]).replace(".fasta",""):float(s.split("\t")[1]) for s in summary}
        group_summary.update(summary)
        
    
    line_n = len(lines)
    desc = []
    row_ids = []
    dom_table = {h:[0 for i in range(line_n)] for h in header}
    dom_mtx = [["HEADER","DESC"]+header]
    for i, l in enumerate(lines):
        l = l.split("\t")
        row_ids.append(l[0])
        desc.append(l[1])
        row_items = [l[0], l[1]]
        l = l[2:]
        for j, id in enumerate(header):
            count = l[j]
            abund = group_summary[id]
            dom_table[id][i] = abund * int(count)
            row_items.append(str(abund * int(count)))
        dom_mtx.append(row_items)
        
    with open(out_fn, "w") as OUT:
        for row in dom_mtx:
            OUT.write("\t".join(row) + "\n")
            
        
    return dom_mtx
    

"""
# do Prodigal
for d in *_5000;do
    id=${d/_5000/}
    for f in $d/*.fasta;do
        bin_id=${f##*/}
        bin_id=${bin_id/.fasta/}
        cmd="~/tools/Prodigal/prodigal -p meta -a $d/Prodigal/bins/$bin_id.faa -d $d/Prodigal/bins/$bin_id.fna -o $d/Prodigal/bins/$bin_id.out -i $f"
        echo $cmd
    done
done


hmm_id="dbCAN"
hmm_fn="/home/siukinng/db/Markers/"$hmm_id"/"$hmm_id".hmm"
for d in *_5000;do
    #id=${d/_5000/}
    outdir="$d/Markers/bins/$hmm_id"
    for f in $d/Prodigal/bins/*.faa;do
        id=${f##*/}
        id=${id/.faa/}
        cmd="~/tools/hmmer/bin/hmmsearch -o $outdir/$id.$hmm_id.out -A $outdir/$id.$hmm_id.faa --tblout $outdir/$id.$hmm_id.tbl --domtblout $outdir/$id.$hmm_id.dom.tbl --pfamtblout $outdir/$id.$hmm_id.pfam.tbl $hmm_fn $f"
        echo $cmd
    done
done

"""


"""
hmm_ids=("specI" "dbCAN" "TIGRFAM" "Pfam")

for f in *_5000/Prodigal/*.faa;do
    id=${f##*/}
    id=${id/.faa/}
    
    outdir=$id"_5000/Markers/combined/bins"
    
    mkdir -p $outdir
    
    #faa_fn="$d/Prodigal/$id.faa"
    faa_fn=$f
    
    for hmm_id in ${hmm_ids[*]};do
        hmm_fn="/home/siukinng/db/Markers/"$hmm_id"/"$hmm_id".hmm"
        cmd="~/tools/hmmer/bin/hmmsearch -o $outdir/$id.$hmm_id.out -A $outdir/$id.$hmm_id.faa --tblout $outdir/$id.$hmm_id.tbl --domtblout $outdir/$id.$hmm_id.dom.tbl --pfamtblout $outdir/$id.$hmm_id.pfam.tbl $hmm_fn $faa_fn"    
        echo $cmd
    done
done

"""