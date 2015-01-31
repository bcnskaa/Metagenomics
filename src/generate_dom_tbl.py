import glob
#import process_HMM
import os
import mg_pipeline


def generate_dom_tbl(out_fn, hmm_id = "dbCAN", is_bin=False):
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



def append_abund_to_dom_tbl(tbl_fn, out_fn, list_ofn=None, summary_dir="."):
    print("Processing " + tbl_fn)
    
    with open(tbl_fn) as IN:
        lines = IN.read().splitlines()
    header = lines[0]
    header = header.split("\t")[2:]
    del lines[0]
    
    # Abundance information
    nr_headers = list(set([h.split(".")[0] for h in header]))
    group_summary = {}
    for nr_header in nr_headers:
        summary_fn = "/disk/rdisk08/siukinng/MG/scaffolds_5000/" + nr_header + "_5000/"+ nr_header+ ".summary"
        with open(summary_fn) as IN:
            summary = IN.read().splitlines()
        del summary[0]
        summary = {(s.split("\t")[0]).replace(".fasta",""):float(s.split("\t")[1]) for s in summary}
        group_summary.update(summary)
        
    
    # Parse the dom table
    line_n = len(lines)
    desc = []
    row_ids = []
    dom_table = {h:[0 for i in range(line_n)] for h in header}
    dom_mtx = [["HEADER","DESC"]+header]
    dom_list = []
    for i, l in enumerate(lines):
        l = l.split("\t")
        row_id = l[0]
        row_ids.append(row_id)
        desc.append(l[1])
        row_items = [l[0], l[1]]
        l = l[2:]
        for j, id in enumerate(header):
            count = l[j]
            abund = group_summary[id]
            dom_table[id][i] = abund * int(count)
            row_items.append(str(abund * int(count)))
            dom_list.append([id, row_id, str(count), str(abund * int(count))])
        dom_mtx.append(row_items)
        
    with open(out_fn, "w") as OUT:
        for row in dom_mtx:
            OUT.write("\t".join(row) + "\n")
         
    if list_ofn is not None:
        with open(list_ofn, "w") as OUT:
            for dom in dom_list:
                OUT.write("\t".join(dom) + "\n")
                   
        
    return dom_mtx
    

"""
import generate_dom_tbl
tbl_fn = "dbCAN.bin.abund.tbl"
generate_dom_tbl.consolidate_dom_tbl_by_taxon(tbl_fn, melt_tbl=True)

tbl_fn = "dbCAN.bin.tbl"
generate_dom_tbl.consolidate_dom_tbl_by_taxon(tbl_fn, melt_tbl=True)

"""
def consolidate_dom_tbl_by_taxon(tbl_fn, consolidated_tbl_ofn=None, melt_tbl=False, melt_consolidated_tbl_ofn=None, taxon_fn="/home/siukinng/MG/scaffolds_5000/TaxonomicAssignment/bins/blast/all.updated.tax.finalized.f.abund"):
    print("Processing " + tbl_fn)
    
    with open(tbl_fn) as IN:
        lines = IN.read().splitlines()
    header = lines[0]
    header = header.split("\t")[2:]
    del lines[0]
    row_id = [l.split("\t")[0] for l in lines]
    row_desc = [l.split("\t")[1] for l in lines]
    tbl_data = [l.split("\t")[2:] for l in lines]
    
    # taxonomical information
    nr_headers = list(set([h.split(".")[0] for h in header]))
    
    taxon = {}
    with open(taxon_fn) as IN:
        taxon = IN.read().splitlines()
    taxon = {t.split("\t")[0]: t.split("\t")[1] for t in taxon}
    
    list_taxon = {}
    group_taxon = {}
    group_taxon_long = {}
    # Generate a mapping table
    for tid in taxon.keys():
        group_id = tid.split(".")[0]
        taxon_id = taxon[tid]
        if group_id not in group_taxon.keys():
            group_taxon[group_id] = {}
            group_taxon_long[group_id] = {}
        if taxon_id not in group_taxon[group_id].keys():
            group_taxon[group_id][taxon_id] = []
            group_taxon_long[group_id][taxon_id] = []
        #group_taxon[group_id][taxon_id].append(tid)
        group_taxon_long[group_id][taxon_id].append(tid)
        group_taxon[group_id][taxon_id].append(header.index(tid))
    
    # Compute the number of columns remained after consolidation
    col_n = sum([len(group_taxon[gid]) for gid in group_taxon.keys()])    
    
    # Generate an empty matrix for storing the merge data
    consolidated_tbl_data = [[0.0 for j in range(col_n)] for i in range(len(tbl_data))]
    
    # Parse the tbl data
    for tbl_idx, tbl in enumerate(tbl_data):
        cur_consolidated_tbl_data_pos = 0
        # Select group
        for gid in group_taxon.keys():
            for taxon_id in group_taxon[gid].keys():
                group_sum = 0.0
                for col_idx in group_taxon[gid][taxon_id]:
                    group_sum = group_sum + float(tbl[col_idx])
                #consolidated_tbl_data[tbl_idx][cur_consolidated_tbl_data_pos] = str(group_sum)
                consolidated_tbl_data[tbl_idx][cur_consolidated_tbl_data_pos] = group_sum
                cur_consolidated_tbl_data_pos = cur_consolidated_tbl_data_pos + 1
                
    if consolidated_tbl_ofn is None:
        consolidated_tbl_ofn = tbl_fn + ".consolidated"
    
    print("Exporting results to " + consolidated_tbl_ofn)
    
    # Export the result
    with open(consolidated_tbl_ofn, "w") as OUT:
        # Group header:
        group_header = ["HEADER", "DESC"]     
        # Taxon header:
        taxon_header = ["HEADER", "DESC"]
        for gid in group_taxon.keys():
            for tid in group_taxon[gid].keys():
                group_header.append(gid)
                taxon_header.append(tid) 
        
        OUT.write("\t".join(group_header) + "\n")
        OUT.write("\t".join(taxon_header) + "\n")
        for i, tbl in enumerate(consolidated_tbl_data):
            OUT.write(row_id[i] + "\t" + row_desc[i] + "\t" + "\t".join([str(t) for t in tbl]) + "\n")   
    
    # Melt form: group_id\ttaxon_id\trow_id\tvalue   
    if melt_tbl:
        if melt_consolidated_tbl_ofn is None:
            melt_consolidated_tbl_ofn = consolidated_tbl_ofn + ".melt"

        print("Melted consolidated dom tbl is exporting to " + melt_consolidated_tbl_ofn)


        with open(melt_consolidated_tbl_ofn, "w") as OUT:
            for k, row_id in enumerate(row_id):
                col_idx = 0
                for i, gid in enumerate(group_taxon.keys()):    
                    for j, tid in enumerate(group_taxon[gid].keys()):
                        OUT.write("\t".join([gid, tid, row_id, str(consolidated_tbl_data[k][col_idx])]) + "\n")
                        col_idx = col_idx + 1


"""
#Merge two tbl files for ease of post-processings
import generate_dom_tbl
tbl_1_fn = "dbCAN.bin.abund.tbl.consolidated.melt"
tbl_2_fn = "dbCAN.bin.tbl.consolidated.melt"
tbl_merged_fn = "dbCAN.bin.tbl.consolidated.melt.abund_merged"
generate_dom_tbl.merge_melt_tbl(tbl_1_fn, tbl_2_fn, tbl_merged_fn)

    
"""
def merge_melt_tbl(tbl_1_fn, tbl_2_fn, out_fn):
    print("Merging " + tbl_1_fn + " and " + tbl_2_fn)
    
    with open(tbl_1_fn) as IN:
        tbl_1 = IN.read().splitlines()
    tbl_1 = {"\t".join(t.split("\t")[0:3]):t.split("\t")[3] for t in tbl_1}
    
    with open(tbl_2_fn) as IN:
        tbl_2 = IN.read().splitlines()
    #tbl_2 = [t.split("\t") for t in tbl_2]
    tbl_2 = {"\t".join(t.split("\t")[0:3]):t.split("\t")[3] for t in tbl_2}
    
    
    with open(out_fn, "w") as OUT:
        for tid in tbl_1.keys():
            tbl_2_val = "-"
            if tid in tbl_2.keys():
                tbl_2_val = tbl_2[tid]
            tbl = [tid, tbl_1[tid], tbl_2_val]
            
            OUT.write("\t".join(tbl) + "\n")

    

def count_gene(out_fn=None, fasta_fn_dir=".", fasta_fn_ext="renamed_faa"):
    from Bio import SeqIO
    import glob
    import os
    
    fasta_fns = glob.glob(fasta_fn_dir + "/*."+fasta_fn_ext)
    if len(fasta_fns) == 0:
        print("No file can be read.")
        return None
    
    row_ids_fns = {(os.path.basename(f)).replace(fasta_fn_ext, ""):f for f in fasta_fns}
    row_ids = list(row_ids_fns.keys())
    
    bin_ids = {}
    # Construct a list of all bin ids
    for fasta_fn in fasta_fns:
        seqs = SeqIO.index(fasta_fn, "fasta")    
        seq_ids = [(id[::-1].split("_",1)[1])[::-1] for id in list(seqs.keys())]
        for id in seq_ids:
            if id not in bin_ids.keys():
                bin_ids[id] = [0 for i in range(len(row_ids))]
        
    for i, row_id in enumerate(row_ids):
        seqs = SeqIO.index(row_ids_fns[row_id], "fasta")
        seq_ids = list(seqs.keys())
        for seq_id in seq_ids:
            seq_id = (seq_id[::-1].split("_",1)[1])[::-1]
            bin_ids[seq_id][i] = bin_ids[seq_id][i] + 1
        
    if out_fn is None:
        out_fn = "count_genes"
       
    with open(out_fn, "w") as OUT:
        bids = list(bin_ids.keys())
        OUT.write("ID\t" + "\t".join(bids) + "\n")
        for i, row_id in enumerate(row_ids):
            OUT.write(row_id)
            for bid in bids:
                OUT.write("\t" + str(bin_ids[bid][i]))   
            OUT.write("\n") 
    
    
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


"""
# Generate scaffold2bin_id table
for d in *_5000;do
    id=${d/_5000/}
    out_fn="$d/$id.scaffold2bin_id"
    rm $out_fn
    for f in $d/*.fasta;do
        bin_id=${f##*/}
        bin_id=${bin_id/.fasta/}
        cat $f | grep -e "^>" | sed -e "s/>//" | sed -e "s/$/\t$bin_id/" >> $out_fn
    done
done

"""