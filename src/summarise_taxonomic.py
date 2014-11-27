"""
 Steps:
     1. run_blastx
     2. map_bla2specI
     3. make_summary
     4. summarise
     5. pick_otu
"""
import os
import operator
import glob


def summarise(summary_ofn="all.tax", in_dir=".", summary_fn_ext="summary"):
    #summary_ofn = "all.tax"
    summary_fns = glob.glob(in_dir + "/*." + summary_fn_ext)
    
    summary_list = {}
    for summary_fn in summary_fns:
        bin_id = (os.path.basename(summary_fn)).replace("." + summary_fn_ext, "")
        print("Processing " + bin_id)   
        with open(summary_fn) as IN:
            summary = IN.read().splitlines()      
        genus = [s.split(" ")[0] for s in summary]        
        genus_ids = list(set(genus))
        genus_summary = {id: genus.count(id) for id in genus_ids}     
        sorted_genus = sorted(genus_summary.items(), key=operator.itemgetter(1))
        sorted_genus.reverse() 
        
        species = summary
        species_ids = list(set(species))
        species_summary = {id: species.count(id) for id in species_ids}
        sorted_species = sorted(species_summary.items(), key=operator.itemgetter(1))
        sorted_species.reverse()
        print("Dominant species=")
        print(sorted_species[0])
        
        summary_list[bin_id] = sorted_genus[0:10]
        summary_list[bin_id].append(sorted_species[0])

    bin_ids = summary_list.keys()
    bin_ids = sorted(bin_ids)
    
    with open(summary_ofn, "w") as OUT:
        for bin_id in bin_ids:
            bin_summary = summary_list[bin_id]
            summary = [str(g[0]) + ":" +  str(g[1]) for g in bin_summary]
            OUT.write(bin_id + "\t" + "\t".join(summary) + "\n")



"""
"""
def pick_otu(gg_out="/home/siukinng/db/Taxanomy/otu_id_to_greengenes.txt"):
    print("ok")




"""
 
"""
def run_blastx(in_dir="../../DistMtx/bins/specI", out_dir="blast", specI_dir="/home/siukinng/db/Markers/specI/sequences", fasta_fn_ext="faa"):
    faa_fns = glob.glob(in_dir + "/*." + fasta_fn_ext)
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    print("Number of fasta files to be processed=" + str(len(faa_fns)))

    for faa_fn in faa_fns:
        hmm_id = (os.path.basename(faa_fn)).replace(fasta_fn_ext, "")
        cmd = "~/tools/blast/bin/tblastn -query " + faa_fn + " -db " + specI_dir + "/" + hmm_id + ".fna -outfmt 6 -out " + out_dir + "/" + hmm_id + ".bla -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 3 -num_threads 16"
        print(cmd)
        os.system(cmd)
        
  
"""
 
"""     
def map_bla2specI(path=".", mapped_bla_fn_ext="mapped", map_fn="/home/siukinng/db/Markers/specI/data/RefOrganismDB.v9.taxids2names"):
    fns = glob.glob(path + "/*.bla")
    
    for fn in fns:
        print("Processing " + fn)
        with open(fn) as IN:
            bla = IN.read().splitlines()
            
        with open(map_fn) as IN:
            map = IN.read().splitlines()
        map = {m.split("\t")[0]: m.split("\t")[1] for m in map}
        
        out_fn = fn + ".mapped"
        with open(out_fn, "w") as OUT:
            for b in bla:
                b = b.split("\t")
                id = b[1].split(".")[0]
                if id in map.keys():
                    b[1] = b[1] + ":" + map[id]
                else:
                    b[1] = b[1] + ":" + "Unknown"
                OUT.write("\t".join(b) + "\n")



"""
 only works for bin-group 
""" 
def make_summary(in_dir=".", fasta_fn_ext="fasta", mapped_bla_fn_ext="mapped", summary_fn_ext="summary"):
    mapped_fns = glob.glob(in_dir + "/*." + mapped_bla_fn_ext)
    
    mapped_list = {}
    for mapped_fn in mapped_fns:
        with open(mapped_fn) as IN:
            bla = IN.read().splitlines()
        bla = [b.split("\t") for b in bla]
        for b in bla:
            qid = (b[0][::-1].split("_",2)[2])[::-1]
            if qid not in mapped_list.keys():
                mapped_list[qid] = []
            sid = b[1].split(":")[1]
            mapped_list[qid].append(sid)
            
    for id in mapped_list.keys():
        with open(id + "." + summary_fn_ext, "w") as OUT:
            sids = mapped_list[id]
            OUT.write("\n".join(sids))

        
#     bla_fns = glob.glob(in_dir + "/*." + mapped_bla_fn_ext)
#     for bla_fn in bla_fns:
#         id = ((os.path.basename(bla_fn)).split(mapped_bla_fn_ext)[0]).replace(".bla", "")
#         cmd = "cat " + bla_fn + " | grep -e \"^" + id + "\" | cut -f2 | cut -d':' -f2 >> "  + bla_fn + ".summary"
#         print(cmd)
#         os.system(cmd)



"""
# run_blastx.sh 
# Generate command for running blastx
specI_dir="~/db/Markers/specI/sequences"
out_dir="blast"
for f in  ../../DistMtx/bins/specI/*.faa;do
    hmm_id=${f##*/}
    hmm_id=${hmm_id/.faa/}

    #echo $hmm_id

    cmd="~/tools/blast/bin/tblastn -query $f -db $specI_dir/$hmm_id.fna -outfmt 6 -out $out_dir/$hmm_id.bla -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -ma
x_target_seqs 3 -num_threads 16"
    echo $cmd

done



# summary.sh
for f in ../../../*_5000/*.fasta;do
    id=${f##*/};
    id=${id/.fasta/};
    echo $id;
    rm $id.summary
    for bla_fn in *.bla.mapped;do
        cat $bla_fn | grep -e "^$id" | cut -f2 | cut -d':' -f2 >> $id.summary
    done
done
"""



"""

python extract_hmm_seq.py


# Combined individual bin sequences
hmm_ids=("specI" "dbCAN" "TIGRFAM")

for hmm_id in ${hmm_ids[@]};do
    for d in *_5000;do
            sample_id=${d##*/}
            sample_id=${sample_id/_5000/}
    
            combined_seq_dir=$d/DistMtx/bins/$sample_id/$hmm_id
            mkdir -p $combined_seq_dir
    
            rm $combined_seq_dir/*.faa
    
            for bin_id in $d/DistMtx/bins/$sample_id.*;do
                    echo "Processing $bin_id"
                    fn_n=`ls $bin_id/$hmm_id/*.faa | wc -l`
    
                    if [ "$fn_n" -gt "0" ];then
                            for hmm_fn in $bin_id/$hmm_id/*.faa;do
                                    touch $combined_seq_dir/${hmm_fn##*/}
                                    cat $hmm_fn >> $combined_seq_dir/${hmm_fn##*/}
                            done
                    fi
            done
    done
done



"""
