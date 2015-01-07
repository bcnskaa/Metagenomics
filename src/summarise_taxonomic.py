"""
 Steps:
     1. run_blastx: bash run_blastx.sh
     2. map_bla2specI: 
     3. make_summary:
     4. summarise:
     5. pick_otu:
     
     
     import summarise_taxonomic
     summarise_taxonomic.map_bla2specI()
     summarise_taxonomic.make_summary()
     summarise_taxonomic.summarise()
     all_tax = summarise_taxonomic.pick_otu()
     
"""
from __future__ import print_function
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
        #summary_list[bin_id].append(sorted_species[0])
        summary_list[bin_id].insert(0,sorted_species[0])

    bin_ids = summary_list.keys()
    bin_ids = sorted(bin_ids)
    
    with open(summary_ofn, "w") as OUT:
        for bin_id in bin_ids:
            bin_summary = summary_list[bin_id]
            summary = [str(g[0]) + ":" +  str(g[1]) for g in bin_summary]
            OUT.write(bin_id + "\t" + "\t".join(summary) + "\n")



def get_tax_val(tax, tax_level="s__"):
    tax_delim = ";"
    tax_val = ""
    val_sidx = tax.find(tax_level)
    if val_sidx != -1:
        val_eidx = tax.find(tax_delim, val_sidx)
        if val_eidx == -1:
            val_eidx = len(tax)
        tax_val = tax[val_sidx : val_eidx]
        tax_val = tax_val.replace(tax_level, "")
    #print(tax_val)
    return tax_val



"""
gg_otu_fn="/home/siukinng/db/Taxanomy/otu_id_to_greengenes.txt"
with open(gg_otu_fn) as IN:
    otu = IN.read().splitlines()
o = otu[12121] 
summarise_taxonomic.get_tax_val(o, "g__")
otu_s = {summarise_taxonomic.get_tax_val(o, tax_level="s__"): o.split("\t")[1] for o in otu if len(summarise_taxonomic.get_tax_val(o, tax_level="s__")) > 0}

otu_g = {summarise_taxonomic.get_tax_val(o, tax_level="g__"): o.split("\t")[1] for o in otu if len(summarise_taxonomic.get_tax_val(o, tax_level="g__")) > 0}


"""
def pick_otu(summary_ifn="all.tax", summary_ofn="all.updated.tax", gg_otu_fn="/home/siukinng/db/Taxanomy/otu_id_to_greengenes.txt"):
    with open(summary_ifn) as IN:
        all_tax = IN.read().splitlines()
    all_tax = [t.split("\t") for t in all_tax]
    
    # Import Greengenes
    with open(gg_otu_fn) as IN:
        otu = IN.read().splitlines()
    otu_s = {get_tax_val(o, tax_level="s__"): o.split("\t")[1] for o in otu if len(get_tax_val(o, tax_level="s__")) > 0}
    #otu = {(o.split(";")[len(o.split(";")) - 1]).replace("s__",""): o.split("\t")[1] for o in otu if len(o.split(";")[len(o.split(";")) - 1]) > 3}
    #otu_g = {(o.split(";")[len(o.split(";")) - 2]).replace("g__",""): o.split("\t")[1] for o in otu if len(o.split(";")[len(o.split(";")) - 2]) > 3}
    otu_g = {get_tax_val(o, tax_level="g__"): o.split("\t")[1] for o in otu if len(get_tax_val(o, tax_level="g__")) > 0}
    
    for i in range(0, len(all_tax)):
        tax = all_tax[i]
        max_idx = len(tax) - 1
        confidence = float(tax[1].split(":")[1]) / float(tax[2].split(":")[1])
        #confidence = float(tax[max_idx].split(":")[1]) / float(tax[1].split(":")[1])
        sp = " ".join(tax[1].split(" ")[0:2])
        #sp = " ".join(tax[max_idx].split(" ")[0:2])
        
        o = ""
        if sp in otu_s.keys():
            o = otu_s[sp]
            
        if len(o) == 0:
            g = tax[1].split(" ")[0]
            if g in otu_g.keys():
                o = otu_g[g]
        
        all_tax[i].insert(1, o)    
        all_tax[i].insert(1, str(confidence))
        
    with open(summary_ofn, "w") as OUT:
        for tax in all_tax:
            OUT.write("\t".join(tax) + "\n")

    return all_tax


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
# all.updated.tax.finalized is a finalized version of all.updated.tax

"""
def apply_abundance(summary_ifn="all.updated.tax.finalized", out_fn=None, normalized_abund=True, abund_dir="/disk/rdisk08/siukinng/MG/scaffolds_5000/", dir_suffix="_5000", bin_summary_fn_suffix=".summary"):
    import glob
    import os 
    
    with open(summary_ifn) as IN:
        all_tax = IN.read().splitlines()
    all_tax = [t.split("\t") for t in all_tax if len(t.split("\t")) > 1] 
    
    abund = {}
    abund_fns = glob.glob(abund_dir + "*" + dir_suffix + "/*" + bin_summary_fn_suffix)
    for abund_fn in abund_fns:
        sample_id = (os.path.basename(abund_fn)).replace(bin_summary_fn_suffix, "")
        print(sample_id)
        with open(abund_fn) as IN:
            aa = IN.read().splitlines()
        del aa[0]
        
        if normalized_abund:
            aa = {(a.split("\t")[0]).replace(".fasta",""): float(a.split("\t")[1]) for a in aa}
            total_aa = sum([aa[a] for a in aa.keys()])
            print(list(aa.keys())[0])
            aa = {a:str(aa[a] / total_aa) for a in aa.keys()}
            abund.update(aa)
        else:
            abund.update({(a.split("\t")[0]).replace(".fasta",""): a.split("\t")[1] for a in aa})
    
    if out_fn is None:
        out_fn = summary_ifn + ".abund"
    
    with open(out_fn, "w") as OUT:
        for tax in all_tax:
            a = abund[tax[0]]
            g = (tax[2]).split("f__")[1]
            g = g.split(";")[0]
            OUT.write(tax[0] + "\t" + g + "\t" + a + "\t" + "\t".join(tax[1:]) + "\n")

        


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

"""
# Generate stacked bar plot
tax_fn = "all.updated.tax.finalized.f.abund"
tax <- read.table("all.updated.tax.finalized.f.abund", sep="\t", header=F, stringsAsFactors=F)
tax$group <- sapply(1:length(tax$V1), function(i) { strsplit(tax$V1[i],".", fixed=T)[[1]][1]})

tax <- tax[which(tax$group != "SWH-Cell55_Y2"), ]
tax <- tax[which(tax$group != "SWH-Seed_Y0"), ]
tax <- tax[which(tax$group != "GZ-Seed_Y0"), ]

tax <- tax[,c("group", "V2", "V3")]
tax <- aggregate(tax["V3"], by=tax[c("group","V2")], FUN=sum)

library(lattice)
#ggplot(tax, aes(x=group, y=V3, fill=V2)) + geom_bar(stat="identity", aes(group=tax$group)) + coord_polar(theta="y")
pdf(paste(tax_fn, ".spectral.pdf" , sep=""))
ggplot(tax, aes(x=group, y=V3, fill=V2)) + geom_bar(stat="identity", aes(group=tax$group)) + scale_fill_brewer(palette="Spectral") + theme(axis.text.x = element_text(angle = -45, hjust = 0)) + ylab("Relative Abundance") + xlab("Enrichment Group")
dev.off()
pdf(paste(tax_fn, ".all.pdf" , sep=""))
ggplot(tax, aes(x=group, y=V3, fill=V2)) + geom_bar(stat="identity", aes(group=tax$group)) + theme(axis.text.x = element_text(angle = -45, hjust = 0)) + ylab("Relative Abundance") + xlab("Enrichment Group")
dev.off()

#ggplot(tax, aes(y=rep(1, nrow(tax)), x=V3, fill=V2)) + geom_bar(stat="identity", aes(group=tax$group)) + coord_polar(theta="y")

"""
