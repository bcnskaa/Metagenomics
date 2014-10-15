import os


"""
    Upload to MG-Rast
    
    webkey=EWbN2FHyJ5M3zrxgb2FDBNRSk
    fn=P1_contigs.fa
    curl -H "auth: $webkey"  -X POST -F "upload=@./$fn" "http://api.metagenomics.anl.gov/1/inbox/" > curl_output.txt
    
"""

def retrieve(argv): 
    project_id = argv[0]
    
    cwd = os.getcwd()
    os.mkdir(project_id)
    os.chdir(project_id)
    
    # Obtain id
    with open("../" + str(project_id) +".lst") as IN:
        meta_ids = IN.read().splitlines()
        
    meta_ids = [l for l in meta_ids if len(l) > 1]
    
    print("Meta_ids=" + str(len(meta_ids)))
    
    for meta_id in meta_ids:
        print("Processing " + meta_id)
        
        cmd = "rm .listing"
        os.system(cmd)
        
        cmd = "wget --no-remove-listing ftp://ftp.metagenomics.anl.gov/projects/" + str(project_id) + "/" + meta_id + "/raw/"
        print(cmd)
        os.system(cmd)
       
        with open(".listing") as IN:
            lines = IN.read().splitlines()
        
        fn = [l.split()[8] for l in lines if ".fna.gz" in l][0]
        
        print("Retrieving " + fn)
        
        cmd = "wget ftp://ftp.metagenomics.anl.gov/projects/" + str(project_id) + "/" + meta_id + "/raw/" + fn
        print(cmd)
        os.system(cmd)
        
    os.chdir(cwd)




"""
 Given an abundance summary produced by maxbin, abundance information
 of sequence header of contig file (exported from idba_ud) can be included.
"""
def prepare_seq_cov(maxbin_abund_fn, contig_fn, output_fn):
    contig_fn = "../../idba_ud/contig.fa"
    maxbin_abund_fn = "GZ_.abund"
    
    with open(maxbin_abund_fn) as IN:
        lines = IN.read().splitlines()     
    abunds = {v.split("\t")[0]:v.split("\t")[1] for v in lines}
    
    mod_seqs = []
    for seq in SeqIO.parse(contig_fn, "fasta"):    
        if seq.id in abund.keys():
            a = str(int(round(float(abunds[seq.id]))))
            seq.id = seq.id + "_[cov=" + a + "]"
            seq.description = ""
            mod_seqs.append(seq)
    
    OUT = open(output_fn, "w")
    SeqIO.write(mod_seqs, OUT, "fasta")
    OUT.close()




# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    