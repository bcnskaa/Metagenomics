import os


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



# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    