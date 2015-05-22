
largest_gi = 815659664 + 1
gi2taxid = [0] * (largest_gi + 1)
gi2seed = [0] * (largest_gi + 1) 
gi2taxid_map_fn = "gi_taxid_prot.dmp"
gi2seed_map_fn="/home/siukinng/db/Markers/ncbi_nr/mapping/gi2seed.map"


## Import gi2seed_map
with open(gi2seed_map_fn) as IN:
    s_vals = IN.read().splitlines()

print("Number of rows imported: " + str(len(s_vals)))
for v in s_vals:
    [gi, seed_id] = map(int, v.split("\t"))
    gi2seed[gi] = seed_id


## Import gi2taxid_map
with open(gi2taxid_map_fn) as IN:
    vals = IN.read().splitlines()


print("Number of rows imported: " + str(len(vals)))
for v in vals:
    [gi, taxid] = map(int, v.split("\t"))
    gi2taxid[gi] = taxid


import glob

m8_fns = glob.glob("/home/siukinng/MG/m8/scaffolds/*.m8")


for m8_fn in m8_fns:
    print("Processing " + m8_fn)
    with open(m8_fn) as IN:
        m8 = IN.read().splitlines()
    m8_map = [0] * len(m8)


    #
    skipped_n = 0
    processed_n = 0
    for i, m in enumerate(m8):
        m = m.split("\t")
        # Query ID
        qid = m[0]
        # GI
        gi = int(m[1].split("|")[1])
        tax_id = ""
        seed_id = ""
        try:
            # tax_id
            tax_id = gi2taxid[gi]
            # seed_id
            seed_id = gi2seed[gi]
            processed_n += 1
        except:
            skipped_n += 1
        m8_map[i] = [qid, gi, tax_id, seed_id]

    out_fn = m8_fn + ".map"
    print("Exporting results to " + out_fn )
    with open(out_fn , "w") as OUT:
        OUT.write("#read_id\tgi\ttax_id\tseed_id\n")    
        for m in m8_map:
            OUT.write("\t".join(map(str, m)) + "\n")
    print("Number of processed: " + str(processed_n) + " (" + str(skipped_n) + " skipped)")

