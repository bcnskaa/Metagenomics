
"""

After mapping reads to NR database, GIs need to be extracted from subject ids

"""
def filter_gi_from_sid(bla_fn, out_fn=None):
    print("Reading from " + bla_fn)
    with open(bla_fn) as IN:
        bla = IN.read().splitlines()
    bla = [b.split("\t") for b in bla]
    
    if out_fn is None:
        out_fn = bla_fn + ".gi"
    
    n = 0
    print("Exporting results to " + out_fn)
    with open(out_fn,"w") as OUT:
        for b in bla:
            n += 1
            gi = b[1].split("|")[1]
            b[1] = gi
            OUT.write("\t".join(b) + "\n")
    
    print(str(n) + " exported.")