



"""

len_fns = ["SWH-Seed_Y0.GZ-Seed_Y0", "SWH-Xyl_Y1.SWH-Seed_Y0", "SWH-Seed_Y0.SWH-Cell_Y1", "GZ-Seed_Y0.GZ-Cell_Y1", "GZ-Xyl_Y1.GZ-Seed_Y0"]

len_fns = ["SWH-Cell_Y2.SWH-Cell_Y1", "SWH-Xyl_Y2.SWH-Xyl_Y1", "GZ-Cell_Y2.GZ-Cell_Y1", "GZ-Xyl_Y2.GZ-Xyl_Y1"]

len_fns = ["SWH-Cell_Y2.GZ-Cell_Y2", "SWH-Xyl_Y2.GZ-Xyl_Y2", "SWH-Cell_Y1.GZ-Cell_Y1", "SWH-Xyl_Y1.GZ-Xyl_Y1"]

len_fns = ["SWH-Xyl_Y2.SWH-Cell_Y2", "SWH-Xyl_Y1.SWH-Cell_Y1", "GZ-Xyl_Y2.GZ-Cell_Y2", "GZ-Xyl_Y1.GZ-Cell_Y1"]

for len_fn in len_fns:
    print(len_fn + "\t" + str(sum_len(len_fn + ".bla.b1000.0.bla2.len")))

    
"""


def sum_len(len_fn):
    with open(len_fn) as IN:
        len = IN.read().splitlines()
    len = [int(l) for l in len]
    return sum(len)