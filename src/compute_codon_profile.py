from __future__ import division
from Bio import SeqIO
import calculate_tetranucleotide_freq
import glob
import os
import operator

"""
import compute_codon_profile

codon_tbls = compute_codon_profile.process_all()
compute_codon_profile.export_codon_tbl(codon_tbls, "clostridium.codons")

"""

def generate_codon_table():
    codons = calculate_tetranucleotide_freq.get_patterns(kmer_len=3)
    codon_tbl = {c: 0 for c in codons}
    return codon_tbl


def compute_codon_profile(seq_fn):
    n = 3
    codon_error_n = 0
    print("Extracting the codon profile from " + seq_fn)
    seqs = SeqIO.index(seq_fn, "fasta")
    seq_ids = list(seqs.keys())
    codon_table = generate_codon_table()
    for seq_id in seq_ids:
        seq = str(seqs[seq_id].seq)
        codons = [seq[i : i + n] for i in range(0, len(seq), n)]
        for codon in codons:
            if codon not in codon_table.keys():
                codon_error_n = codon_error_n + 1
                continue
            codon_table[codon] = codon_table[codon] + 1
        
    if codon_error_n > 0:
        print("Warning: " + str(codon_error_n) + " malformed/truncated codon(s) found\n")
    return codon_table


def process_all(dir="./*"):
    ffn_dirs = glob.glob(dir)
    
    codon_tbls = {}
    for ffn_dir in ffn_dirs:
        if not os.path.isdir(ffn_dir):
            continue
        ffns = glob.glob(ffn_dir + "/*.ffn")
        if len(ffns) == 0:
            print(ffn_dir + " skpped")
            continue
        species_id = os.path.basename(ffn_dir)
        print("Processing " + species_id)
        
        ffn_size = {f:(os.stat(f)).st_size for f in ffns}
        sorted_ffn_size = sorted(ffn_size.items(), key=operator.itemgetter(1), reverse=True)
        codon_tbls[species_id] = compute_codon_profile(sorted_ffn_size[0][0])
        
    return codon_tbls



def export_codon_tbl(codon_tbls, out_fn):
    codons = calculate_tetranucleotide_freq.get_patterns(kmer_len=3)
    species_ids = list(codon_tbls.keys())
    with open(out_fn, "w") as OUT:
        OUT.write("SPECIES\t" + "\t".join(codons) + "\n")
        for species_id in species_ids:
            c_vals = ["0" for c in range(len(codons))]
            codon_tbl = codon_tbls[species_id]
            codon_total = sum([codon_tbl[c] for c in codon_tbl.keys()])
            for i, codon in enumerate(codons):
                c_vals[i] = str((codon_tbl[codon] / codon_total) * 100.0)
            OUT.write(species_id + "\t" + "\t".join(c_vals) + "\n")

