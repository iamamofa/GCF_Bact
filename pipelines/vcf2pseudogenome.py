#!/usr/bin/env python

import sys
import pysam

def create_consensus_fasta(vcf_file, reference_file):
    # Open the reference genome
    reference = pysam.FastaFile(reference_file)

    # Open the VCF file
    vcf = pysam.VariantFile(vcf_file)

    # Prepare to write the consensus FASTA
    output = {}
    for rec in vcf.fetch():
        chrom = rec.chrom
        pos = rec.pos
        ref = rec.ref
        alt = rec.alts[0] if rec.alts else ref
        if chrom not in output:
            output[chrom] = list(reference.fetch(chrom))

        if ref != alt:
            output[chrom][pos - 1] = alt

    # Write to FASTA format
    for chrom, seq in output.items():
        print(f">{chrom}")
        print("".join(seq))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: vcf2pseudogenome.py <input_vcf> <reference>", file=sys.stderr)
        sys.exit(1)

    input_vcf = sys.argv[1]
    reference = sys.argv[2]

    create_consensus_fasta(input_vcf, reference)
