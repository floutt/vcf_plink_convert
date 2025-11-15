import argparse
from EIGENTOOLS import PackedAncestryMap
from math import nan


def flatten(xss):
    """
    flattens a list of lists into a list of elements
    """
    return [x for xs in xss for x in xs]


def unique(sequence):
    """
    Gets all of the unique elements of a list while maintianing the original order
    """
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("--file_prefix", help="PACKEDANCESTRYMAP file prefix")
    parser.add_argument("--pops", help="File with newline-separated list of populations to include")

    args = parser.parse_args()

    print("##fileformat=VCFv4.2")
    print("##FILTER=<ID=PASS,Description=\"All filters passed\">")

    pam = PackedAncestryMap(file_prefix=args.file_prefix)
    # get pops to include
    pops = []
    with open(args.pops) as f:
        for line in f:
            pops.append(line.strip())

    inds = flatten([pam.ind_info.get_label_indices(pop) for pop in pops])

    for chrom in unique(pam.snp_info.chrom):
        print("##contig=<ID=%s>" % chrom)

    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    # print header
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", end="")
    print("\t".join([pam.ind_info.ind_name[i] for i in inds]))
    map_dosage = {2: "\t0/0", 1: "\t0/1", 0: "\t1/1", nan: "\t./."}

    for e in pam:
        sinfo = e.get_SNP_Info()
        print("\t".join([sinfo.chrom[0], "%i" % sinfo.pos[0], sinfo.var_name[0],
                         sinfo.ref[0], sinfo.alt[0], "999", "PASS", ".", "GT"]), end="")
        for c in [e.geno[i] for i in inds]:
            print(map_dosage[c], end="")
        print("")
