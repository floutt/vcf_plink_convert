from EIGENTOOLS import SNP_Info, Ind_Info, PackedAncestryMapWriter
from math import nan
import gzip
import argparse


# PACKEDANCESTRYMAP files are encode with the REFERENCE dosage not the
# alternative
GT_MAP = {"0/0": 2, "1/0": 1, "0/1": 1, "1/1": 0, "./.": nan, ".": nan}


def write_snp_ind(vcf_file, snp_out, ind_out, ind_label, ind_sex=None):
    f = None
    if vcf_file.endswith(".gz"):
        f = gzip.open(vcf_file, "rt")
    else:
        f = open(vcf_file)

    f_out_snp = open(snp_out, "w+")
    f_out_ind = open(ind_out, "w+")

    ind_name = None
    ind_sex = ["U"] * len(ind_label) if ind_sex is None else ind_sex
    for line in f:
        if line.startswith("#CHROM"):
            # write ind file
            ind_name = line.split()[9:]
            if not (len(ind_sex) == len(ind_label) == len(ind_name)):
                raise Exception("Length mismatch between individual file information")
            for name, sex, label in zip(ind_name, ind_sex, ind_label):
                f_out_ind.write("%s\t%s\t%s\n" % (name, sex, label))
            f_out_ind.close()
            continue
        elif line.startswith("#"):
            continue

        elems = line.split()
        chrom = elems[0].replace("chr", "")  # remove chr symbol
        pos = elems[1]
        ref = elems[3]
        alt = elems[4]
        # skip if multillelic
        if "," in alt:
            continue
        snp_id = elems[2] if elems[2].startswith("rs") else "_".join(["snp", chrom, pos])
        f_out_snp.write("\t".join([snp_id, chrom, "0.0", pos, ref, alt]) + "\n")
    f_out_snp.close()
    f.close()


def vcf_to_pam(vcf_file, pam_prefix, ind_label, ind_sex=None):
    # write snp and ind files from the VCF
    write_snp_ind(vcf_file, pam_prefix + ".snp", pam_prefix + ".ind",
                  ind_label, ind_sex)

    snp_info = SNP_Info(pam_prefix + ".snp")
    ind_info = Ind_Info(pam_prefix + ".ind")

    paw = PackedAncestryMapWriter(snp_info, ind_info, file_prefix=pam_prefix,
                                  write_snp=False, write_ind=False)

    f = None
    if vcf_file.endswith(".gz"):
        f = gzip.open(vcf_file, "rt")
    else:
        f = open(vcf_file)

    for line in f:
        if line.startswith("#"):
            continue
        elems = line.split()
        alt = elems[4]
        # skip if multiallelic
        if "," in alt:
            continue
        fmt_info = elems[8].split(":")
        gt_idx = fmt_info.index("GT")
        genotypes = [x.split(":")[gt_idx].replace("|", "/") for x in elems[9:]]
        dosages = [GT_MAP[x] for x in genotypes]
        paw.write_record(dosages)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("--vcf", help="VCF file input")
    parser.add_argument("--out_prefix", help="output PACKEDANCESTRYMAP file prefix")
    parser.add_argument("--ind_pop", nargs="*", default=None, help="Population labels for each of the individuals in the VCF file")
    parser.add_argument("--ind_sex", nargs="*", default=None, help="Sex for each of the individuals in the VCF file")

    args = parser.parse_args()

    vcf_to_pam(args.vcf, args.out_prefix, args.ind_pop, args.ind_sex, args.flip)
