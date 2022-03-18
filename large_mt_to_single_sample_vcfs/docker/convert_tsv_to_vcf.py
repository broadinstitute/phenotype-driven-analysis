import argparse
import gzip
import json

#%%
p = argparse.ArgumentParser()
p.add_argument("--add-info-field", help="Add information to the info field", action="store_true")
p.add_argument("--vcf-header", help="VCF header file path", default="/vcf_header.txt")
p.add_argument("tsv_path", help="One of the .tsv.bgz files generated by hl.experimental.export_entries_by_col")
args = p.parse_args()

#%%
"""
args = argparse.Namespace()
args.add_info_field = True
root = "/Users/weisburd/code/phenotype-driven-analysis/large_mt_to_single_sample_vcfs"
args.tsv_path = f"{root}/part-0058.tsv.bgz"
args.vcf_header = f"{root}/docker/vcf_header.txt"
"""

#%%

# Read TSV header lines
print(f"Parsing {args.tsv_path}")

# parse sample id
f = gzip.open(args.tsv_path, "rt")
sample_id_json = json.loads(next(f).lstrip("#"))
sample_id = sample_id_json["s"]

# parse header line
# locus   alleles filters rsid    info    AD      DP      GQ      GT      MIN_DP  PGT     PID     PL      RGQ     SB
header = next(f).strip().split("\t")

# read VCF header
with open(args.vcf_header, "rt") as f_header:
    vcf_header_contents = f_header.read()

vcf_header_contents = vcf_header_contents.replace("[sample_id]", sample_id)

# convert TSV to VCF
output_path = f"{sample_id}.vcf"
fo = open(output_path, "wt")
fo.write(vcf_header_contents)

format = ["GT", "AD", "DP", "GQ", "PL"]

output_line_counter = 0
for i, line in enumerate(f, start=2):
    fields = line.strip().split("\t")
    row = dict(zip(header, fields))
    row["alleles"] = json.loads(row["alleles"])
    if len(row["alleles"]) != 2:
        raise ValueError(f"tsv row #{i+1} contains {len(row['alleles'])} alleles: {line}")

    if row["GT"] in {"0/0", "./.", "0|0", "0\\0"}:
        continue

    row["chrom"], row["pos"] = row["locus"].split(":")
    row["filters"] = ",".join(json.loads(row["filters"])) or "PASS"

    if args.add_info_field:
        row["info"] = {k: v for k, v in json.loads(row["info"]).items() if v is not None and v != "nul"}
        info_field = ";".join([f"{key}={value}" for key, value in row["info"].items()])
    else:
        info_field = "."

    ad = row["AD"] = ",".join(map(str, json.loads(row["AD"])))
    pl = row["PL"] = ",".join(map(str, json.loads(row["PL"])))

    genotype = [row[key] for key in format]

    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	[sample_id]
    vcf_line = "\t".join(map(str, [
        row["chrom"],
        row["pos"],
        ".",
        row["alleles"][0],
        row["alleles"][1],
        ".",
        row["filters"],
        info_field,
        ":".join(format),
        ":".join(genotype),
    ])) + "\n"

    output_line_counter += 1
    fo.write(vcf_line)

    # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12877
print(f"Wrote {output_line_counter} out of {i-1} ({100*output_line_counter/(i-1):0.1f}%) of records to {output_path} (after filtering out GT=0/0)")
f.close()
fo.close()
