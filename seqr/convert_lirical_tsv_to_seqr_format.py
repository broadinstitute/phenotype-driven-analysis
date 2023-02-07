"""
This script takes the LIRICAL output TSV and converts it to a format that can be uploaded
to seqr (pending implementation of https://github.com/broadinstitute/seqr/issues/2742).
"""

import argparse
import pandas as pd
import re

from utils.gene_ids import get_entrez_to_ensembl_id_map


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--max-results", default=200, type=int, help="Output at most this many LIRICAL results")
    p.add_argument("--output-path", help="Output table path")
    p.add_argument("lirical_tsv", nargs="+", help="One or more LIRICAL output .tsv files")
    args = p.parse_args()

    output_rows = []
    gene_id_map = get_entrez_to_ensembl_id_map()
    for lirical_tsv in args.lirical_tsv:
        print(f"Parsing {lirical_tsv}")
        with open(lirical_tsv) as f:
            next(f)
            sample_id_line = next(f)
            sample_id = sample_id_line.rstrip().replace("! Sample: ", "")

        df = pd.read_table(lirical_tsv, comment="!", dtype=str)
        #print(f"Read {len(df)} rows from {lirical_tsv}")

        for _, row in df.iterrows():
            entrez_gene_id = str(int(re.sub("^NCBIGene:", "", row.entrezGeneId)))
            ensemble_gene_id = gene_id_map.get(entrez_gene_id)
            if not ensemble_gene_id:
                print(f"WARNING: Entrez gene Id {row.entrezGeneId} doesn't have a matching Ensembl gene Id")

            if int(row["rank"]) > args.max_results:
                #print(f"Skipping row(s) after row #{args.max_results}: {row.to_dict()}")
                continue

            output_rows.append({
                "tool": "lirical",
                "sampleId": sample_id,
                "rank": int(row["rank"]),
                "geneId": ensemble_gene_id,
                "diseaseId": row.diseaseCurie,   # "OMIM:130720"
                "diseaseName": row.diseaseName,  # "Lateral meningocele syndrome"
                "scoreName1": "post_test_probability",
                "score1": float(row.posttestprob.strip("%")),
                "scoreName2": "compositeLR",
                "score2": float(row.compositeLR.replace(",", "")),
                "scoreName3": None,
                "score3": None,
                "project": "Rare Genomes Project_Genomes_HMB" if sample_id not in "RGP_2035_3, RGP_2040_3, RGP_2058_1, RGP_2058_3, RGP_2061_3, RGP_2066_3, RGP_2069_3, RGP_2075_3, RGP_2095_3, RGP_2098_3, RGP_2107_3, RGP_2109_3, RGP_2113_3, RGP_2120_3, RGP_2122_3, RGP_2133_3".split(", ") else "Rare Genomes Project_Genomes_GRU", # TODO resolve this some other way
            })

    if not args.output_path:
        args.output_path = f"combined_lirical_results_for_{len(args.lirical_tsv)}_samples.for_seqr.tsv",

    output_df = pd.DataFrame(output_rows)
    output_df.to_csv(args.output_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(output_df)} rows to {args.output_path}")


if __name__ == "__main__":
    main()
