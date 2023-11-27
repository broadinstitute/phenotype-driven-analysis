"""
This script takes the LIRICAL output TSV and converts it to a format that can be uploaded
to seqr (pending implementation of https://github.com/broadinstitute/seqr/issues/2742).
"""

import argparse
import pandas as pd
import re
import os
from glob import glob

from utils.gene_ids import get_entrez_to_ensembl_id_map


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--max-results", default=200, type=int,
                   help="Output at most this many LIRICAL results")
    p.add_argument("--output-path", help="Output table path")
    p.add_argument("--sample-metadata",
                   help="Local path to a metadata text file that matches project IDs "
                        "to sample IDs for seqr upload.",
                   required=True)
    p.add_argument("lirical_tsv", nargs="+",
                   help="Local path(s) to one or more LIRICAL output .tsv files. Can "
                        "use regex wildcard expression \"*\" to denote multiple files "
                        "with similar paths.")
    args = p.parse_args()

    metadata = pd.read_table(args.sample_metadata)
    metadata = metadata.rename(columns={"individual_id": "sample_id"})

    output_rows = []
    gene_id_map = get_entrez_to_ensembl_id_map()

    # Check for wildcard insertion.
    if len(args.lirical_tsv) == 1 and "*" in args.lirical_tsv[0]:
        lirical_tsv_all = glob(args.lirical_tsv[0])
        print(lirical_tsv_all)
    else:
        lirical_tsv_all = args.lirical_tsv

    for lirical_tsv in lirical_tsv_all:
        print(f"Parsing {lirical_tsv}")
        with open(lirical_tsv) as f:
            next(f)
            sample_id_line = next(f)
            sample_id = sample_id_line.rstrip().replace("! Sample: ", "")

        df = pd.read_table(lirical_tsv, comment="!", dtype=str)
        # print(f"Read {len(df)} rows from {lirical_tsv}")

        for _, row in df.iterrows():
            entrez_gene_id = str(int(re.sub("^NCBIGene:", "", row.entrezGeneId)))
            ensemble_gene_id = gene_id_map.get(entrez_gene_id)
            if not ensemble_gene_id:
                print(
                    f"WARNING: Entrez gene Id {row.entrezGeneId} doesn't have a matching Ensembl gene Id")

            if int(row["rank"]) > args.max_results:
                # print(f"Skipping row(s) after row #{args.max_results}: {row.to_dict()}")
                continue

            output_rows.append({
                "tool": "lirical",
                "sampleId": sample_id,
                "rank": int(row["rank"]),
                "geneId": ensemble_gene_id,
                "diseaseId": row.diseaseCurie,  # "OMIM:130720"
                "diseaseName": row.diseaseName,  # "Lateral meningocele syndrome"
                "scoreName1": "post_test_probability",
                "score1": float(row.posttestprob.strip("%")),
                "scoreName2": "compositeLR",
                "score2": float(row.compositeLR.replace(",", "")),
                "scoreName3": None,
                "score3": None,
                "project": metadata[metadata["sample_id"] == sample_id][
                    "project_name"].values[0],
                # TODO resolve this some other way
            })

    if not args.output_path:
        args.output_path = f"combined_lirical_results_for_{len(args.lirical_tsv)}_samples.for_seqr.tsv",

    output_df = pd.DataFrame(output_rows)
    output_df.to_csv(args.output_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(output_df)} rows to {args.output_path}")


if __name__ == "__main__":
    main()
