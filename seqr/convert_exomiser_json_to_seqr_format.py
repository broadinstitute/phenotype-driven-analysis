"""
This script takes the Exomiser output json file and converts it to a format that can be uploaded to seqr
(pending implementation of https://github.com/broadinstitute/seqr/issues/2742)
"""

import argparse
import json
import os
import pandas as pd
import re


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--max-results", default=200, type=int, help="Output at most this many Exomiser results")
    p.add_argument("--output-suffix", default="for_seqr", help="This suffix will be added to the input filename to make the output filename")
    p.add_argument("exomiser_json", nargs="+", help="One or more Exomiser output .tsv files")
    args = p.parse_args()

    for exomiser_json_path in args.exomiser_json:
        with open(exomiser_json_path) as f:
            exomiser_results = json.load(f)

        sample_id = re.sub("_exomiser.json$", "", os.path.basename(exomiser_json_path))

        output_rows = []
        for i, result in enumerate(exomiser_results):

            # result keys: [
            #   'geneSymbol', 'geneIdentifier', 'combinedScore', 'priorityScore', 'variantScore',
            #   'priorityResults', 'compatibleInheritanceModes', 'geneScores', 'variantEvaluations',
            # ]

            # result["geneIdentifier"] == { 'geneId': 'ENSG00000133612', 'geneSymbol': 'AGAP3', 'entrezId': '116988'}
            exomiser_rank = i + 1
            if exomiser_rank > args.max_results:
                #print(f"Skipping row(s) after row #{args.max_results}: {row.to_dict()}")
                continue

            ensembl_gene_id = result["geneIdentifier"]["geneId"]
            exomiser_score = result.get("combinedScore")
            if exomiser_score is None:
                print(f"WARNING: skipping result #{exomiser_rank} because it's missing an Exomiser score.")
                continue

            variant_score = result.get("variantScore")
            #priority_score = result["priorityScore"]

            hiphive_priority = result["priorityResults"]["HIPHIVE_PRIORITY"]
            phenotype_evidence_list = []
            phenotype_evidence_list += hiphive_priority.get("phenotypeEvidence", [])
            phenotype_evidence_list += hiphive_priority.get("ppiEvidence", [])
            phenotype_evidence_list += hiphive_priority.get("diseaseMatches", [])

            #phenotype_evidence_list = [r for r in phenotype_evidence_list if r["model"]["organism"] == "HUMAN"]

            phenotype_evidence_results = []
            for phenotype_evidence in phenotype_evidence_list:

                try:
                    phenotype_score = phenotype_evidence["score"]
                    disease_id = phenotype_evidence["model"]["disease"]["diseaseId"]
                    disease_name = phenotype_evidence["model"]["disease"]["diseaseName"]
                    phenotype_evidence_results.append((phenotype_score, disease_id, disease_name))
                except KeyError:
                    continue

            if not phenotype_evidence_results:
                phenotype_evidence_results.append((None, None, None))
                #pprint(result["priorityResults"]["HIPHIVE_PRIORITY"])

            for phenotype_score, disease_id, disease_name in sorted(list(set(phenotype_evidence_results)), reverse=True):
                output_rows.append({
                    "tool": "exomiser",
                    "sampleId": sample_id,
                    "rank": exomiser_rank,
                    "geneId": ensembl_gene_id,
                    "diseaseId": disease_id,   # "OMIM:130720"
                    "diseaseName": disease_name,  # "Lateral meningocele syndrome"
                    "scoreName1": "exomiser_score",
                    "score1": exomiser_score,
                    "scoreName2": "phenotype_score",
                    "score2": phenotype_score,
                    "scoreName3": "variant_score",
                    "score3": variant_score,
                })

        output_path = re.sub(".json$", f".{args.output_suffix}.tsv", exomiser_json_path)
        output_df = pd.DataFrame(output_rows)
        output_df.to_csv(output_path, sep="\t", index=False, header=True)
        print(f"Wrote {len(output_df)} rows to {output_path}")


if __name__ == "__main__":
    main()