import argparse
import hail as hl
import json
import os
import pandas as pd
import re

from jycm.jycm import YouchamaJsonDiffer
from pprint import pformat

hl.init(log="/dev/null", quiet=True, idempotent=True)

PHENOPACKET_TEMPLATE = """
  "id": "%(individual_id)s",
  "subject": {
    "id": "%(individual_id)s",
    "sex": "%(sex)s"
  },
  "phenotypicFeatures": %(phenotypic_features_list)s,
  "metaData": {
    "createdBy": "export_phenopackets.py",
    "resources": [{
      "id": "hp",
      "name": "human phenotype ontology",
      "url": "http://purl.obolibrary.org/obo/hp.owl",
      "version": "unknown HPO version",
      "namespacePrefix": "HP",
      "iriPrefix": "http://purl.obolibrary.org/obo/HP_"
    }, {
      "id": "eco",
      "name": "Evidence and Conclusion Ontology",
      "url": "http://purl.obolibrary.org/obo/eco.owl",
      "version": "2019-10-16",
      "namespacePrefix": "ECO",
      "iriPrefix": "http://purl.obolibrary.org/obo/ECO_"
    }],
    "phenopacketSchemaVersion": "1.0.0"
  }
"""

PHENOPACKET_TEMPLATE_WITH_VCF = PHENOPACKET_TEMPLATE.strip() + """,
 "htsFiles": [{
    "uri": "%(vcf_path)s",
    "htsFormat": "VCF",
    "genomeAssembly": "%(genome_version)s"
  }]
"""

"""
    "ageAtCollection": {
      "age": "P17Y"
    },
"""


def parse_args():
    p = argparse.ArgumentParser(description="This script generates phenopackets")
    p.add_argument("-s", "--sample-id", nargs="+", help="Process specific sample id(s)")
    p.add_argument("--vcf-dir",
                   default="gs://bw-proj/seqr-bw/single_sample_vcfs/RDG_WGS_Broad_Internal/filtered_gnomad_popmax_0.01/single_sample_vcfs/",
                   help="Directory that contains VCFs for each sample. The VCF filename is expected to be [sample_id].vcf.*")

    #p.add_argument("--output-dir", default=".", help="Phenopacket output directory")
    p.add_argument("sample_table", help=".tsv sample table with columns ''")
    args = p.parse_args()

    #if not os.path.isdir(args.output_dir):
    #    p.error(f"{args.output_dir} is not a valid directory")

    #if args.vcf_dir and not os.path.isdir(args.vcf_dir):
    #    p.error(f"{args.vcf_dir} is not a valid directory")

    if not os.path.isfile(args.sample_table):
        p.error(f"{args.sample_table} doesn't exist")

    return args


class Phenotype:
    def __init__(self, hpo_id, category, label, present_or_absent):
        self.hpo_id = hpo_id
        self.category = category
        self.label = label
        self.present_or_absent = present_or_absent

    def __str__(self):
        return f"HP:{self.hpo_id}: {self.label}{' ABSENT' if self.present_or_absent == 'ABSENT' else ''}"

    def __repr__(self):
        return str(self)


def parse_phenotypes_from_string(phenotypes, is_absent):
    for i, phenotype in enumerate(phenotypes.split(", HP")):
        match = re.match("H?P?:([0-9]+):([^:]+):(.+)", phenotype)
        if not match:
            print(f"ERROR: couldn't parse '{phenotype}' from phenotype list: {phenotypes}. Skipping this term...")
            continue
        hpo_id, hpo_category, hpo_label = match.groups()
        yield {
            "type": {
                "id": f"HP:{hpo_id}",
                "label": hpo_label,
            },
            "negated": is_absent,
            "evidence": [{
                "evidenceCode": {
                    "id": "ECO:0000302",
                    "label": "author statement used in manual assertion"
                }
            }]
        }

    # double-check that phenotypes were parsed correctly
    hpo_id_count = phenotypes.count("HP:")
    if i+1 != hpo_id_count:
        raise ValueError(f"Parsed only {i+1} out of {hpo_id_count} phenotypes from phenotype list {phenotypes}")


def main():
    args = parse_args()
    df = pd.read_table(args.sample_table)
    df = df.rename(columns={"individual_id": "sample_id"})
    if args.sample_id:
        df = df.loc[df["sample_id"].isin(args.sample_id)]
        # df = df[df.sample_id.isin(args.sample_id)]

    vcf_paths = hl.hadoop_ls(os.path.join(args.vcf_dir, "*.vcf.*gz"))

    sample_ids_with_vcf = {
        os.path.basename(p["path"]).replace(".vcf.gz", "").replace(".vcf.bgz", ""): p["path"] for p in vcf_paths
    }
    print(f"Found {len(sample_ids_with_vcf)} VCF paths")
    print(f"Processing {len(df)} rows from {args.sample_table}")

    existing_phenopacket_paths_to_check = {
        os.path.join(os.path.dirname(os.path.dirname(vcf_path)), "phenopackets/*.phenopacket.json")
        for vcf_path in sample_ids_with_vcf.values()
    }
    existing_phenopacket_paths = set()
    for p in existing_phenopacket_paths_to_check:
        existing_phenopacket_paths |= {p["path"] for p in hl.hadoop_ls(p)}

    for _, row in df.iterrows():
        if row.sample_type != "WGS":
            continue

        vcf_path = None
        if row.sample_id in sample_ids_with_vcf:
            vcf_path = sample_ids_with_vcf[row.sample_id]

        if not vcf_path:
            continue

        phenopacket_output_path = os.path.join(
            os.path.dirname(os.path.dirname(vcf_path)),
            "phenopackets", f"{row.sample_id}.phenopacket.json")

        #if row.affected != "Affected":
        #    print("   Skipping unaffected")
        #    continue

        print(f"Processing: {row.sample_id}")

        #vcf_path = os.path.join(args.vcf_dir, f"{row.sample_id}.vcf")
        #if not os.path.isfile(vcf_path):
        #    print(f"WARNING: {vcf_path} not found. Skipping {row.sample_id}...")
        #    #continue

        phenotypic_features_list = []
        if not pd.isna(row.phenotypes):
            present_phenotypic_features = list(parse_phenotypes_from_string(row.phenotypes, is_absent=False))
            phenotypic_features_list += present_phenotypic_features
            print(f"Parsed {len(present_phenotypic_features)} phenotypes")
        if not pd.isna(row.absent_phenotypes):
            absent_phenotypic_features = list(parse_phenotypes_from_string(row.absent_phenotypes, is_absent=True))
            phenotypic_features_list += absent_phenotypic_features
            print(f"Parsed {len(absent_phenotypic_features)} absent phenotypes")

        if len(phenotypic_features_list) == 0:
            print(f"No phenotypes found for {row.sample_id}. Skipping...")
            continue

        if row.sex.upper().startswith("F"):
            sex = "FEMALE"
        elif row.sex.upper().startswith("M"):
            sex = "MALE"
        else:
            sex = "UNKNOWN_SEX"

        parameters = {
            "individual_id": row.sample_id,
            "sex": sex,
            "genome_version": "hg38",
            "phenotypic_features_list": json.dumps(phenotypic_features_list, indent=2),
            "vcf_path": f"file:///{os.path.basename(vcf_path).replace('.gz', '').replace('.bgz', '')}",
        }

        phenopacket_string = PHENOPACKET_TEMPLATE_WITH_VCF % parameters
        phenopacket_json = json.loads("{" + phenopacket_string + "}")

        if phenopacket_output_path in existing_phenopacket_paths:
            # check whether phenotypes changed:
            with hl.hadoop_open(phenopacket_output_path, "r") as f:
                previous_phenopacket = json.load(f)

            ycm = YouchamaJsonDiffer(previous_phenopacket, phenopacket_json)
            json_diff = ycm.get_diff(no_pairs=True)

            if len(json_diff) == 0:
                print(f"Skipping existing phenopacket path: {phenopacket_output_path}")
                continue

            print(f"{row.sample_id} phenopacket ({phenopacket_output_path}) has changed:")
            for value_change in json_diff.get("value_changes", []):
                left_path = value_change.get('left_path', '.')
                right_path = value_change.get('right_path', '.')
                print(f"Previous [{left_path}]: ", pformat(value_change.get("old", '.')))
                print(f"Current [{right_path}]: ", pformat(value_change.get("new", '.')))
                print("------")

        with hl.hadoop_open(phenopacket_output_path, "w") as f:
            json.dump(phenopacket_json, f, indent=3)
        print(f"Wrote out {phenopacket_output_path}")

        # PHENOPACKET_TEMPLATE_WITH_VCF


if __name__ == "__main__":
    main()


"""
['sample_id',
 'genome_version',
 'project_guid',
 'project_name',
 'family_guid',
 'individual_guid',
 'father_id',
 'father_guid',
 'mother_id',
 'mother_guid',
 'sex',
 'affected',
 'analysis_status',
 'phenotypes',
 'absent_phenotypes',
 'disorders',
 'notes',
 'filter_flags',
 'assigned_analyst',
 'population',
 'coded_phenotype',
 'bucket',
 'cram_path',
 'crai_path',
 'cram_and_crai_exist',
 'sample_type',
 'sample_id_from_filename']
 """
#%%

{
    "individual_id": "",
    "sex": "FEMALE",
    "phenotypic_features_list": [],
    "vcf_uri": "file://Users/weisburd/project_phenotypes/2022_02_24__test_lirical/vcfs/RGP_1186_3.vcf",
    "genome_version": "hg38",
}

# get RGP stats by family
#df_RGP = df[df.sample_id.str.startswith("RGP_")]
#df.analysis_status.isin({'Analysis in Progress', 'Reviewed, currently pursuing candidates', 'Reviewed, no clear candidate',})
#df.analysis_status.str.contains("Solved")
#df_RGP_solved = df_RGP[df.analysis_status.str.contains("Solved")]
