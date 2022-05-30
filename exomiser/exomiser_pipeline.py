import json
import hail as hl
import os
import re

import pandas as pd
from step_pipeline import pipeline, Backend, Localize, Delocalize

"""
root@02045285786b:/# java -jar exomiser-cli-13.0.1/exomiser-cli-13.0.1.jar
usage: java -jar exomiser-cli-{build.version}.jar [...]
    --analysis <file>         Path to analysis script file. This should be
                              in yaml format.
    --analysis-batch <file>   Path to analysis batch file. This should be
                              in plain text file with the path to a single
                              analysis script file in yaml format on each
                              line.
    --assembly <assembly>     Genome assembly of sample VCF file. Either
                              'GRCh37' or 'GRCh38'
    --batch <file>            Path to cli batch file. This should be in
                              plain text file with the cli input for an
                              analysis on each line.
 -h,--help                    Shows this help
    --job <file>              Path to job file. This should be in JSON or
                              YAML format.
    --output <string>         Path to outputOptions file. This should be
                              in JSON or YAML format.
    --ped <file>              Path to sample PED file.
    --preset <string>         The Exomiser analysis preset for the input
                              sample. One of 'exome', 'genome' or
                              'phenotype-only'
    --sample <file>           Path to sample or phenopacket file. This
                              should be in JSON or YAML format.
    --vcf <file>              Path to sample VCF file.
"""

DOCKER_IMAGE = "weisburd/exomiser@sha256:e617537f10cfd1a4ff85308cea880fd1a2f8d89672aafcf4e227ec4da4c48415"

def define_args(pipeline):
    """Define command-line args for the Exomiser pipeline.

    Args:
        pipeline (step_pipeline._Pipeline): The step_pipeline pipeline object.
    """
    parser = pipeline.get_config_arg_parser()
    grp = parser.add_argument_group("Exomiser")
    #grp.add_argument("--assembly", choices={"GRCh37", "GRCh38"}, default="GRCh38")  # pipeline currently supports only GRCh38
    grp.add_argument("-e", "--exomiser-data-dir",
                     help="Google Storage path of the Exomiser reference data directory. "
                          "For example 'gs://reference-data-bucket/exomiser-cli-13.0.0",
                     default="gs://lirical-reference-data/exomiser-cli-13.0.0")

    grp.add_argument("-o", "--output-dir",
                     help="Google Storage directory where to write Exomiser output",
                     required=True)

    grp.add_argument("--prioritiser",
                     choices={"phive", "phenix", "exomewalker", "hiphive"},
                     default="hiphive",
                     help="Phenotype similarity algorithm")

    grp.add_argument("--vcf",
                     help="Google Storage path that contains single-sample VCFs referenced by the phenopackets. More "
                          "than one path can be provided by specifying this argument more than once. Also each path "
                          "can optionally contain wildcards (*).",
                     required=True,
                     action="append")

    grp.add_argument("phenopacket_paths",
                     nargs="+",
                     help="Google Storage path of phenopacket JSON files to process. More than one path can be "
                          "specified. Also each path can optionally contain wildcards (*).")

    grp.add_argument("-s", "--sample-id", help="Optionally, process only this sample id. Useful for testing.")


def parse_args(pipeline):
    """Define and parse command-line args.

    Args:
        pipeline (step_pipeline._Pipeline): The step_pipeline pipeline object.

    Return:
         argparse.Namespace: parsed command-line args
         pandas.DataFrame: DataFrame with 1 row per phenopacket and columns: "sample_id", "phenopacket_path", "vcf_path"
    """

    define_args(pipeline)
    args = pipeline.parse_args()

    parser = pipeline.get_config_arg_parser()

    # initialize hail with workaround for Hadoop bug involving requester-pays buckets:
    # https://discuss.hail.is/t/im-encountering-bucket-is-a-requester-pays-bucket-but-no-user-project-provided/2536/2
    def get_bucket(path):
        if not path.startswith("gs://"):
            parser.error(f"{path} must start with gs://")
        return re.sub("^gs://", "", path).split("/")[0]

    all_buckets = {
        get_bucket(path) for path in [args.exomiser_data_dir] + args.phenopacket_paths + args.vcf
    }

    hl.init(log="/dev/null", quiet=True, idempotent=True, spark_conf={
        "spark.hadoop.fs.gs.requester.pays.mode": "CUSTOM",
        "spark.hadoop.fs.gs.requester.pays.buckets": ",".join(all_buckets),
        "spark.hadoop.fs.gs.requester.pays.project.id": args.gcloud_project,
    })

    # validate input paths
    def check_paths(paths):
        checked_paths = []
        for path in paths:
            if not path.startswith("gs://"):
                parser.error(f"Path must start with gs:// {path}")
            current_paths = [r["path"] for r in hl.hadoop_ls(path)]
            if not current_paths:
                parser.error(f"{path} not found")
            checked_paths += current_paths
        return checked_paths

    phenopacket_paths = check_paths(args.phenopacket_paths)
    vcf_paths = check_paths(args.vcf)

    # create DataFrame of phenopackets to process, with columns: "sample_id", "phenopacket_path", "vcf_path"
    rows = []
    requested_sample_id_found = False
    print(f"Processing {len(phenopacket_paths)} phenopacket(s)")
    for phenopacket_path in phenopacket_paths:
        print(f"Parsing {phenopacket_path}")
        with hl.hadoop_open(phenopacket_path, "r") as f:
            phenopacket_json = json.load(f)
            sample_id = phenopacket_json.get("subject", {}).get("id")
            if args.sample_id:
                if args.sample_id != sample_id:
                    continue
                else:
                    requested_sample_id_found = True

            if sample_id is None:
                parser.error(f"{phenopacket_path} is missing a 'subject' section")

            if ("htsFiles" not in phenopacket_json or not isinstance(phenopacket_json["htsFiles"], list) or
                    "uri" not in phenopacket_json["htsFiles"][0]):
                parser.error(f"{phenopacket_path} is missing an 'htsFiles' section with a VCF uri")
            vcf_filename = phenopacket_json["htsFiles"][0]["uri"].replace("file:///", "")

            matching_vcf_paths = [vcf_path for vcf_path in vcf_paths if vcf_filename in vcf_path]
            if not matching_vcf_paths:
                print(f"WARNING: Couldn't find {vcf_filename} referred to by {phenopacket_path}. Skipping...")
                continue
            vcf_path = matching_vcf_paths[0]

            hpo_ids = [d['type']['id'] for d in phenopacket_json["phenotypicFeatures"]]

        rows.append({
            "sample_id": sample_id,
            "phenopacket_path": phenopacket_path,
            "vcf_path": vcf_path,
            "hpo_ids": str(hpo_ids),
        })

        if requested_sample_id_found:
            break

    metadata_df = pd.DataFrame(rows)

    return args, metadata_df


def main():
    bp = pipeline("Exomiser", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    args, metadata_df = parse_args(bp)

    for _, row in metadata_df.iterrows():
        s1 = bp.new_step(f"Exomiser: {row.sample_id}", arg_suffix="exomiser",
                         image=DOCKER_IMAGE, cpu=1, storage="85Gi", memory="highmem",
                         localize_by=Localize.COPY, delocalize_by=Delocalize.COPY,
                         output_dir=args.output_dir,
        )

        s1.switch_gcloud_auth_to_user_account()
        #phenopacket_input = s1.input(row.phenopacket_path)

        vcf_input = s1.input(row.vcf_path)
        exomiser_hg38_dir_input = s1.input(os.path.join(args.exomiser_data_dir, "2109_hg38"))
        exomiser_phenotype_dir_input = s1.input(os.path.join(args.exomiser_data_dir, "2109_phenotype"))

        s1.command("cd /exomiser-cli-13.0.1/")
        s1.command("set -ex")

        s1.command(f"ln -s {exomiser_hg38_dir_input} {exomiser_hg38_dir_input.filename}")
        s1.command(f"ln -s {exomiser_phenotype_dir_input} {exomiser_phenotype_dir_input.filename}")

        hpo_ids = row.hpo_ids.replace("'", '\\"')
        s1.command(f"""python3 <<- EOF
with open("/exome_analysis.yml", "rt") as f: contents = f.read()
contents = contents.replace("VCF_PATH", "{vcf_input}")
contents = contents.replace("HPO_IDS", "{hpo_ids}")
with open("/exome_analysis.yml", "wt") as f: f.write(contents)             
EOF""")

        s1.command("cat /exome_analysis.yml")
        s1.command("mkdir -p results")
        s1.command("java -jar exomiser-cli-13.0.1.jar --analysis /exome_analysis.yml")
        s1.command("ls -ltr")
        s1.command("ls -ltr results")
        s1.output(f"results/{row.sample_id}_exomiser.json")
        s1.output(f"results/{row.sample_id}_exomiser.html")

    # run the pipeline
    bp.run()


if __name__ == "__main__":
    main()
