import hail as hl
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/slivar@sha256:b97f4b57dd58d0f3ef9824fc209a17665073d97a5af37129ba09b6e3605edd51"


def parse_args(pipeline):
    """Define and parse command-line args.

    Args:
        pipeline (step_pipeline._Pipeline): The step_pipeline pipeline object.

    Return:
         argparse.Namespace: parsed command-line args
         pandas.DataFrame: DataFrame with 1 row per phenopacket and columns: "sample_id", "phenopacket_path", "vcf_path"
    """

    parser = pipeline.get_config_arg_parser()
    grp = parser.add_group("Prefilter VCF")
    grp.add_argument("-g", "--gnomad-af", type=float, default=0.01, help="Filter VCF to variants with gnomAD POPMAX "
                        "AF below this threshold.")
    grp.add_argument("-o", "--output-dir", help="Google Storage directory where to write the filtered VCF. "
                        "If not specified, the result vcf will be copied to the same directory as the input vcf")
    grp.add_argument("vcf_path", nargs="+", help="Google Storage path of VCF file(s) to filter")

    args = parser.parse_args()

    # initialize hail with workaround for Hadoop bug involving requester-pays buckets:
    # https://discuss.hail.is/t/im-encountering-bucket-is-a-requester-pays-bucket-but-no-user-project-provided/2536/2
    def get_bucket(path):
        if not path.startswith("gs://"):
            parser.error(f"{path} must start with gs://")
        return re.sub("^gs://", "", path).split("/")[0]

    all_buckets = {get_bucket(path) for path in args.vcf_path}

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

    check_paths(args.vcf_path)

    return args


def main():
    sp = pipeline("prefilter_vcf", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    args = parse_args(sp)

    # create DataFrame of phenopackets to process, with columns: "sample_id", "phenopacket_path", "vcf_path"
    if len(args.vcf_path) == 1:
        print(f"Processing {len(args.vcf_path[0])}")
    for vcf_path in args.vcf_path:
        print(f"Parsing {vcf_path}")

    for vcf_path in args.vcf_path:
        s1 = sp.new_step(f"prefilter vcf: {os.path.basename(vcf_path)}",
                         cpu=1, memory="standard", image=DOCKER_IMAGE,
                         localize_by=Localize.COPY, delocalize_by=Delocalize.COPY,
                         output_dir=args.output_dir or os.path.dirname(vcf_path))

        #s1.switch_gcloud_auth_to_user_account()
        input_vcf = s1.input(vcf_path)

        s1.command("cd /io/")
        s1.command("set -ex")
        s1.command(f"wget --quiet https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.fix.zip")

        cat_command = "zcat" if input_vcf.filename.endswith("gz") else "cat"
        s1.command(f"{cat_command} {input_vcf} | grep -v :NA: > without_NA.vcf")

        s1.command(f"java -jar /gatk.jar FixVcfHeader -I without_NA.vcf -O fixed_header.vcf")

        output_vcf_filename = input_vcf.filename.replace(".vcf", f".filtered_gnomad_popmax_{args.gnomad_af}.vcf")

        s1.command(f"/slivar expr "
            f"--js /slivar-functions.js "
            "-g gnomad.hg38.genomes.v3.fix.zip "
            f"--info 'INFO.gnomad_popmax_af < {args.gnomad_af}' "
            f"--vcf fixed_header.vcf "
            f"-o {output_vcf_filename} "
        )

        s1.command(f"tabix {output_vcf_filename}")

        s1.output(output_vcf_filename)
        s1.output(f"{output_vcf_filename}.tbi")

    # run the pipeline
    sp.run()


if __name__ == "__main__":
    main()