"""After single_sample_vcf_pipeline_for_large_mt.py is done, this script can be used to convert the TSVs to VCFs"""

import hail as hl
import os
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/large_mt_to_single_sample_vcfs@sha256:392c600d7346b8ab65819239da7ab7657b9b101c00e283e9053397693c198aa0"


def parse_args(pipeline):
    """Define command-line args for the pipeline.

    Args:
        pipeline (step_pipeline._Pipeline): The step_pipeline pipeline object.

    Return:
        3-tuple: list of tsv paths to process, output directory, add_info_field arg value
    """

    parser = pipeline.get_config_arg_parser()
    parser.add_argument("--add-info-field", help="Add information to the info field", action="store_true")
    parser.add_argument("-o", "--output-dir", help="Where to write output vcfs. If not specified, they will be written "
                                                   "to the directory that contains the input tsvs.")
    parser.add_argument("-s", "--sample-id", help="Optionally, process only this sample id. Useful for testing.")
    parser.add_argument("tsv_path", help="Path of one or more tsv files generated by hl.experimental.export_entries_by_col")
    args = pipeline.parse_args()

    # initialize hail with workaround for Hadoop bug involving requester-pays buckets:
    # https://discuss.hail.is/t/im-encountering-bucket-is-a-requester-pays-bucket-but-no-user-project-provided/2536/2
    def get_bucket(path):
        if not path.startswith("gs://"):
            parser.error(f"{path} must start with gs://")
        return re.sub("^gs://", "", path).split("/")[0]

    all_buckets = {
        get_bucket(path) for path in [args.tsv_path, args.output_dir] if path is not None
    }

    hl.init(log="/dev/null", quiet=True, idempotent=True, spark_conf={
        "spark.hadoop.fs.gs.requester.pays.mode": "CUSTOM",
        "spark.hadoop.fs.gs.requester.pays.buckets": ",".join(all_buckets),
        "spark.hadoop.fs.gs.requester.pays.project.id": args.gcloud_project,
    })

    tsv_paths = [p["path"] for p in hl.hadoop_ls(args.tsv_path) if p["path"].endswith(".tsv.bgz")]
    if not tsv_paths:
        parser.error(f" No .tsv.bgz files found at {args.tsv_path}")

    if args.sample_id:
        tsv_paths = [p for p in tsv_paths if get_sample_id(p) == args.sample_id]
        if not tsv_paths:
            parser.error(f"No .tsv.bgz file found for sample id {args.sample_id}")

    if not args.output_dir:
        args.output_dir = os.path.dirname(tsv_paths[0])

    return tsv_paths, args.output_dir, args.add_info_field


def get_sample_id(tsv_path):
    return re.sub(".tsv.bgz$", "", os.path.basename(tsv_path))


def main():
    bp = pipeline("convert_tsvs_to_vcfs", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    tsv_paths, output_dir, add_info_field = parse_args(bp)

    for tsv_path in tsv_paths:
        sample_id = get_sample_id(tsv_path)
        s1 = bp.new_step(
            f"convert_tsv_to_vcf: {sample_id}",
            image=DOCKER_IMAGE,
            cpu=0.25,
            storage=5,
            localize_by=Localize.HAIL_BATCH_GCSFUSE,
            delocalize_by=Delocalize.COPY,
            output_dir=output_dir)

        input_tsv = s1.input(tsv_path)

        s1.command("set -x")
        cmd = f"time python3 convert_tsv_to_vcf.py --vcf-header /vcf_header.txt {input_tsv}"
        if add_info_field:
            cmd += "--add-info-field"
        s1.command(cmd)

        s1.command(f"bgzip {sample_id}.vcf")
        s1.command(f"tabix {sample_id}.vcf.gz")

        s1.output(f"{sample_id}.vcf.gz")
        s1.output(f"{sample_id}.vcf.gz.tbi")

    # run the pipeline
    bp.run()


if __name__ == "__main__":
    main()