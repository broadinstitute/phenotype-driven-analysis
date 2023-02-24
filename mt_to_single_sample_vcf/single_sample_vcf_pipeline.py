"""This script exports single-sample VCF(s) from a matrix table (mt) by creating Hail Batch jobs each of which filters
the matrix table (mt) to a single sample and then exports a VCF.
"""

import configargparse
import hail as hl
import hailtop.batch as hb
import os
import re


def parse_args():
    p = configargparse.ArgumentParser(default_config_files=["~/.step_pipeline"], ignore_unknown_config_file_keys=True)
    p.add_argument("-o", "--output-dir", help="Google Storage path where to write the single-sample VCF(s)",
                   required=True)
    p.add_argument("-s", "--sample-id", help="Sample id(s) to export from the Matrix Table. If not specified, all "
                   "samples will be exported.", action="append")
    p.add_argument("--add-info-field", help="Add INFO field to VCF based on Matrix Table fields generated by the "
                   "seqr loading pipeline", action="store_true")
    p.add_argument("--cpu", help="Number of CPUs to request for each job", default=8)
    p.add_argument("--batch-billing-project", help="Hail Batch billing project", required=True)
    p.add_argument("--batch-remote-tmpdir", help="Google Storage path for Hail Batch to use as a temp directory",
                   required=True)
    p.add_argument("matrix_table_paths", nargs="+", help="Google Storage path of one or more Hail Matrix Tables from "
                   "which to extract the single sample VCF(s)")
    args = p.parse_args()
    return p, args


def main():
    parser, args = parse_args()

    hl.init(log="/dev/null", idempotent=True, quiet=True)

    path_to_mt = {}
    for mt_path in args.matrix_table_paths:
        mt = hl.read_matrix_table(mt_path)
        path_to_mt[mt_path] = mt

    b = hb.Batch(
        default_python_image="hailgenetics/hail:0.2.109",   # from https://hub.docker.com/r/hailgenetics/hail/tags
        backend=hb.ServiceBackend(
            billing_project=args.batch_billing_project,
            remote_tmpdir=args.batch_remote_tmpdir,
        )
    )

    for mt_path, mt in path_to_mt.items():
        sample_ids = mt.s.collect()
        print(f"Found {len(sample_ids)} total samples in {mt_path}")
        if args.sample_id:
            sample_ids = set(sample_ids) & set(args.sample_id)
            print(f"Found {len(sample_ids)} out of {len(args.sample_id)} requested sample ids")
            if len(sample_ids) < len(set(args.sample_id)) and len(args.matrix_table_paths) == 1:
                parser.error(", ".join(set(args.sample_id) - set(sample_ids)) + f" sample id(s) not found in {mt_path}")

        existing_files = {p["path"] for p in hl.hadoop_ls(args.output_dir)}
        for i, sample_id in enumerate(sample_ids):
            sample_id = sample_id.replace(" ", "_")
            output_vcf_path = os.path.join(args.output_dir, f"{sample_id}.vcf")

            if f"{output_vcf_path}.gz" in existing_files or f"{output_vcf_path}.bgz" in existing_files:
                print(f"   {output_vcf_path}*gz already exists. Skipping {sample_id}...")
                continue

            print(f"   {i}: Processing {sample_id}: {output_vcf_path}.gz")
            j = b.new_python_job(name=f"{sample_id} vcf")
            j.cpu(args.cpu)
            j.call(export_vcf, mt_path, sample_id, f"{output_vcf_path}.bgz", args.cpu, args.add_info_field)

        b.run()


def export_vcf(mt_path, sample_id, output_vcf_path, num_cpu=8, add_info_field=False):
    """Runs hail to read in the matrix table and export a single sample to a vcf. Optionally, it computes INFO field
    annotations from Matrix Table fields generated by the seqr loading pipeline.

    Args:
        mt_path (str): Google Storage path of Matrix Table
        sample_id (str): The sample id to export to VCF
        output_vcf_path (str): Google Storage path where to write the VCF
        num_cpu (int): Number of CPUs available to the execution container.
        add_info_field (bool): Whether to compute an INFO field based on annotations in the Matrix Table.
    """
    hl.init(master=f"local[{num_cpu}]")
    mt = hl.read_matrix_table(mt_path)
    mt.describe()
    mt = mt.filter_cols(mt.s == sample_id)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref() & (mt.DP > 0)))

    if add_info_field:
        mt = mt.annotate_rows(info=hl.struct(
            cohort_AC=mt.AC,
            cohort_AF=hl.format("%.3f", mt.AF),
            cohort_AN=mt.AN,
            hgmd_class=mt.hgmd['class'],
            clinvar_allele_id=mt.clinvar.allele_id,
            clinvar_clinsig=mt.clinvar.clinical_significance,
            clinvar_gold_stars=mt.clinvar.gold_stars,
            consequence=mt.mainTranscript.major_consequence,
            gene_id=mt.mainTranscript.gene_id,
            gene=mt.mainTranscript.gene_symbol,
            CADD=hl.format("%.3f", mt.cadd.PHRED),
            eigen=hl.format("%.3f", mt.eigen.Eigen_phred),
            revel=hl.format("%.3f", hl.float(mt.dbnsfp.REVEL_score)),
            splice_ai=mt.splice_ai.delta_score,
            primate_ai=mt.primate_ai.score,
            exac_AF=hl.format("%.3f", mt.exac.AF_POPMAX),
            gnomad_exomes_AF=hl.format("%.3f", mt.gnomad_exomes.AF_POPMAX_OR_GLOBAL),
            gnomad_genomes_AF=hl.format("%.3f", mt.gnomad_genomes.AF_POPMAX_OR_GLOBAL),
            topmed_AF=hl.format("%.3f", mt.topmed.AF),
        ))

    hl.export_vcf(mt, output_vcf_path, tabix=True)


if __name__ == "__main__":
    main()