"""
This script gets all readviz file paths from the seqr database.
It checks that the files exist, and that confirms the genome version
based on the bam/cram header.

Then, it writes out a table with one row per seqr (individual, sample type) with the following columns:

        'genome_version',
        'project_guid',
        'project_name',
        'family_guid',
        'individual_id',
        'sex',
        'affected',
        'analysis_status',
        'population',
        'coded_phenotype',
        'bucket',
        'cram_path',
        'crai_path',
        "cram_and_crai_exist",
        'sample_type',
"""

#%%

import os
import pandas as pd
import re
import sys
from tqdm import tqdm
import hail as hl

hl.init(log="/dev/null", idempotent=True)

sys.path.append(os.path.expanduser("~/code/seqr"))

import django
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
django.setup()

from seqr.models import Project, Family, Individual, Sample, IgvSample, SavedVariant, VariantTag, VariantNote

from reference_data.models import Omim, HumanPhenotypeOntology
from metadata.gcloud_api_utils import get_genome_version_from_bam_or_cram_header, hadoop_exists_using_cache
from metadata.seqr_utils import impute_sample_type

AFFECTED_STATUS_LOOKUP = {
  'A': 'Affected',
  'N': 'Not Affected',
  'U': 'Unknown',
}

ANALYSIS_STATUS_LOOKUP = dict(Family.ANALYSIS_STATUS_CHOICES)

hpo_name_map = {hpo.hpo_id: hpo.name for hpo in HumanPhenotypeOntology.objects.all()}
hpo_name_map = {hpo.hpo_id: (hpo_name_map[hpo.category_id] + ": " + hpo.name) if hpo.category_id else hpo.name
                for hpo in HumanPhenotypeOntology.objects.all()}


print("Starting step1__export_cram_paths_from_seqr.py")

import datetime
date_string = datetime.datetime.now().strftime("%Y_%m_%d")
output_path = os.path.expanduser(f"~/code/sample_metadata/output_tables/seqr_cram_paths__{date_string}.txt_unfinished")
print("Writing to", output_path)
with open(output_path, "wt") as f:
    f.write("\t".join([
        'genome_version',
        'project_guid',
        'project_name',
        'family_id',
        'family_guid',
        'individual_id',
        'individual_guid',
        'father_id',
        'father_guid',
        'mother_id',
        'mother_guid',
        'sex',
        'affected',
        'analysis_status',

        'phenotypes',  # aka. features
        'absent_phenotypes',  # aka. absent_features
        'disorders',
        'notes',
        'filter_flags',
        'assigned_analyst',

        'population',
        'coded_phenotype',

        'post_discovery_omim_number',
        'discovery_gene_ids',
        'discovery_gene_symbols',

        'bucket',
        'cram_path',
        'crai_path',
        "cram_and_crai_exist",
        'sample_type',
    ]) + "\n")

    # go through all Individual records in seqr (there will be multiple records for each Individual if WES and WGS are in different projects)
    total_individuals = Individual.objects.count()
    missing_crams_counter = error_crams_counter = 0
    for individual in tqdm(Individual.objects.all(), total=total_individuals, unit=" individuals"):
        family = individual.family
        project = family.project

        project_guid_lowercase = project.guid.lower()
        project_name_lowercase = project.name.lower().replace(' ', '_')

        # skip some projects
        if "test_project" in project.name or "demo" in project_name_lowercase or "test-run" in project.name \
                or "_rna" in project_guid_lowercase or project.guid in (
                "R0326_engle_moebius", "R0270_marnero_acc_16_samples_r", "R0272_sample_1413_1_muscled",
                "R0486_cmg_gcnv", "R0391_wintermoran_june_2017", "R0224_ibd_273_samples_march_20",
                "R0283_sweetser_trio_1558955_li", "R0305_sigma_extremes", "R0417_wintermoran_november_201",
                "R0287_meei_pierce_samples", "R0548_bwh_pina_aguilar_dgap", "R0218_cdh_89s", "R0209_cdh2_292s",
                "R0551_seqr_loading_anvil_test", "R0585_210406_210506_210601_joi", "R0555_seqr_demo",
                "R0577_misc", "R0595_ibs_sequencing_project", "R0596_ibs_seq_merged",
        ):
            continue

        #if sample_type != "WGS":
        #    continue
        sample_type, file_path = impute_sample_type(individual)
        if sample_type is None:
            print("Guessing that sample type is WES")
            sample_type = "WES"

        gs_bucket = ""
        cram_path = ""
        crai_path = ""
        if file_path and len(file_path) > 5:
            gs_bucket = file_path.split('/')[2]
            cram_path = file_path

        if gs_bucket in ("tgg-rnaseq", ):
            continue

        # figure out crai path
        cram_and_crai_exist = False
        genome_version = project.genome_version
        if cram_path and (cram_path.endswith(".cram") or cram_path.endswith(".bam")):
            try:
                path_exists = hadoop_exists_using_cache(cram_path, double_check_if_cache_says_yes=True)
            except Exception as e:
                print(f"ERROR #{error_crams_counter}: while checking if file exists: {cram_path}: {e}. Skipping this cram...")
                error_crams_counter += 1
                path_exists = False # while requester-pays buckets aren't working

            if not path_exists:
                missing_crams_counter += 1
                print(f"{missing_crams_counter}: {sample_type} cram path not found: {cram_path}")
            else:
                if cram_path.endswith(".cram") and hadoop_exists_using_cache(re.sub(".cram$", ".cram.crai", cram_path), double_check_if_cache_says_yes=True):
                    crai_path = re.sub(".cram$", ".cram.crai", cram_path)
                elif cram_path.endswith(".cram") and hadoop_exists_using_cache(re.sub(".cram$", ".crai", cram_path), double_check_if_cache_says_yes=True):
                    crai_path = re.sub(".cram$", ".crai", cram_path)
                elif cram_path.endswith(".bam") and hadoop_exists_using_cache(re.sub(".bam$", ".bam.bai", cram_path), double_check_if_cache_says_yes=True):
                    crai_path = re.sub(".bam$", ".bam.bai", cram_path)
                elif cram_path.endswith(".bam") and hadoop_exists_using_cache(re.sub(".bam$", ".bai", cram_path), double_check_if_cache_says_yes=True):
                    crai_path = re.sub(".bam$", ".bai", cram_path)

                if crai_path:
                    cram_and_crai_exist = True
                else:
                    print(f"{missing_crams_counter}: {sample_type} crai path not found for cram: {cram_path}")

        # check for errors in seqr database where a WGS project contains a WES cram or vice versa.
        #if cram_and_crai_exist and (sample_type == "WGS" and any(k in cram_path.lower() for k in ("_wes_", "exome"))) \
        #    or (sample_type == "WES" and any(k in cram_path.lower() for k in ("_wgs_", "genome"))):
        #        print(f"ERROR: {sample_type} project contains cram with different sample type: {project.name}  {cram_path}. Skipping...")
        #        cram_path = crai_path = ""
        #        cram_and_crai_exist = False

        # validate genome_version
        if cram_and_crai_exist:
            true_genome_version = get_genome_version_from_bam_or_cram_header(cram_path)
            if str(true_genome_version) != str(project.genome_version):
                print(f"WARNING: {project.name} {project.guid} cram {cram_path} has genome_version = {true_genome_version} rather than {project.genome_version}")
                genome_version = true_genome_version
        else:
            missing_crams_counter += 1

        discovery_gene_ENSG_ids = []
        discovery_gene_symbols = []
        saved_variants = SavedVariant.objects.filter(family=family, varianttag__variant_tag_type__category='CMG Discovery Tags')
        for saved_variant in saved_variants:
            if not saved_variant.saved_variant_json or len(saved_variant.saved_variant_json) == 0:
                continue

            transcripts = [t for transcripts in saved_variant.saved_variant_json.get('transcripts', {}).values() for t in transcripts]
            if transcripts:
                transcripts.sort(key=lambda d: d.get("transcriptRank", 10000))
                if transcripts[0].get("geneId"):
                    discovery_gene_ENSG_ids.append(transcripts[0]["geneId"])
                if transcripts[0].get("geneSymbol"):
                    discovery_gene_symbols.append(transcripts[0]["geneSymbol"])

        f.write("\t".join(map(lambda s: re.sub("[\t\r\n]", " --- ", str(s)), [
            genome_version,
            project.guid,
            project.name,
            family.family_id,
            family.guid,
            individual.individual_id,
            individual.guid,
            individual.father.individual_id if individual.father else "",
            individual.father.guid if individual.father else "",
            individual.mother.individual_id if individual.mother else "",
            individual.mother.guid if individual.mother else "",
            individual.sex,
            AFFECTED_STATUS_LOOKUP[individual.affected],
            ANALYSIS_STATUS_LOOKUP[family.analysis_status],
            ", ".join(sorted({
                feature['id'] + ":" + hpo_name_map.get(feature['id'], "") for feature in (individual.features or [])
            })),  # aka. features
            ", ".join(sorted({
                feature['id'] + ":" + hpo_name_map.get(feature['id'], "") for feature in (individual.absent_features or [])
            })),  # aka. absent_features
            individual.disorders[0][0] if individual.disorders and len(individual.disorders[0]) > 0 else "",
            individual.notes or "",
            ", ".join(f"{k}: {float(v):0.1f}" for k, v in (individual.filter_flags or {}).items()),
            (
                    f"{family.assigned_analyst.username}:" + (
                    f"{family.assigned_analyst.first_name}:" if family.assigned_analyst.first_name else "") + (
                    f"{family.assigned_analyst.email}" if family.assigned_analyst.email else "")
            ) if family.assigned_analyst else "",
            individual.population or "",
            family.coded_phenotype or "",

            family.post_discovery_omim_number or "",
            ",".join(sorted(set(discovery_gene_ENSG_ids))),
            ",".join(sorted(set(discovery_gene_symbols))),

            gs_bucket,
            cram_path,
            crai_path,
            cram_and_crai_exist,
            sample_type,
        ])) + "\n")


# TODO get mean coverage, contamination, chimera rate from
#gs://seqr-datasets/v02/GRCh38/RDG_WES_Broad_Internal/v12/sample_qc/resources/callset_seq_metrics.txt
#gs://seqr-datasets/v02/GRCh38/RDG_WGS_Broad_Internal/v21/sample_qc/resources/callset_seq_metrics.txt,
#gs://seqr-datasets/v02/GRCh38/RDG_WGS_Broad_Internal/v22/sample_qc/resources/callset_seq_metrics.txt

print(f"{missing_crams_counter} cram/crai pairs missing")
print(f"{error_crams_counter} files skipped due to error while checking if file exists")

#%%
# read the table, sort it, and write it back out
output_path2 = re.sub("_unfinished$", "", output_path)
df = pd.read_table(output_path).sort_values(by=[
    'sample_type', 'genome_version', 'project_guid', 'individual_id',
])

df.to_csv(
    output_path2,
    index=False,
    header=True,
    sep="\t",
)

os.remove(output_path)
print("Done")
