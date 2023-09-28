from seqr.models import IgvSample


def impute_sample_type(individual):
    file_paths = []
    for igv_sample in IgvSample.objects.filter(individual=individual).distinct('file_path'):
        if igv_sample.file_path and igv_sample.file_path.startswith("gs://") \
            and not igv_sample.file_path.endswith(".bed.gz") \
            and not igv_sample.file_path.endswith(".bigWig"):  # ignore .bed.gz gCNV files and .bigWig RNA-seq files
            file_paths.append(igv_sample.file_path)

    if len(file_paths) > 1:
        print(f"Multiple file paths for {individual.guid}: {individual.individual_id}:")
        print(f"   " + "\n".join(file_paths))

    file_path = None
    if len(file_paths) > 0:
        file_path = file_paths[0]

    family = individual.family
    project = family.project

    project_guid_lowercase = project.guid.lower()
    project_name_lowercase = project.name.lower().replace(' ', '_')

    if 'genome' in project_guid_lowercase or 'genome' in project_name_lowercase or '_wgs' in project_guid_lowercase or '_wgs' in project_name_lowercase:
        # and "_wes_" not in file_path_lowercase and "_wes/" not in file_path_lowercase and "exome" not in file_path_lowercase
        sample_type = 'WGS'
    elif 'exome' in project_guid_lowercase or 'exome' in project_name_lowercase or '_wes' in project_guid_lowercase or '_wes' in project_name_lowercase or project.guid in (
            "R0294_myoseq_v20", "R0208_kang_v11", "R0230_bonnemann_jan2016_107s", "R0418_manton_qi_project",
            "R0540_mgh_pathways_probands_on", "R0310_manton_pcrfree_4s", "R0261_bonnemann_v9"):
        # or "_wes_" in file_path_lowercase or "_wes/" in file_path_lowercase or "exome" in file_path_lowercase
        sample_type = 'WES'
    else:
        #  try to infer sample type from adjacent samples
        adjacent_samples = individual.sample_set.all()
        sample_types = set([sample.sample_type for sample in adjacent_samples if sample.sample_type != "RNA"])

        try:
            assert len(sample_types) == 1, (
                "sample types", sample_types,
                "project guid", individual.family.project.guid,
                "project name", individual.family.project.name,
                "individual id", individual.individual_id,
            )
            sample_type = next(iter(sample_types))
        except Exception as e:
            print("Unable to infer sample type. ", e)
            return None, file_path

    return sample_type, file_path
