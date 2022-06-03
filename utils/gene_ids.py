import os
import pandas as pd


def get_entrez_to_ensembl_id_map():
    """Read ensembl_ncbi_gene_id_map.tsv.gz and return dict. mapping NCBI gene id to Ensembl gene id"""

    ensembl_id_column = "Ensembl ID(supplied by Ensembl)"
    ensembl_id_column2 = "Ensembl gene ID"
    ncbi_id_column = "NCBI Gene ID(supplied by NCBI)"
    ncbi_id_column2 = "NCBI Gene ID"

    mapping_tsv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ensembl_ncbi_gene_id_map.tsv.gz")
    ncbi_to_ensembl_gene_ids_df = pd.read_table(mapping_tsv, dtype=str)
    ncbi_to_ensembl_gene_ids_df = ncbi_to_ensembl_gene_ids_df.fillna("")

    ncbi_to_ensembl_gene_ids = {}
    for _, row in ncbi_to_ensembl_gene_ids_df.iterrows():
        if row[ensembl_id_column] and row[ncbi_id_column]:
            ncbi_to_ensembl_gene_ids[str(int(row[ncbi_id_column]))] = row[ensembl_id_column]
        if row[ensembl_id_column2] and row[ncbi_id_column2]:
            ncbi_to_ensembl_gene_ids[str(int(row[ncbi_id_column2]))] = row[ensembl_id_column2]
    return ncbi_to_ensembl_gene_ids
