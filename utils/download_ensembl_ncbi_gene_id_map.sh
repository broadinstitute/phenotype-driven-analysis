# Download the gene id mapping from https://www.genenames.org/download/custom/
# by selecting "NCBI Gene ID", "Ensembl gene ID", and also
# "NCBI Gene ID(supplied by NCBI)" and "Ensembl ID(supplied by Ensembl)"
# These

N=$(
    curl 'https://www.genenames.org/cgi-bin/download/custom?col=md_eg_id&col=md_ensembl_id&col=gd_pub_ensembl_id&col=gd_pub_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit' \
    | grep -v '^\s*$'  \
    | tee >(gzip -c - >  ensembl_ncbi_gene_id_map.tsv.gz) \
    | wc -l
 )
echo Downloaded $N Ensembl Id to NCBI Gene Id mappings from https://www.genenames.org/download/custom/
