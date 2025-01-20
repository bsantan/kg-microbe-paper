
rhea_chebi_competencies:
    awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} ($1 ~ /NCBITaxon:/ && $3 ~ /CHEBI:/) || ($3 ~ /NCBITaxon:/ && $1 ~ /CHEBI:/) || ($3 ~ /NCBITaxon:/ && $1 ~ /UniprotKB:/) || ($1 ~ /RHEA:/ && $3 ~ /CHEBI:/) || ($1 ~ /UniprotKB:/ && $3 ~ /RHEA:/)' merged-kg_edges.tsv > merged-kg_edges_noEC.tsv


ec_competencies:
    awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} ($3 ~ /NCBITaxon:/ && $1 ~ /UniprotKB:/) || ($1 ~ /UniprotKB:/ && $3 ~ /EC:/)' merged-kg_edges.tsv > merged-kg_edges_competency_only_ec.tsv

taxonomy_competencies:
    awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} ($3 ~ /NCBITaxon:/ && $1 ~ /NCBITaxon:/)' merged-kg_edges.tsv > merged-kg_edges_ncbitaxon.tsv

include kg-microbe-paper.Makefile