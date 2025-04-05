
all: download_kg rhea_chebi_competencies ec_competencies taxonomy_competencies

download_kg:
	mkdir -p Input_Files &&  wget -nc -P Input_Files https://kghub.io/kg-microbe/KGMicrobe-biomedical-function-20250222.tar.gz  && mkdir -p Input_Files/kg-microbe-biomedical-function-cat && /usr/bin/tar  --strip-components=3 -xzf Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz -C Input_Files/kg-microbe-biomedical-function-cat

rhea_chebi_competencies:
	cd Input_Files/kg-microbe-biomedical-function-cat && awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} ($$1 ~ /NCBITaxon:/ && $$3 ~ /CHEBI:/) || ($$3 ~ /NCBITaxon:/ && $$1 ~ /CHEBI:/) || ($$3 ~ /NCBITaxon:/ && $$1 ~ /UniprotKB:/) || ($$1 ~ /RHEA:/ && $$3 ~ /CHEBI:/) || ($$1 ~ /UniprotKB:/ && $$3 ~ /RHEA:/) || ($$1 ~ /RHEA:/ && $$3 ~ /RHEA:/)' merged-kg_edges.tsv > merged-kg_edges_noEC.tsv

ec_competencies:
	cd Input_Files/kg-microbe-biomedical-function-cat && awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} ($$3 ~ /NCBITaxon:/ && $$1 ~ /UniprotKB:/) || ($$1 ~ /UniprotKB:/ && $$3 ~ /EC:2./) || ($$1 ~ /UniprotKB:/ && $$3 ~ /EC:1./) || ($$1 ~ /UniprotKB:/ && $$3 ~ /EC:4./) || ($$1 ~ /UniprotKB:/ && $$3 ~ /EC:5./)' merged-kg_edges.tsv > merged-kg_edges_competency_specific_ec.tsv

taxonomy_competencies:
	cd Input_Files/kg-microbe-biomedical-function-cat && awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} ($$3 ~ /NCBITaxon:/ && $$1 ~ /NCBITaxon:/)' merged-kg_edges.tsv > merged-kg_edges_ncbitaxon.tsv

