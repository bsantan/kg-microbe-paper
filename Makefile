.PHONY: all download_kg download_ontologies rhea_chebi_competencies ec_competencies taxonomy_competencies setup_gold_standard clean

all: download_kg download_ontologies rhea_chebi_competencies ec_competencies taxonomy_competencies

download_kg:
	@echo "Downloading knowledge graph from NERSC..."
	mkdir -p src/Input_Files && \
	wget -nc -P src/Input_Files https://portal.nersc.gov/project/m4689/KGMicrobe-biomedical-function-20250222.tar.gz && \
	mkdir -p src/Input_Files/kg-microbe-biomedical-function-cat && \
	tar --strip-components=3 -xzf src/Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz -C src/Input_Files/kg-microbe-biomedical-function-cat
	@echo "Knowledge graph downloaded and extracted successfully."

download_ontologies:
	@echo "Downloading NCBI Taxonomy ontologies..."
	wget -nc -P src/Input_Files https://github.com/Knowledge-Graph-Hub/kg-microbe/releases/download/2025-03-07/ontologies.tar.gz
	@echo "Extracting ncbitaxon_nodes.tsv..."
	cd src/Input_Files && \
	mkdir -p tmp_ontologies && \
	tar -xzf ontologies.tar.gz -C tmp_ontologies && \
	find tmp_ontologies -name 'ncbitaxon_nodes.tsv' -exec mv {} . \; && \
	rm -rf tmp_ontologies
	@echo "Ontologies downloaded and extracted successfully."

rhea_chebi_competencies:
	@echo "Creating RHEA-CHEBI competency edge file..."
	cd src/Input_Files/kg-microbe-biomedical-function-cat && \
	awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} \
	{s=$$1; o=$$3; \
	 if ((index(s,"NCBITaxon:")==1 && index(o,"CHEBI:")==1) || \
	     (index(o,"NCBITaxon:")==1 && index(s,"CHEBI:")==1) || \
	     (index(o,"NCBITaxon:")==1 && index(s,"UniprotKB:")==1) || \
	     (index(s,"RHEA:")==1 && index(o,"CHEBI:")==1) || \
	     (index(s,"UniprotKB:")==1 && index(o,"RHEA:")==1) || \
	     (index(s,"RHEA:")==1 && index(o,"RHEA:")==1)) print}' merged-kg_edges.tsv > merged-kg_edges_noEC.tsv
	@echo "Created merged-kg_edges_noEC.tsv"

ec_competencies:
	@echo "Creating EC competency edge file..."
	cd src/Input_Files/kg-microbe-biomedical-function-cat && \
	awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} \
	{s=$$1; o=$$3; \
	 if ((index(o,"NCBITaxon:")==1 && index(s,"UniprotKB:")==1) || \
	     (index(s,"UniprotKB:")==1 && (index(o,"EC:2.")==1 || index(o,"EC:1.")==1 || index(o,"EC:4.")==1 || index(o,"EC:5.")==1))) print}' merged-kg_edges.tsv > merged-kg_edges_competency_specific_ec.tsv
	@echo "Created merged-kg_edges_competency_specific_ec.tsv"

taxonomy_competencies:
	@echo "Creating taxonomy competency edge file..."
	cd src/Input_Files/kg-microbe-biomedical-function-cat && \
	awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} \
	{if (index($$1,"NCBITaxon:")==1 && index($$3,"NCBITaxon:")==1) print}' merged-kg_edges.tsv > merged-kg_edges_ncbitaxon.tsv
	@echo "Created merged-kg_edges_ncbitaxon.tsv"

setup_gold_standard:
	@echo "Setting up Gold Standard file for low-memory systems..."
	@if [ ! -f "src/Intermediate_Files_Competencies/butyrate_produces/Gold_Standard_Species_Overlap_butyrate_produces.csv" ]; then \
		mkdir -p src/Intermediate_Files_Competencies/butyrate_produces && \
		cp data/Gold_Standard_Species_Overlap_butyrate_produces.csv src/Intermediate_Files_Competencies/butyrate_produces/ && \
		echo "Copied Gold_Standard_Species_Overlap_butyrate_produces.csv to src/Intermediate_Files_Competencies/butyrate_produces/"; \
	else \
		echo "Gold Standard file already exists, skipping."; \
	fi

clean:
	@echo "Cleaning downloaded files..."
	rm -rf src/Input_Files/kg-microbe-biomedical-function-cat
	rm -rf src/Input_Files/hmp_supplementary
	rm -f src/Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz
	rm -f src/Input_Files/ontologies.tar.gz
	rm -f src/Input_Files/ncbitaxon_nodes.tsv
	rm -rf Input_Files
	@echo "Clean complete."
