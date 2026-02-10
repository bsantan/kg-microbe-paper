.PHONY: all setup_input_files download_ontologies extract_ncbitaxon restore_vital_file download_kg rhea_chebi_competencies ec_competencies taxonomy_competencies setup_gold_standard verify clean clean-all help

all: setup_input_files download_kg rhea_chebi_competencies ec_competencies taxonomy_competencies setup_gold_standard

# Setup Input_Files directory with necessary data files
setup_input_files: download_ontologies extract_ncbitaxon restore_vital_file

# Download ontologies archive from kg-microbe releases
download_ontologies:
	@echo "Checking for ontologies.tar.gz..."
	@if [ ! -f "data/Input_Files/ontologies.tar.gz" ]; then \
		mkdir -p data/Input_Files && \
		echo "Downloading ontologies.tar.gz..." && \
		(wget -nc -P data/Input_Files https://github.com/Knowledge-Graph-Hub/kg-microbe/releases/download/2025-03-07/ontologies.tar.gz || \
		(echo "Error: Could not download ontologies.tar.gz" && exit 1)); \
	else \
		echo "ontologies.tar.gz already exists"; \
	fi

# Extract ncbitaxon_nodes.tsv from ontologies archive
extract_ncbitaxon:
	@echo "Checking for ncbitaxon_nodes.tsv..."
	@if [ ! -f "data/Input_Files/ncbitaxon_nodes.tsv" ]; then \
		if [ -f "data/Input_Files/ontologies.tar.gz" ]; then \
			echo "Extracting ncbitaxon_nodes.tsv from ontologies.tar.gz..." && \
			cd data/Input_Files && tar -xzf ontologies.tar.gz ./ncbitaxon_nodes.tsv && \
			echo "Extracted ncbitaxon_nodes.tsv"; \
		else \
			echo "Warning: ontologies.tar.gz not found, cannot extract ncbitaxon_nodes.tsv"; \
			echo "Run 'make download_ontologies' first"; \
		fi; \
	else \
		echo "ncbitaxon_nodes.tsv already exists"; \
	fi

# Restore Vital et al. gold standard file from git history
# Note: This depends on git history being available (fails in shallow clones/archives)
restore_vital_file:
	@echo "Checking for Vital et al. gold standard file..."
	@if [ ! -f "data/Input_Files/Vital_etal_butyrate+producing_microbes.csv" ]; then \
		mkdir -p data/Input_Files && \
		echo "Restoring Vital_etal_butyrate_producing_microbes.csv from git history..." && \
		(git show f760bf0:src/Input_Files/Vital_etal_butyrate_producing_microbes.csv > data/Input_Files/Vital_etal_butyrate_producing_microbes.csv 2>/dev/null && \
		cp data/Input_Files/Vital_etal_butyrate_producing_microbes.csv data/Input_Files/Vital_etal_butyrate+producing_microbes.csv && \
		echo "Vital et al. file restored") || \
		(echo "Error: Could not restore Vital et al. file from git history (commit f760bf0 not available)" && \
		 echo "This may fail in shallow clones or source archives. Consider sourcing from release artifacts." && \
		 exit 1); \
	else \
		echo "Vital et al. file already exists"; \
	fi

# Download KG archive from NERSC
download_kg:
	@echo "Downloading knowledge graph from NERSC..."
	mkdir -p data/Input_Files && \
	wget -nc -P data/Input_Files https://portal.nersc.gov/project/m4689/KGMicrobe-biomedical-function-20250222.tar.gz && \
	mkdir -p data/Input_Files/kg-microbe-biomedical-function-cat && \
	tar --strip-components=3 -xzf data/Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz -C data/Input_Files/kg-microbe-biomedical-function-cat
	@echo "Knowledge graph downloaded and extracted successfully."

# Generate RHEA-CHEBI competency edge subset
rhea_chebi_competencies:
	@echo "Creating RHEA-CHEBI competency edge file..."
	cd data/Input_Files/kg-microbe-biomedical-function-cat && \
	awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} \
	{s=$$1; o=$$3; \
	 if ((index(s,"NCBITaxon:")==1 && index(o,"CHEBI:")==1) || \
	     (index(o,"NCBITaxon:")==1 && index(s,"CHEBI:")==1) || \
	     (index(o,"NCBITaxon:")==1 && index(s,"UniprotKB:")==1) || \
	     (index(s,"RHEA:")==1 && index(o,"CHEBI:")==1) || \
	     (index(s,"UniprotKB:")==1 && index(o,"RHEA:")==1) || \
	     (index(s,"RHEA:")==1 && index(o,"RHEA:")==1)) print}' merged-kg_edges.tsv > merged-kg_edges_noEC.tsv
	@echo "Created merged-kg_edges_noEC.tsv"

# Generate EC competency edge subset
ec_competencies:
	@echo "Creating EC competency edge file..."
	cd data/Input_Files/kg-microbe-biomedical-function-cat && \
	awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} \
	{s=$$1; o=$$3; \
	 if ((index(o,"NCBITaxon:")==1 && index(s,"UniprotKB:")==1) || \
	     (index(s,"UniprotKB:")==1 && (index(o,"EC:2.")==1 || index(o,"EC:1.")==1 || index(o,"EC:4.")==1 || index(o,"EC:5.")==1))) print}' merged-kg_edges.tsv > merged-kg_edges_competency_specific_ec.tsv
	@echo "Created merged-kg_edges_competency_specific_ec.tsv"

# Generate NCBITaxon taxonomy edge subset
taxonomy_competencies:
	@echo "Creating taxonomy competency edge file..."
	cd data/Input_Files/kg-microbe-biomedical-function-cat && \
	awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} \
	{if (index($$1,"NCBITaxon:")==1 && index($$3,"NCBITaxon:")==1) print}' merged-kg_edges.tsv > merged-kg_edges_ncbitaxon.tsv
	@echo "Created merged-kg_edges_ncbitaxon.tsv"

# Setup Gold Standard file for low-memory systems
setup_gold_standard:
	@echo "Setting up Gold Standard file for low-memory systems..."
	@if [ ! -f "data/Intermediate_Files_Competencies/butyrate_produces/Gold_Standard_Species_Overlap_butyrate_produces.csv" ]; then \
		mkdir -p data/Intermediate_Files_Competencies/butyrate_produces && \
		if [ -f "data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv" ]; then \
			cp data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv data/Intermediate_Files_Competencies/butyrate_produces/ && \
			echo "Copied Gold_Standard_Species_Overlap_butyrate_produces.csv to data/Intermediate_Files_Competencies/butyrate_produces/"; \
		else \
			echo "Warning: Source file data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv not found, skipping."; \
		fi; \
	else \
		echo "Gold Standard file already exists, skipping."; \
	fi

# Verify all required files are in place
verify:
	@echo "=========================================="
	@echo "Verifying Input Files Setup"
	@echo "=========================================="
	@echo ""
	@echo "Checking ncbitaxon_nodes.tsv:"
	@if [ -f "data/Input_Files/ncbitaxon_nodes.tsv" ]; then \
		ls -lh data/Input_Files/ncbitaxon_nodes.tsv | awk '{print "  ✓ File exists: " $$9 " (" $$5 ")"}'; \
		head -1 data/Input_Files/ncbitaxon_nodes.tsv | awk '{print "  ✓ Header: " substr($$0,1,80) "..."}'; \
	else \
		echo "  ✗ Missing: data/Input_Files/ncbitaxon_nodes.tsv"; \
	fi
	@echo ""
	@echo "Checking Vital et al. gold standard file:"
	@if [ -f "data/Input_Files/Vital_etal_butyrate+producing_microbes.csv" ]; then \
		ls -lh data/Input_Files/Vital_etal_butyrate+producing_microbes.csv | awk '{print "  ✓ File exists: " $$9 " (" $$5 ")"}'; \
	else \
		echo "  ✗ Missing: data/Input_Files/Vital_etal_butyrate+producing_microbes.csv"; \
	fi
	@echo ""
	@echo "Checking Gold Standard file:"
	@if [ -f "data/Intermediate_Files_Competencies/butyrate_produces/Gold_Standard_Species_Overlap_butyrate_produces.csv" ]; then \
		ls -lh data/Intermediate_Files_Competencies/butyrate_produces/Gold_Standard_Species_Overlap_butyrate_produces.csv | awk '{print "  ✓ File exists: " $$9 " (" $$5 ")"}'; \
	else \
		echo "  ✗ Missing: data/Intermediate_Files_Competencies/butyrate_produces/Gold_Standard_Species_Overlap_butyrate_produces.csv"; \
	fi
	@echo ""
	@echo "Checking KG files:"
	@if [ -d "data/Input_Files/kg-microbe-biomedical-function-cat" ]; then \
		FILE_COUNT=$$(ls data/Input_Files/kg-microbe-biomedical-function-cat 2>/dev/null | wc -l | xargs); \
		echo "  ✓ Directory exists: data/Input_Files/kg-microbe-biomedical-function-cat ($$FILE_COUNT files)"; \
		for f in merged-kg_edges_ncbitaxon.tsv merged-kg_edges_noEC.tsv merged-kg_edges_competency_specific_ec.tsv merged-kg_nodes.tsv merged-kg_edges.tsv; do \
			if [ -f "data/Input_Files/kg-microbe-biomedical-function-cat/$$f" ]; then \
				SIZE=$$(ls -lh "data/Input_Files/kg-microbe-biomedical-function-cat/$$f" | awk '{print $$5}'); \
				echo "    ✓ $$f ($$SIZE)"; \
			else \
				echo "    ✗ Missing: $$f"; \
			fi; \
		done; \
	else \
		echo "  ✗ Missing: data/Input_Files/kg-microbe-biomedical-function-cat"; \
	fi
	@echo ""
	@echo "=========================================="

# Clean up generated files (keeps downloaded archives)
clean:
	@echo "Cleaning generated competency files..."
	rm -f data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_noEC.tsv
	rm -f data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_competency_specific_ec.tsv
	rm -f data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_ncbitaxon.tsv
	@echo "Done"

# Clean everything including downloads (use with caution)
clean-all: clean
	@echo "Removing all downloaded files..."
	rm -rf data/Input_Files/
	rm -rf data/Intermediate_Files_Competencies/
	@echo "Done"

# Display help information
help:
	@echo "=========================================="
	@echo "KG-Microbe Paper Analysis - Makefile Help"
	@echo "=========================================="
	@echo ""
	@echo "Main targets:"
	@echo "  make all                - Setup all input files and process KG data"
	@echo "  make setup_input_files  - Download and setup required input files"
	@echo "  make verify             - Check that all required files are present"
	@echo ""
	@echo "Setup sub-targets:"
	@echo "  make download_ontologies  - Download ontologies.tar.gz"
	@echo "  make extract_ncbitaxon    - Extract ncbitaxon_nodes.tsv from archive"
	@echo "  make restore_vital_file   - Restore Vital et al. gold standard file from git"
	@echo "  make setup_gold_standard  - Setup gold standard file for analysis"
	@echo ""
	@echo "KG processing targets:"
	@echo "  make download_kg              - Download KG archive from NERSC"
	@echo "  make rhea_chebi_competencies  - Generate RHEA-CHEBI edge subset"
	@echo "  make ec_competencies          - Generate EC edge subset"
	@echo "  make taxonomy_competencies    - Generate NCBITaxon edge subset"
	@echo ""
	@echo "Cleanup targets:"
	@echo "  make clean      - Remove generated files (keeps downloads)"
	@echo "  make clean-all  - Remove all files including downloads"
	@echo ""
	@echo "Usage examples:"
	@echo "  make              # Run all setup and processing"
	@echo "  make verify       # Check current setup status"
	@echo "  make help         # Show this help message"
	@echo ""
