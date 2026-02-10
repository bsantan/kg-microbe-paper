
all: setup_input_files download_kg rhea_chebi_competencies ec_competencies taxonomy_competencies

# Setup Input_Files directory with necessary data files
setup_input_files: download_ontologies extract_ncbitaxon restore_vital_file

# Download ontologies archive from kg-microbe releases
download_ontologies:
	@echo "Checking for ontologies.tar.gz..."
	@if [ ! -f "Input_Files/ontologies.tar.gz" ]; then \
		mkdir -p Input_Files && \
		echo "Downloading ontologies.tar.gz..." && \
		(wget -nc -P Input_Files https://github.com/Knowledge-Graph-Hub/kg-microbe/releases/download/2025-03-07/ontologies.tar.gz || \
		(echo "Error: Could not download ontologies.tar.gz" && exit 1)); \
	else \
		echo "ontologies.tar.gz already exists"; \
	fi

# Extract ncbitaxon_nodes.tsv from ontologies archive
extract_ncbitaxon:
	@echo "Checking for ncbitaxon_nodes.tsv..."
	@if [ ! -f "Input_Files/ncbitaxon_nodes.tsv" ]; then \
		if [ -f "Input_Files/ontologies.tar.gz" ]; then \
			echo "Extracting ncbitaxon_nodes.tsv from ontologies.tar.gz..." && \
			cd Input_Files && tar -xzf ontologies.tar.gz ./ncbitaxon_nodes.tsv && \
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
	@if [ ! -f "Input_Files/Vital_etal_butyrate+producing_microbes.csv" ]; then \
		mkdir -p Input_Files && \
		echo "Restoring Vital_etal_butyrate_producing_microbes.csv from git history..." && \
		(git show f760bf0:src/Input_Files/Vital_etal_butyrate_producing_microbes.csv > Input_Files/Vital_etal_butyrate_producing_microbes.csv 2>/dev/null && \
		cp Input_Files/Vital_etal_butyrate_producing_microbes.csv Input_Files/Vital_etal_butyrate+producing_microbes.csv && \
		echo "Vital et al. file restored") || \
		(echo "Error: Could not restore Vital et al. file from git history (commit f760bf0 not available)" && \
		 echo "This may fail in shallow clones or source archives. Consider sourcing from release artifacts." && \
		 exit 1); \
	else \
		echo "Vital et al. file already exists"; \
	fi

download_kg:
	mkdir -p Input_Files &&  wget -nc -P Input_Files https://portal.nersc.gov/project/m4689/KGMicrobe-biomedical-function-20250222.tar.gz  && mkdir -p Input_Files/kg-microbe-biomedical-function-cat && /usr/bin/tar  --strip-components=3 -xzf Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz -C Input_Files/kg-microbe-biomedical-function-cat

rhea_chebi_competencies:
	cd Input_Files/kg-microbe-biomedical-function-cat && awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} ($$1 ~ /NCBITaxon:/ && $$3 ~ /CHEBI:/) || ($$3 ~ /NCBITaxon:/ && $$1 ~ /CHEBI:/) || ($$3 ~ /NCBITaxon:/ && $$1 ~ /UniprotKB:/) || ($$1 ~ /RHEA:/ && $$3 ~ /CHEBI:/) || ($$1 ~ /UniprotKB:/ && $$3 ~ /RHEA:/) || ($$1 ~ /RHEA:/ && $$3 ~ /RHEA:/)' merged-kg_edges.tsv > merged-kg_edges_noEC.tsv

ec_competencies:
	cd Input_Files/kg-microbe-biomedical-function-cat && awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} ($$3 ~ /NCBITaxon:/ && $$1 ~ /UniprotKB:/) || ($$1 ~ /UniprotKB:/ && $$3 ~ /EC:2./) || ($$1 ~ /UniprotKB:/ && $$3 ~ /EC:1./) || ($$1 ~ /UniprotKB:/ && $$3 ~ /EC:4./) || ($$1 ~ /UniprotKB:/ && $$3 ~ /EC:5./)' merged-kg_edges.tsv > merged-kg_edges_competency_specific_ec.tsv

taxonomy_competencies:
	cd Input_Files/kg-microbe-biomedical-function-cat && awk -F '\t' 'BEGIN {print "subject\tpredicate\tobject"} ($$3 ~ /NCBITaxon:/ && $$1 ~ /NCBITaxon:/)' merged-kg_edges.tsv > merged-kg_edges_ncbitaxon.tsv

# Verify all required files are in place
verify:
	@echo "=========================================="
	@echo "Verifying Input Files Setup"
	@echo "=========================================="
	@echo ""
	@echo "Checking ncbitaxon_nodes.tsv:"
	@if [ -f "Input_Files/ncbitaxon_nodes.tsv" ]; then \
		ls -lh Input_Files/ncbitaxon_nodes.tsv | awk '{print "  ✓ File exists: " $$9 " (" $$5 ")"}'; \
		head -1 Input_Files/ncbitaxon_nodes.tsv | awk '{print "  ✓ Header: " substr($$0,1,80) "..."}'; \
	else \
		echo "  ✗ Missing: Input_Files/ncbitaxon_nodes.tsv"; \
	fi
	@echo ""
	@echo "Checking Vital et al. gold standard file:"
	@if [ -f "Input_Files/Vital_etal_butyrate+producing_microbes.csv" ]; then \
		ls -lh Input_Files/Vital_etal_butyrate+producing_microbes.csv | awk '{print "  ✓ File exists: " $$9 " (" $$5 ")"}'; \
	else \
		echo "  ✗ Missing: Input_Files/Vital_etal_butyrate+producing_microbes.csv"; \
	fi
	@echo ""
	@echo "Checking KG files:"
	@if [ -d "Input_Files/kg-microbe-biomedical-function-cat" ]; then \
		FILE_COUNT=$$(ls Input_Files/kg-microbe-biomedical-function-cat 2>/dev/null | wc -l | xargs); \
		echo "  ✓ Directory exists: Input_Files/kg-microbe-biomedical-function-cat ($$FILE_COUNT files)"; \
		for f in merged-kg_edges_ncbitaxon.tsv merged-kg_edges_noEC.tsv merged-kg_edges_competency_specific_ec.tsv merged-kg_nodes.tsv; do \
			if [ -f "Input_Files/kg-microbe-biomedical-function-cat/$$f" ]; then \
				ls -lh "Input_Files/kg-microbe-biomedical-function-cat/$$f" | awk '{print "    ✓ " $$9}'; \
			else \
				echo "    ✗ Missing: $$f"; \
			fi; \
		done; \
	else \
		echo "  ✗ Missing: Input_Files/kg-microbe-biomedical-function-cat"; \
	fi
	@echo ""
	@echo "=========================================="

# Clean up generated files (keeps downloaded archives)
clean:
	@echo "Cleaning generated competency files..."
	rm -f Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_noEC.tsv
	rm -f Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_competency_specific_ec.tsv
	rm -f Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_ncbitaxon.tsv
	@echo "Done"

# Clean everything including downloads (use with caution)
clean-all: clean
	@echo "Removing all downloaded files..."
	rm -rf Input_Files/
	@echo "Done"

# Display help information
help:
	@echo "=========================================="
	@echo "KG-Microbe Paper Analysis - Makefile Help"
	@echo "=========================================="
	@echo ""
	@echo "Main targets:"
	@echo "  make all              - Setup all input files and process KG data"
	@echo "  make setup_input_files - Download and setup required input files"
	@echo "  make verify           - Check that all required files are present"
	@echo ""
	@echo "Setup sub-targets:"
	@echo "  make download_ontologies - Download ontologies.tar.gz"
	@echo "  make extract_ncbitaxon   - Extract ncbitaxon_nodes.tsv from archive"
	@echo "  make restore_vital_file  - Restore Vital et al. gold standard file from git"
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

.PHONY: all setup_input_files download_ontologies extract_ncbitaxon restore_vital_file download_kg rhea_chebi_competencies ec_competencies taxonomy_competencies verify clean clean-all help

