
# Read in edges and nodes files
create or replace table edges as 
select 
    subject,
    predicate,
    object,
from read_csv("/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Input_Files/merged-kg/merged-kg_edges.tsv", filename=true, union_by_name=true);

create or replace table nodes as
select
    id,
    name
from read_csv("/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Input_Files/merged-kg/merged-kg_nodes.tsv", filename=true, union_by_name=true);


## To query organismal traits by species and strain
# Get all bugs that have organismal trait of produces/consumes metabolite
## NCBITaxon - produces/consumes - metabolite
## species_metabolite_edges
CREATE OR REPLACE TABLE species_metabolite_edges AS
SELECT * FROM edges e
WHERE (subject LIKE 'NCBITaxon%' OR subject LIKE 'strain%')
  AND object = 'CHEBI:17968'; -- 'butyrate' 'tryptophan' 'CHEBI:17968'

COPY (SELECT subject, object FROM species_metabolite_edges) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_organismal_traits.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

# Instead, replace the strains with their NCBITaxon associated species
## strain - produces/consumes - metabolite
## strain_to_ncbi
CREATE OR REPLACE TABLE strain_to_ncbi AS
SELECT e.subject AS strain_subject,
       e.object AS ncbi_object
FROM edges e
WHERE split_part(e.subject, ':', 1) = 'strain'
  AND e.predicate = 'biolink:subclass_of'
  AND split_part(e.object, ':', 1) = 'NCBITaxon';

-- Step 2: Create the updated strain_metabolite_edges table
CREATE OR REPLACE TABLE strain_metabolite_edges AS
SELECT e.*, 
       split_part(subject, ':', 1) AS subject_prefix,
       COALESCE(s2n.ncbi_object, e.subject) AS updated_subject
FROM edges e
LEFT JOIN strain_to_ncbi s2n ON e.subject = s2n.strain_subject
WHERE (split_part(subject, ':', 1) = 'NCBITaxon' OR split_part(subject, ':', 1) = 'strain')
  AND object = 'CHEBI:17968';

-- Step 3: Export the updated results to a TSV file
COPY metabolite_edges TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_organismal_traits_strains.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

--------------------
#### To query genomic traits through GO
## UniprotKB - derives_from - Proteome - derives_from - NCBITaxon
## uniprot_to_ncbi
-- Step 1: Create the uniprot_to_ncbi table
CREATE OR REPLACE TABLE uniprot_to_ncbi AS
SELECT e1.subject AS uniprotkb, e2.object AS ncbi_taxon
FROM edges e1
JOIN edges e2
  ON e1.object = e2.subject
WHERE split_part(e1.subject, ':', 1) = 'UniprotKB'
  AND e1.predicate = 'biolink:derives_from'
  AND split_part(e2.object, ':', 1) = 'NCBITaxon'
  AND e2.predicate = 'biolink:derives_from';

-- Step 2: Perform the main query and include the replacement of the subject column, Use # NCBITaxon if available, else original subject and also don't include uniprotkb without proteome
#! Including statement to get rid of proteins without proteome, which is bug in this build only- NOT NULL
## NCBITaxon - participates_in - GO
## final_result
CREATE OR REPLACE TABLE ncbi_to_go AS
SELECT u.ncbi_taxon AS subject,
        e.predicate,
        e.object
FROM edges e
LEFT JOIN uniprot_to_ncbi u ON e.subject = u.uniprotkb
WHERE split_part(subject, ':', 1) = 'UniprotKB'
  AND (e.object = 'GO:0046358' OR e.object = 'GO:0046359' OR e.object = 'GO:0047761' OR e.object = 'GO:1903544' OR e.object = 'GO:0019605' OR e.object = 'GO:1903545' OR e.object = 'GO:0047371' OR e.object = 'GO:0030645' OR e.object = 'GO:0033508' OR e.object = 'GO:0044581' OR e.object = 'GO:0052638' OR e.object = 'GO:0019672' OR e.object = 'GO:1900500' OR e.object = 'GO:1900502' OR e.object = 'GO:1900501' OR e.object = 'GO:0034338' OR e.object = 'GO:1900751')
  AND u.ncbi_taxon IS NOT NULL;
    --tryptophan: object = 'GO:0009034' OR object = 'GO:0006569' OR object = 'GO:0006568' OR object = 'GO:0036469' or object = 'GO:0000162' OR e.object = 'GO:0004830' OR e.object = 'GO:0006436');
  --butyrate: object = 'GO:0046358' OR object = 'GO:0046359' OR object = 'GO:0047761' OR object = 'GO:1903544' OR object = 'GO:0019605' OR object = 'GO:1903545' OR object = 'GO:0047371' OR object = 'GO:0030645' OR object = 'GO:0033508' OR object = 'GO:0044581' OR object = 'GO:0052638' OR object = 'GO:0019672' OR object = 'GO:1900500' OR object = 'GO:1900502' OR object = 'GO:1900501' OR object = 'GO:0034338' OR object = 'GO:1900751'

-- Step 3: Export the ncbi_to_go to a TSV file
COPY (SELECT subject, object FROM ncbi_to_go) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_genomic_traits_GO.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

# Get all NCBITaxon IDs in organismal trait species and genomic results
## NCBITaxon - metabolite
## matching_species_organismal_go
CREATE OR REPLACE TABLE matching_species_organismal_go AS
SELECT sme.subject AS subject_id, sme.object AS metabolite_object, ncbi_to_go.subject
FROM species_metabolite_edges sme
JOIN ncbi_to_go ncbi_to_go ON sme.subject = ncbi_to_go.subject;

-- Step 3: Export the matching rows to a TSV file
COPY (SELECT subject_id, metabolite_object FROM matching_species_organismal_go) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_organismal_genomic_go_comparison_species.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

# Get all NCBITaxon IDs in organismal trait species and genomic results, updated_subject column is associated ncbitaxon
## NCBITaxon - metabolite
## matching_strain_organismal_go
CREATE OR REPLACE TABLE matching_strain_organismal_go AS
SELECT sme.updated_subject AS subject_id, sme.object AS metabolite_object, ncbi_to_go.subject
FROM strain_metabolite_edges sme
JOIN ncbi_to_go ncbi_to_go ON sme.subject = ncbi_to_go.subject;

-- Step 3: Export the matching rows to a TSV file
COPY (SELECT subject_id, metabolite_object FROM matching_strain_organismal_go) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_organismal_genomic_go_comparison_strain.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

-------------------

#### To query genomic traits through RHEA
# Create table of RHEA to CHEBI relevant terms
## rhea - has_input/output - chebi
## rhea_to_chebi
-- Step 1: Create the rhea_to_chebi table
CREATE OR REPLACE TABLE rhea_to_chebi AS
SELECT *
FROM edges e
WHERE split_part(e.subject, ':', 1) = 'RHEA'
  AND (e.object = 'CHEBI:17968' OR e.object = 'CHEBI:85867' OR e.object = 'CHEBI:88764' OR e.object = 'CHEBI:87684' OR e.object = 'CHEBI:64103' OR e.object = 'CHEBI:88806' OR e.object = 'CHEBI:90150' OR e.object = 'CHEBI:87683' OR e.object = 'CHEBI:89719' OR e.object = 'CHEBI:87318' OR e.object = 'CHEBI:453' OR e.object = 'CHEBI:32097' OR e.object = 'CHEBI:50477' OR e.object = 'CHEBI:87422' OR e.object = 'CHEBI:31415');
  ---  butyrate: e.object = 'CHEBI:85867'
  --- tryptophan: e.object = 'CHEBI:27897' OR e.object = 'CHEBI:16296' OR e.object = 'CHEBI:16828' OR e.object = 'CHEBI:57719'

# Create table of NCBITaxon nodes that have edges to those Chebi results
## ncbitaxon - rhea
## ncbitaxon_to_rhea
#! Including statement to get rid of proteins without proteome, which is bug in this build only- NOT NULL
-- Step 1: Create the ncbitaxon_to_rhea table
CREATE OR REPLACE TABLE ncbitaxon_to_rhea AS
SELECT 
  COALESCE(u.ncbi_taxon, e.subject) AS subject, 
  e.predicate, 
  e.object
FROM edges e
LEFT JOIN uniprot_to_ncbi u
  ON e.subject = u.uniprotkb
WHERE split_part(e.subject, ':', 1) = 'UniprotKB'
  AND split_part(e.object, ':', 1) = 'RHEA'
  AND u.ncbi_taxon IS NOT NULL; 


# Create table of ncbitaxon to Chebi
## ncbitaxon - chebi
## ncbitaxon_to_chebi
CREATE OR REPLACE TABLE ncbitaxon_to_chebi AS
SELECT ncbitaxon_to_rhea.subject AS ncbitaxon, 
       rhea_to_chebi.object AS chebi
FROM ncbitaxon_to_rhea
JOIN rhea_to_chebi ON ncbitaxon_to_rhea.object = rhea_to_chebi.subject;

COPY (SELECT ncbitaxon, chebi FROM ncbitaxon_to_chebi) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_genomic_traits_RHEA.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

# Get all NCBITaxon IDs in organismal trait species and RHEA genomic results
## NCBITaxon - metabolite
## matching_species_organismal_rhea
CREATE OR REPLACE TABLE matching_species_organismal_rhea AS
SELECT sme.subject AS subject_id, sme.object AS metabolite_object, ncbitaxon_to_chebi.ncbitaxon
FROM strain_metabolite_edges sme
JOIN ncbitaxon_to_chebi ON sme.subject = ncbitaxon_to_chebi.ncbitaxon;

-- Step 3: Export the matching rows to a TSV file
COPY (SELECT subject_id, metabolite_object FROM matching_species_organismal_rhea) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_organismal_genomic_rhea_comparison_species.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

# Get all NCBITaxon IDs in organismal trait species and RHEA genomic results, updated_subject column is associated ncbitaxon
## NCBITaxon - metabolite
## matching_strain_organismal_rhea
CREATE OR REPLACE TABLE matching_strain_organismal_rhea AS
SELECT sme.updated_subject AS subject_id, sme.object AS metabolite_object, ncbitaxon_to_chebi.ncbitaxon
FROM strain_metabolite_edges sme
JOIN ncbitaxon_to_chebi ON sme.subject = ncbitaxon_to_chebi.ncbitaxon;

-- Step 3: Export the matching rows to a TSV file
COPY (SELECT subject_id, metabolite_object FROM matching_strain_organismal_rhea) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_organismal_genomic_rhea_comparison_strain.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

-------------------