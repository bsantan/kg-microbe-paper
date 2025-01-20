
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

# Get all bugs that have organismal trait of produces/consumes metabolite
CREATE OR REPLACE TABLE metabolite_edges AS
SELECT * FROM edges e
WHERE (subject LIKE 'NCBITaxon%' OR subject LIKE 'strain%')
  AND object = 'CHEBI:17968'; -- 'butyrate' 'tryptophan' 'CHEBI:17968'

COPY (SELECT subject, object FROM metabolite_edges) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_organismal_traits.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

# Instead, replace the strains with their NCBITaxon associated species
CREATE OR REPLACE TABLE strain_to_ncbi AS
SELECT e.subject AS strain_subject,
       e.object AS ncbi_object
FROM edges e
WHERE split_part(e.subject, ':', 1) = 'strain'
  AND e.predicate = 'biolink:subclass_of'
  AND split_part(e.object, ':', 1) = 'NCBITaxon';

-- Step 2: Create the updated metabolite_edges table
CREATE OR REPLACE TABLE metabolite_edges AS
SELECT e.*, 
       split_part(subject, ':', 1) AS subject_prefix,
       COALESCE(s2n.ncbi_object, e.subject) AS updated_subject
FROM edges e
LEFT JOIN strain_to_ncbi s2n ON e.subject = s2n.strain_subject
WHERE (split_part(subject, ':', 1) = 'NCBITaxon' OR split_part(subject, ':', 1) = 'strain')
  AND object = 'CHEBI:17968';

-- Step 3: Export the updated results to a TSV file
COPY metabolite_edges TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_organismal_traits_species.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);

-- (SELECT updated_subject, object_name FROM metabolite_edges) 

# Get all uniprot-bug pairs
CREATE OR REPLACE TABLE uniprot_ncbi_edges AS
SELECT e.*,
       split_part(subject, ':', 1) AS subject_prefix,
       split_part(object, ':', 1) AS object_prefix,
       predicate
FROM edges e
WHERE split_part(subject, ':', 1) = 'UniprotKB'
  AND predicate = 'biolink:derives_from'
  AND split_part(object, ':', 1) = 'NCBITaxon';

# Look at all uniprot-GO edges
SELECT e.*, 
       split_part(subject, ':', 1) AS subject_prefix,
       split_part(object, ':', 1) AS object_prefix,
       predicate
FROM edges e
WHERE split_part(subject, ':', 1) = 'UniprotKB'
  AND (object = 'GO:0046358' OR object = 'GO:0046359' OR object = 'GO:0047761' OR object = 'GO:1903544' OR object = 'GO:0019605' OR object = 'GO:1903545' OR object = 'GO:0047371' OR object = 'GO:0030645' OR object = 'GO:0033508' OR object = 'GO:0044581' OR object = 'GO:0052638' OR object = 'GO:0019672' OR object = 'GO:1900500' OR object = 'GO:1900502' OR object = 'GO:1900501' OR object = 'GO:0034338' OR object = 'GO:1900751');
  --tryptophan: object = 'GO:0009034' OR object = 'GO:0006569' OR object = 'GO:0006568' OR object = 'GO:0036469' or object = 'GO:0000162' OR e.object = 'GO:0004830' OR e.object = 'GO:0006436');
  --butyrate: object = 'GO:0046358' OR object = 'GO:0046359' OR object = 'GO:0047761' OR object = 'GO:1903544' OR object = 'GO:0019605' OR object = 'GO:1903545' OR object = 'GO:0047371' OR object = 'GO:0030645' OR object = 'GO:0033508' OR object = 'GO:0044581' OR object = 'GO:0052638' OR object = 'GO:0019672' OR object = 'GO:1900500' OR object = 'GO:1900502' OR object = 'GO:1900501' OR object = 'GO:0034338' OR object = 'GO:1900751'


# Create table of NCBITaxon to GO edges, according to Uniprot associations
-- Step 1: Create the uniprot_to_ncbi table
CREATE OR REPLACE TABLE uniprot_to_ncbi AS
SELECT e.subject AS uniprotkb_subject,
       e.object AS ncbi_object
FROM edges e
WHERE split_part(e.subject, ':', 1) = 'UniprotKB'
  AND e.predicate = 'biolink:derives_from'
  AND split_part(e.object, ':', 1) = 'NCBITaxon';

-- Step 2: Perform the main query and include the replacement of the subject column, Use NCBITaxon if available, else original subject
CREATE OR REPLACE TABLE final_result AS
SELECT e.*, 
       split_part(subject, ':', 1) AS subject_prefix,
       split_part(object, ':', 1) AS object_prefix,
       COALESCE(u.ncbi_object, e.subject) AS new_subject,
       predicate
FROM edges e
LEFT JOIN uniprot_to_ncbi u ON e.subject = u.uniprotkb_subject
WHERE split_part(subject, ':', 1) = 'UniprotKB'
  AND (e.object = 'GO:0046358' OR e.object = 'GO:0046359' OR e.object = 'GO:0047761' OR e.object = 'GO:1903544' OR e.object = 'GO:0019605' OR e.object = 'GO:1903545' OR e.object = 'GO:0047371' OR e.object = 'GO:0030645' OR e.object = 'GO:0033508' OR e.object = 'GO:0044581' OR e.object = 'GO:0052638' OR e.object = 'GO:0019672' OR e.object = 'GO:1900500' OR e.object = 'GO:1900502' OR e.object = 'GO:1900501' OR e.object = 'GO:0034338' OR e.object = 'GO:1900751'); --GO:0047761

-- Step 3: Export the final results to a TSV file
COPY (SELECT new_subject, object FROM final_result) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_genomic_traits_GO.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);


# Get all NCBITaxon IDs in organismal trait strains and genomic results
# CREATE OR REPLACE TABLE overlapping_taxa AS
# SELECT be.subject AS ncbi_subject_butyrate,
#        fr.new_subject AS ncbi_subject_final_result,
#        fr.object
# FROM metabolite_edges be
# JOIN final_result fr ON be.subject = fr.new_subject;

# COPY overlapping_taxa TO 'NCBI_in_both_trp.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);



# Get all NCBITaxon IDs in organismal trait species and genomic results
CREATE OR REPLACE TABLE matching_subjects AS
SELECT me.updated_subject AS subject_id, me.object AS metabolite_object, fr.new_subject
FROM metabolite_edges me
JOIN final_result fr ON me.updated_subject = fr.new_subject;

-- Step 3: Export the matching rows to a TSV file
COPY (SELECT subject_id, metabolite_object FROM matching_subjects) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_organismal_genomic_comparison_species.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);


# Get all NCBITaxon IDs with Uniprot proteome
CREATE OR REPLACE TABLE ncbi_to_proteome AS
SELECT *
FROM edges
WHERE object LIKE 'NCBITaxon:%'
  AND subject LIKE 'Proteomes:%'
  AND predicate = 'biolink:derives_from';

COPY (SELECT * FROM ncbi_to_proteome) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/Uniprot_proteome_NCBITaxa.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);


# Get all NCBITaxon IDs in organismal trait species and uniprot
CREATE OR REPLACE TABLE matching_proteomes AS
SELECT metabolite_edges.updated_subject AS subject_id, 
       metabolite_edges.object_name AS metabolite_object, 
       ncbi_to_proteome.object
FROM metabolite_edges
JOIN ncbi_to_proteome ON metabolite_edges.updated_subject = ncbi_to_proteome.object;

-- Step 3: Export the matching rows to a TSV file
COPY (SELECT subject_id FROM matching_proteomes) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/Uniprot_proteome_organismal_comparison_species.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);


#### To query tryptophan trait via RHEA
# Create table of RHEA to CHEBI trp relevant terms
-- Step 1: Create the rhea_to_chebi table
CREATE OR REPLACE TABLE rhea_to_chebi AS
SELECT e.subject AS rhea_subject,
       e.object AS chebi_object
FROM edges e
-- JOIN nodes obj ON e.object = obj.id
WHERE split_part(e.subject, ':', 1) = 'RHEA'
  AND (e.object = "CHEBI:27897" OR e.object = "CHEBI:16296" OR e.object = "CHEBI:16828" OR e.object = "CHEBI:57719");

# Create table of NCBITaxon nodes that have edges to those Chebi results
 # Create table of RHEA to CHEBI trp relevant terms
-- Step 1: Create the uniprot_to_rhea table
CREATE OR REPLACE TABLE uniprot_to_rhea AS
SELECT e.subject AS uniprot_subject,
       e.object AS rhea_object
FROM edges e
JOIN nodes obj ON e.object = obj.id
WHERE split_part(e.subject, ':', 1) = 'UniprotKB'
  AND split_part(e.object, ':', 1) = 'Rhea';

# Create table of Uniprot to Rhea
CREATE OR REPLACE TABLE uniprot_to_chebi AS
SELECT uniprot_to_rhea.subject AS subject_id, 
       uniprot_to_rhea.object AS object_id, 
       rhea_to_chebi.subject
FROM uniprot_to_rhea
JOIN rhea_to_chebi ON uniprot_to_rhea.object = rhea_to_chebi.subject;

# Create table of NCBItaxon to Uniprot



-- Step 2: Perform the main query and include the replacement of the subject column
CREATE OR REPLACE TABLE final_result AS
SELECT e.*, 
       COALESCE(obj.name, obj.id) AS object_name,
       COALESCE(sub.name, sub.id) AS subject_name,
       split_part(sub.id, ':', 1) AS subject_prefix,
       split_part(obj.id, ':', 1) AS object_prefix,
       COALESCE(u.ncbi_object, e.subject) AS new_subject, -- Use NCBITaxon if available, else original subject
       predicate
FROM edges e
LEFT JOIN nodes obj ON e.object = obj.id
LEFT JOIN nodes sub ON e.subject = sub.id
LEFT JOIN uniprot_to_chebi u ON e.subject = u.uniprotkb_subject
WHERE split_part(sub.id, ':', 1) = 'UniprotKB'
  AND (e.object = "CHEBI:27897" OR e.object = "CHEBI:16296" OR e.object = "CHEBI:16828" OR e.object = "CHEBI:57719");
  
  -- (e.object = 'GO:0009034' OR e.object = 'GO:0006569' OR e.object = 'GO:0006568' OR e.object = 'GO:0036469' or e.object = 'GO:0000162' OR e.object = 'GO:0004834' OR e.object = 'GO:0004830' OR e.object = 'GO:0006436'); --GO:0047761

-- Step 3: Export the final results to a TSV file
COPY (SELECT new_subject, object FROM final_result) TO '/Users/brooksantangelo/Documents/Repositories/kg-microbe-projects/kg_microbe_paper_2024/Scripts/Intermediate_Files_Competencies/butyrate/NCBI_genomic_traits_GO.tsv' (FORMAT 'csv', DELIMITER '\t', HEADER TRUE);



# # Assessing how many bugs we have from uniprot in this graph
# create or replace table uniprot_edges as 
# select 
#     subject,
#     predicate,
#     object,
# from read_csv('uniprot_genome_features/edges.tsv', filename=true, union_by_name=true);

# SELECT *
# FROM uniprot_edges
# WHERE object LIKE 'NCBITaxon:%'
#   AND subject LIKE 'Proteomes:%'
#   AND predicate = 'biolink:derives_from';


# # See how many edges are in the biomedical-traits graph
# create or replace table base_edges as 
# select 
#     subject,
#     predicate,
#     object,
# from read_csv('merged-kg_kg-microbe-humanasahost-traits_edges.tsv', filename=true, union_by_name=true);

# select count(distinct subject || ',' || predicate || ',' || object) as unique_edges_count
# from base_edges;