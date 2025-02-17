
import csv
import pandas as pd
import csv
import py4cytoscape as p4c
from py4cytoscape import gen_node_color_map, gen_edge_color_map
from py4cytoscape import palette_color_brewer_d_RdBu
from tqdm import tqdm

edges_file = "./Input_Files/kg-microbe-biomedical-function/merged-kg_edges.tsv"

output_lines = []
all_prefixes = []

# cmd = r"grep -o '\b[A-Za-z0-9_-]*:' merged-kg_edges.tsv | sed 's/:$//' | sort -u"
# result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

# # Convert output to a list of strings
# unique_prefixes = result.stdout.strip().split("\n")

# traits = ['oxygen', 'cell_length', 'motility', 'cell_width', 'NaCl_range', 'NaCl_opt', 'temp_delta', 'gc', 'pathways', 'pH_range', 'isolation_source', 'sporulation', 'trophic_type', 'pH_delta', 'strain', 'temp_range', 'pigment', 'gram_stain', 'solution', 'cell_shape', 'salinity', 'temperature', 'bacdive', 'pathogen', 'temp_opt', 'NaCl_delta', 'pH_opt']
# traits = ['assay', 'BSL', 'carbon_substrates', 'cell_length', 'cell_shape', 'cell_width',  'gc', 'gram_stain', 'http', 'isolation_source', 'medium', 'motility', 'NaCl_delta', 'NaCl_opt', 'NaCl_range', 'oxygen', 'pathogen', 'pathways', 'pH_delta', 'pH_opt', 'pH_range', 'pigment', 'production', 'salinity', 'sporulation', 'temp_delta', 'temperature', 'temp_opt', 'temp_range', 'trophic_type']

# with open(edges_file, "r") as file:
#     csv_reader = csv.DictReader(file, delimiter="\t")
#     for row in tqdm(csv_reader):

with open(edges_file, "r") as file:
    total_lines = sum(1 for _ in file) - 1  # Subtract 1 for the header

# Process with tqdm
with open(edges_file, "r") as file:
    csv_reader = csv.DictReader(file, delimiter="\t")
    for row in tqdm(csv_reader, total=total_lines, desc="Processing rows"):
        if not any(row.values()):  # Skip row if all values are empty
            continue
        subject_prefix = row["subject"].split(":")[0]
        # if subject_prefix in traits:
        #     subject_prefix = 'trait'
        predicate = row["predicate"].replace("biolink:","")
        object_prefix = row["object"].split(":")[0]
        # if object_prefix in traits:
        #     object_prefix = 'trait'
        # source = row["primary_knowledge_source"]
        output_lines.append('\t'.join([subject_prefix, predicate, object_prefix])) #, source]))
        all_prefixes.append(subject_prefix)
        all_prefixes.append(object_prefix)

    output_lines = sorted(set(output_lines))
# Add header
output_lines.insert(0,'\t'.join(["S","P","O"])) #,"Source"]))

with open("schema.tsv", 'w') as outfile:
    outfile.write('\n'.join(output_lines))

# Categories prefix types
all_prefixes = list(set(all_prefixes))
ontologies = ['MF', 'RO', 'NBO',  'OBO', 'OGMS', 'owl', 'UBERON', 'UPA', 'ENVO', 'PO', 'CHEBI', 'NCBITaxon', 'SO', 'NCIT', 'PCO', 'RHEA', 'MFOMD', 'IAO', 'GO', 'HP', 'BFO', 'FAO', 'EC', 'FOODON', 'MONDO', 'CARO', 'OBI', 'CHR', 'PATO', 'CL', 'CAS-RN', 'KEGG']

node_cateogories = []
node_cateogories.append('\t'.join(["Node","Category"]))

for prefix in all_prefixes:
    if prefix in ontologies:
        category = "Ontology"
    else: 
        category = "Dataset"
    node_cateogories.append('\t'.join([prefix, category]))

with open("schema.noa", 'w') as outfile:
    outfile.write('\n'.join(node_cateogories))

# Output to Cytoscape
subgraph_df = pd.read_csv("schema.tsv",sep='\t')
subgraph_df = pd.read_csv("schema_manual_update.tsv", sep='\t')
subgraph_attributes_df = pd.read_csv("schema.noa",sep='\t')

png_file = 'Subgraph_Visualization.png'

#Update column names for cytoscape
#Subset columns
subgraph_df = subgraph_df[['S','P','O']]
subgraph_df.columns = ['source','interaction','target']
subgraph_attributes_df.columns = ['id','index']

subgraph_attributes_df = subgraph_attributes_df[subgraph_attributes_df['id'].isin(subgraph_df['source']) | subgraph_attributes_df['id'].isin(subgraph_df['target'])]


p4c.create_network_from_data_frames(subgraph_attributes_df,subgraph_df,title='subgraph')

#Ensure no network exists named subgraph in Cytoscape or you will have to manually override before it can be output
p4c.set_visual_style('BioPAX_SIF',network='subgraph')

p4c.set_node_color_mapping(**gen_node_color_map('index', mapping_type='d',style_name='BioPAX_SIF'))

p4c.set_edge_color_mapping(**gen_edge_color_map('index', mapping_type='d',style_name='BioPAX_SIF'))

p4c.set_edge_label_mapping('interaction')

p4c.export_image(png_file,network='subgraph')
