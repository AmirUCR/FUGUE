import requests
import pandas as pd
import xml.etree.ElementTree as ET

df = pd.read_csv('fourdbs_hi_input_species.csv')


def get_taxonomic_hierarchy(species_name):
    # URL for esearch which searches and retrieves primary IDs (for example, taxonomy IDs) for a particular species name
    esearch_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={species_name}&retmode=json'
    
    # Send a request to the esearch URL
    response = requests.get(esearch_url)
    data = response.json()

    if len(data['esearchresult']['idlist']) == 0:
        return 'Not found'
    
    # Get the Taxonomy ID from esearch response
    taxonomy_id = data['esearchresult']['idlist'][0]

    # Retrieve taxonomy data using a taxonomy ID
    efetch_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxonomy_id}'
    
    response = requests.get(efetch_url)
    
    tree = ET.ElementTree(ET.fromstring(response.content))

    # Find the Lineage element and extract the taxonomic hierarchy
    lineage = None
    for elem in tree.iterfind('Taxon/Lineage'):
        lineage = elem.text

    return lineage

names = list()
sources = list()
found_class = list()

for idx, row in df.iterrows():
    if idx % 50 == 0:
        print(idx)
  
    new_name = row['species_name'].capitalize().replace('_', ' ')
    names.append(new_name)
    sources.append(row['source'])

    tax = get_taxonomic_hierarchy(new_name)
    found_class.append(tax)

species_w_classification = pd.DataFrame.from_dict({
    'species_name': names,
    'classification': found_class,
    'source': sources
})

species_w_classification.to_csv('species_w_classifications.csv', index=False)