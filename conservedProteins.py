#!/usr/bin/env python3

import os
import sys
from collections import defaultdict
import math
from Bio import SeqIO

# list of the input fasta files given to the OrthoFinder 
fasta_files = [
    'input/GCF_000002975.1.faa',
    'input/GCF_000018645.1.faa',
    'input/GCF_000194455.1.faa',
    'input/GCF_000315625.1.faa'
]

# create a mapping from protein IDs to species
protein_to_species = {}
all_species = set()

for fasta_file in fasta_files:
    species = os.path.basename(fasta_file).split('.')[0]
    # map species codes to descriptive names (optional), I used scientific name taken from NCBI (it will be used leater on in the script)
    species_names = {
        'GCF_000002975': 'Guillardia_theta',
        'GCF_000018645': 'Hemiselmis_andersenii',
        'GCF_000194455': 'Cryptomonas_paramecium',
        'GCF_000315625': 'Guillardia_theta_CCMP2712 '
    }
    species_name = species_names.get(species, species)
    all_species.add(species_name)
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            protein_id = record.id
            protein_to_species[protein_id] = species_name

# Step 1: create a mapping from protein IDs to species
protein_to_species = defaultdict(list)  # Use a list in case of duplicate protein IDs across species
all_species = set()

for fasta_file in fasta_files:
    # Use the base name of the file (without extension) as the species code
    species_code = os.path.basename(fasta_file).split('.')[0]
    species_name = species_names.get(species_code, species_code)
    all_species.add(species_name)
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            protein_id = record.id
            protein_to_species[protein_id].append(species_name)

# Step 2: read the Orthogroups.txt file and assign species using the mapping
orthogroups_file = 'input/Orthogroups.txt'
orthogroups = {}

with open(orthogroups_file, 'r') as og_file:
    for line in og_file:
        line = line.strip()
        if not line:
            continue
        og_id, proteins_str = line.split(':')
        proteins = proteins_str.strip().split()
        species_in_group = set()
        protein_entries = []
        for protein in proteins:
            species_list = protein_to_species.get(protein)
            if species_list:
                for species in species_list:
                    species_in_group.add(species)
                    protein_entries.append((species, protein))
            else:
                print(f"Warning: Protein ID {protein} not found in mapping.")
        orthogroups[og_id] = {
            'species': species_in_group,
            'proteins': protein_entries
        }

# Step 3: calculate total species and determine min_species (50% threshold)
total_species = len(all_species)
min_species = math.ceil(total_species * 0.5)
print(f"Total number of species: {total_species}")
print(f"Selecting orthogroups present in at least {min_species} species (50% of total species)")

# Step 4: select orthogroups present in at least min_species
selected_orthogroups = {}
for og_id, og_data in orthogroups.items():
    if len(og_data['species']) >= min_species:
        selected_orthogroups[og_id] = og_data

# Step 5: output the selected orthogroups and their proteins
output_file = 'output/Conserved_Orthogroups.txt'
with open(output_file, 'w') as out_file:
    for og_id, og_data in selected_orthogroups.items():
        proteins_formatted = [f"{species}|{protein}" for species, protein in og_data['proteins']]
        out_file.write(f"{og_id}: {' '.join(proteins_formatted)}\n")

print(f"Selected {len(selected_orthogroups)} orthogroups present in at least {min_species} species.")
print(f"Results written to {output_file}")

# Step 6: prepare protein sequences for BRAKER2 with descriptive file names
# combine all selected proteins into a single FASTA file suitable for BRAKER

extract_sequences = True  # Set to False if you don't need this

if extract_sequences:
    # get all protein IDs to extract, including species
    proteins_to_extract = set()
    for og_data in selected_orthogroups.values():
        for species, protein_id in og_data['proteins']:
            proteins_to_extract.add((species, protein_id))

    # create a combined FASTA file with a descriptive name
    output_fasta = 'output/Conserved_Proteins.faa'
    with open(output_fasta, 'w') as out_handle:
        for fasta_file in fasta_files:
            # get species code and name
            species_code = os.path.basename(fasta_file).split('.')[0]
            species_name = species_names.get(species_code, species_code)
            with open(fasta_file, 'r') as in_handle:
                for record in SeqIO.parse(in_handle, 'fasta'):
                    protein_id = record.id
                    if (species_name, protein_id) in proteins_to_extract:
                        # modify the header to include species name
                        record.id = f"{species_name}|{protein_id}"
                        record.description = ''
                        SeqIO.write(record, out_handle, 'fasta')

    print(f"Extracted protein sequences are saved in '{output_fasta}', ready for BRAKER2/3.")
