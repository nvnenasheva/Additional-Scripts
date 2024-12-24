import pandas as pd
import requests
import os
import matplotlib.pyplot as plt
import numpy as np


'''
This script processes GO (Gene Ontology) annotations for transcripts and categorizes the GO terms into their primary categories:
Biological Process (BP): Processes that contribute to the functioning of cells, tissues, or organisms (e.g., mitosis, photosynthesis).
Molecular Function (MF): The activities performed by a gene product at the molecular level (e.g., enzymatic activity, binding
Cellular Component (CC): The location or structure where a gene product is active (e.g., nucleus, plasma membrane).

The script calculates the distribution of these categories across the proteome to identify transcripts with annotations in one, two, 
or all three categories. 
This is particularly useful for analyzing the distribution of GO annotations in single-exon genes or across different species.

inputs:
- GO annotations file, FANTASIA output (Asterionella_formosa_topgo.txt, for example)
columns: Transcript ID, Comma-separated GO terms for each transcrip
- GO term-to-category mapping file (go_categories.csv)
columns: GO_term, category (BP,MF,CC)

output:
- category destribution summary: number of transcripts with annotations in one, two, or all three categories
- enchanced GO data: transcript ID, GO terms, and teir classification into BP, MF, CC categories
'''

# parse go-basic.obo file and extract GO-term categorires
def parse_obo(file_path):
    go_mapping = {}
    with open(file_path, 'r') as f:
        current_id = None
        current_namespace = None
        for line in f:
            line = line.strip()
            if line.startswith("id: GO:"):
                current_id = line.split(": ")[1]
            elif line.startswith("namespace: "):
                current_namespace = line.split(": ")[1]
                if current_namespace == "biological_process":
                    go_mapping[current_id] = "BP"
                elif current_namespace == "molecular_function":
                    go_mapping[current_id] = "MF"
                elif current_namespace == "cellular_component":
                    go_mapping[current_id] = "CC"
    return go_mapping

# download go-basic.obo file
def download_obo(url, output_file):
    try:
        print("Downloading go-basic.obo file...")
        # Send a GET request to download the file
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Check for HTTP request errors

        # Write the content to a local file
        with open(output_file, "wb") as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

        print(f"Download completed. File saved as: {output_file}")

    except requests.exceptions.RequestException as e:
        print(f"Error during download: {e}")

    return output_file

# categorize GO terms into BP, MF, and CC
def categorize_go_terms(go_terms, go_categories_df):
    go_terms_list = go_terms.split(", ")
    categories = go_categories_df.set_index("GO_Term")["Category"]
    term_categories = [categories.get(go_term, None) for go_term in go_terms_list]
    return {"BP": term_categories.count("BP"), 
            "MF": term_categories.count("MF"), 
            "CC": term_categories.count("CC")}

# visualize the distribution of GO categories
def visualize_distribution_stacked(distribution, title_name):
    distribution_df = distribution.reset_index()
    distribution_df.columns = ["Has_BP", "Has_MF", "Has_CC", "Count"]

    # Categorize the combinations
    def categorize(row):
        categories = []
        if row["Has_BP"]:
            categories.append("BP")
        if row["Has_MF"]:
            categories.append("MF")
        if row["Has_CC"]:
            categories.append("CC")
        
        if len(categories) == 3:
            return "All 3 Categories (BP, MF, CC)"
        elif len(categories) == 2:
            return f"{categories[0]} and {categories[1]}"
        elif len(categories) == 1:
            return f"Only {categories[0]}"
        else:
            return "No Category"

    distribution_df["Category_Combination"] = distribution_df.apply(categorize, axis=1)

    # Aggregate counts by the new categories
    grouped_data = distribution_df.groupby("Category_Combination")["Count"].sum().reset_index()

    # Sort categories logically
    category_order = [
        "All 3 Categories (BP, MF, CC)",
        "BP and MF", "BP and CC", "MF and CC",
        "Only BP", "Only MF", "Only CC",
        "No Category"
    ]
    grouped_data["Category_Combination"] = pd.Categorical(grouped_data["Category_Combination"], categories=category_order, ordered=True)
    grouped_data = grouped_data.sort_values("Category_Combination")

    # Plot the distribution
    fig, ax = plt.subplots(figsize=(3, 10))
    bars = ax.bar(grouped_data["Category_Combination"], grouped_data["Count"], color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'])

    ax.set_xlabel("GO Category Combinations", fontsize=12, fontweight="bold")
    ax.set_ylabel("Number of Transcripts", fontsize=12, fontweight="bold")
    ax.set_title(title_name, fontsize=14, fontweight="bold")
    ax.bar_label(bars, padding=3)
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.tight_layout()
    plt.show()

# ---------------------------------------------------------------------------------------------------------------------
# Step 1: prepare files and data

# Load GO annotation file
input_dir = 'input'
file_name='Asterionella_formosa_topgo.txt'
go_annotation_file = os.path.join(input_dir, file_name)
go_data = pd.read_csv(go_annotation_file, sep="\t", header=None, names=["Transcript", "GO_Terms"])

# get go-basic.obo file
file_name='Skeletonema_potamos_topgo.txt'
url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
go_basic = "/home/natalia/PycharmProjects/addtional_scripts/input/go-basic.obo"
go_categories_file = "/home/natalia/PycharmProjects/addtional_scripts/output/go_categories.csv"

if not os.path.exists(go_basic):
    download_obo(url, go_basic)

# Parse the go-basic.obo file
print("Parsing go-basic.obo file...")
if not os.path.exists(go_categories_file):
    go_mapping = parse_obo(go_basic)
    go_categories_file = "go_categories.csv"
    go_mapping_df = pd.DataFrame(list(go_mapping.items()), columns=["GO_Term", "Category"])
    go_mapping_df.to_csv("/home/natalia/PycharmProjects/addtional_scripts/output/go_categories.csv", index=False)

    print(f"GO term-to-category mapping saved to {go_categories_file}")

# ---------------------------------------------------------------------------------------------------------------------
# Step 2: categorize GO terms

# Load the GO annotation file
print("Loading GO annotation file...")
go_data = pd.read_csv(go_annotation_file, sep="\t", header=None, names=["Transcript", "GO_Terms"])

# Load the GO term-to-category mapping file
print("Loading GO categories...")
go_categories = pd.read_csv(go_categories_file)


# Categorize GO terms for each transcript
print("Categorizing GO terms...")
go_data["GO_Categories"] = go_data["GO_Terms"].apply(
    lambda x: categorize_go_terms(x, go_categories)
)

# Analyze the distribution
go_data["Has_BP"] = go_data["GO_Categories"].apply(lambda x: x["BP"] > 0)
go_data["Has_MF"] = go_data["GO_Categories"].apply(lambda x: x["MF"] > 0)
go_data["Has_CC"] = go_data["GO_Categories"].apply(lambda x: x["CC"] > 0)

# Count transcripts in different combinations
distribution = go_data.groupby(["Has_BP", "Has_MF", "Has_CC"]).size()

'''
example: 

data = {
    (False, True, True): 21,
    (True, True, True): 16745,
}
distribution = pd.Series(data)
'''


print(distribution)

# Visualize the results
title = file_name.split('_topgo.txt')[0].replace('_', ' ')
visualize_distribution_stacked(distribution, title_name=title)
