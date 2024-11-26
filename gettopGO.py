import pandas as pd
import re
import os
import goatools
from goatools.obo_parser import GODag
from urllib.request import urlretrieve
import logging
import matplotlib.pyplot as plt

# Download the Gene Ontology DAG file if it doesn't exist
def download_go_dag(input_dir):
    obo_file = os.path.join(input_dir, 'go-basic.obo')
    print(f"Checking for GO DAG file at: {obo_file}")
    if not os.path.exists(obo_file):
        print("Downloading GO DAG file...")
        urlretrieve('http://purl.obolibrary.org/obo/go/go-basic.obo', obo_file)
        print("GO DAG file downloaded.")
    else:
        print("GO DAG file already exists. Skipping download.")
    return obo_file


# Calculate the total number of genes annotated
def calculate_annotation_percentage(input_file):
    """
    Calculates the percentage of genes annotated via gene family and/or similarity search methods.

    Parameters:
    - input_file: Path to the gene-to-GO mapping file.

    Returns:
    - percentage: Percentage of genes annotated via gene family/similarity search methods.
    - total_genes: Total number of genes in the dataset.
    - annotated_genes: Number of genes with at least one GO annotation.
    """
    total_genes = 0
    annotated_genes = 0

    with open(input_file, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            parts = line.split('\t')
            if len(parts) != 2:
                print(f"Line {line_number} is incorrectly formatted: {line}")
                continue  # Skip lines that don't have exactly two columns
            gene_id, go_terms = parts
            total_genes += 1
            # Clean up GO terms list
            go_terms_list = [go_term.strip() for go_term in go_terms.split(',') if go_term.strip()]
            if go_terms_list:
                annotated_genes += 1

    if total_genes > 0:
        percentage = (annotated_genes / total_genes) * 100
    else:
        percentage = 0.0

    return percentage, total_genes, annotated_genes

# Create a human-readable version of the gene-to-GO mapping file
def create_GO_readable():
    # Read the gene-to-GO mapping file
    # Format: Gene<TAB>GO:ID1,GO:ID2,...
    gene2go = {}
    input_file = os.path.join(input_dir, 'Asterionella_formosa_topgo.txt')
    print("Reading gene-to-GO mapping file...")
    with open(input_file, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            try:
                gene, go_terms = line.strip().split('\t')
                go_ids = go_terms.strip().split(',')
                gene2go[gene] = go_ids
            except ValueError:
                logging.error(f"Line {line_number} is incorrectly formatted: {line.strip()}")

    # Convert GO IDs to human-readable terms
    gene2go_human_readable = {}
    for gene, go_ids in gene2go.items():
        descriptions = []
        for go_id in go_ids:
            if go_id in godag:
                descriptions.append(godag[go_id].name)
            else:
                descriptions.append('Unknown GO term')
                logging.info(f"GO ID {go_id} not found in GO DAG.")
        gene2go_human_readable[gene] = descriptions
    print("Conversion complete.")

    # Convert the dictionary to a long-format DataFrame
    data = []
    for gene, descriptions in gene2go_human_readable.items():
        for desc in descriptions:
            data.append({'Gene': gene, 'Description': desc})

    gene2go_df = pd.DataFrame(data)

    # Save the resulting DataFrame to a file
    # Generate output filename based on input filename using regex
    input_filename = os.path.basename(input_file)
    output_filename = re.sub(r'_topgo(\.txt)?$', '_GOreadable.txt', input_filename)

    output_file = os.path.join(output_dir, output_filename)
    gene2go_df.to_csv(output_file, sep='\t', index=False)
    print(f"File saved to '{output_file}'.")
    return output_file

# Updated plot_pie_chart function with include_unknown parameter
def plot_pie_chart(data, category_column, count_column, top_n=None, title='Pie Chart of Categories',
                   colors=None, explode=None, shadow=False, include_unknown=True):

    df = data.copy()

    # Exclude 'Unknown GO term' category if include_unknown is False
    if not include_unknown:
        df = df[df[category_column] != 'Unknown GO term']

    df = df.sort_values(by=count_column, ascending=False)

    # If top_n is specified, select only the top N categories
    if top_n is not None:
        df = df.head(top_n)
        # Calculate 'Others' category if top_n is less than total categories
        total_others = data[count_column].sum() - df[count_column].sum()
        if total_others > 0:
            df = df.append({category_column: 'Others', count_column: total_others}, ignore_index=True)

    # Recalculate total_count after possibly excluding 'Unknown GO term'
    total_count = df[count_column].sum()

    df['Percentage'] = (df[count_column] / total_count) * 100
    df['Legend_Label'] = df[category_column] + ' (' + df['Percentage'].round(1).astype(str) + '%)'

    labels = df['Legend_Label']
    sizes = df[count_column]

    # Plot pie chart
    fig, ax = plt.subplots(figsize=(8, 8))

    wedges, texts = ax.pie(
        sizes,
        labels=None,  # Remove labels from slices
        startangle=140,
        colors=colors,
        explode=explode,
        shadow=shadow,
        textprops={'fontsize': 10}
    )

    # Add legend with category labels and percentages
    ax.legend(
        wedges,
        labels,
        title=category_column,
        loc='center left',
        bbox_to_anchor=(1, 0, 0.5, 1),
        fontsize=10
    )

    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title(title, fontsize=14)
    plt.tight_layout()
    plt.show()

    return fig, ax


# Logging configuration
logging.basicConfig(filename='go_conversion.log', level=logging.INFO)

# Ensure the input and output directories exist
input_dir = 'input'
output_dir = 'output'
for directory in [input_dir, output_dir]:
    if not os.path.exists(directory):
        os.makedirs(directory)
print("goatools version:", goatools.__version__)

# Load the Gene Ontology DAG
print("Loading Gene Ontology DAG...")
obo_file = download_go_dag(input_dir)
godag = GODag(obo_file)

# Create a human-readable version of the gene-to-GO mapping file
create_GO_readable()
df = pd.read_csv('output/Asterionella_formosa_GOreadable.txt', sep='\t')

# Calculate the percentage of genes assigned to GO categories
input_file = os.path.join(input_dir, 'Asterionella_formosa_topgo.txt')
percentage, total_genes, genes_with_go = calculate_annotation_percentage(input_file)
print(f"Total number of genes: {total_genes}")
print(f"Number of genes with GO annotations: {genes_with_go}")
#print(f"Percentage of genes annotated via gene family/similarity search: {percentage:.2f}%")


# Count the frequency of each GO term
go_term_counts = df['Description'].value_counts().reset_index()
go_term_counts.columns = ['Description', 'Count']


# Plot the pie chart of the top 10 GO terms
print("Plotting pie chart...")
fig, ax = plot_pie_chart(
    data=go_term_counts,
    category_column='Description',
    count_column='Count',
    top_n=10, 
    title='Top 10 GO Terms in Asterionella formosa',
    shadow=True,
    include_unknown=False  # !!! Set to False to exclude 'Unknown GO term'
)

# Save the plot to a file
plot_filename = 'Asterionella_formosa_GOterms_piechart.png'
plot_file = os.path.join(output_dir, plot_filename)
fig.savefig(plot_file, bbox_inches='tight')
print(f"Pie chart saved to '{plot_file}'.")
