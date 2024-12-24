import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


'''Output:
- binary GO term matrix: a CSV file with species as rows and GO terms as columns, where 1 indicates the presence of a term.
- PCA plot: a scatter plot of species based on their GO term profiles.
- Focus categories analysis: a summary of the number of GO terms in each focus category for each species.
- Hierarchical clustering dendrogram: a dendrogram showing the clustering of species based on focus categories.
- Elbow plot: a plot to determine the optimal number of clusters for k-means clustering.
- K-means clusters: a CSV file with species and their assigned clusters.
- Cluster heatmap: a heatmap showing the distribution of focus categories across clusters.
'''

input_folder = "/home/natalia/PycharmProjects/addtional_scripts/input"
output_folder = "/home/natalia/PycharmProjects/addtional_scripts/output/FantasiaAnalysis"

# ==================================FUNCTIONS======================================



def create_go_term_matrix(input_folder, output_file):
    """
    Parses GO terms from all files in the input folder, creates a binary matrix,
    and saves it to a CSV file.

    Args:
        input_folder (str): Path to the folder containing input files.
        output_file (str): Path to save the binary GO term matrix.
    """
    species_go_data = {}

    # Parse GO terms from each file
    for file_name in os.listdir(input_folder):
        if file_name.endswith("_topgo.txt"):
            file_path = os.path.join(input_folder, file_name)
            species_name = file_name.split('_topgo.txt')[0].replace('_', ' ')
            go_terms = parse_go_terms_from_file(file_path)  # Extract GO terms
            species_go_data[species_name] = go_terms

    # Create binary matrix
    all_go_terms = set(term for terms in species_go_data.values() for term in terms)
    matrix = pd.DataFrame(0, index=species_go_data.keys(), columns=sorted(all_go_terms))

    for species, terms in species_go_data.items():
        matrix.loc[species, terms] = 1

    # Save matrix to file
    matrix.to_csv(output_file)
    print(f"Binary GO term matrix saved to {output_file}")

    return matrix


def plot_pca(matrix, output_file=os.path.join(output_folder,"pca_plot.png")):
    # Remove the "Cluster" column if it exists
    go_matrix = matrix.drop(columns=["Cluster"]) if "Cluster" in matrix.columns else matrix
    
    # Perform PCA
    pca = PCA()
    pca_result = pca.fit_transform(go_matrix)
    
    # Compute explained variance
    explained_variance = pca.explained_variance_ratio_
    cumulative_variance = explained_variance.cumsum()
    
    # Scree Plot to show explained variance
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(explained_variance) + 1), cumulative_variance, marker='o', linestyle='--')
    plt.title('Explained Variance by PCA Components')
    plt.xlabel('Number of Components')
    plt.ylabel('Cumulative Explained Variance')
    plt.grid()
    plt.savefig("/home/natalia/PycharmProjects/addtional_scripts/output/components_plot.png")
    plt.show()
    
    # Decide the number of components (e.g., 2 for visualization)
    n_components = 2
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(go_matrix)
    
    # Convert PCA result to DataFrame
    pca_df = pd.DataFrame(pca_result, columns=[f"PC{i+1}" for i in range(n_components)], index=go_matrix.index)
    
    # Plot PCA
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x="PC1", y="PC2", data=pca_df, s=100)
    for i, species in enumerate(pca_df.index):
        plt.annotate(species, (pca_df.iloc[i, 0], pca_df.iloc[i, 1]), fontsize=9)
    plt.title("PCA of GO Terms")
    plt.xlabel("Principal Component 1")
    plt.ylabel("Principal Component 2")
    plt.grid()
    plt.savefig(output_file)
    plt.show()
    return pca_df


def parse_go_terms_from_file(file_path):
    """
    Parses GO terms from a file with the given format.

    Args:
    file_path (str): Path to the file.

    Returns:
    list: A list of unique GO terms found in the file.
    """
    go_terms = set()
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')  # Split by tab
            if len(columns) > 1:  # Ensure there's a second column
                terms = columns[1].split(', ')  # Split GO terms by comma
                go_terms.update(terms)  # Add terms to the set
    return list(go_terms)


def analyze_focus_categories(matrix, focus_categories):
    """
    Analyze GO term focus categories based on the provided focus_categories dictionary.

    Args:
    matrix (pd.DataFrame): Binary GO term matrix (species x GO terms).
    focus_categories (dict): Dictionary with focus categories and their GO terms.

    Returns:
    pd.DataFrame: DataFrame summarizing the counts for each category.
    """
    results = {}
    available_go_terms = set(matrix.columns)  # Extract available GO terms
    
    for category, go_terms in focus_categories.items():
        # Check for matching GO terms in the dataset
        matching_terms = [term for term in go_terms if term in available_go_terms]
        if matching_terms:
            counts = matrix[matching_terms].sum(axis=1)  # Sum matches per species
            results[category] = counts
        else:
            print(f"Category '{category}' skipped: no matching GO terms found in data.")
    
    return pd.DataFrame(results)

# to get the number of clusters to use as n_clusters in kmeans_clustering
def plot_hierarchical_clustering(focus_df, output_file=os.path.join(output_folder,"hierarchical_clustering_dendrogram.png")):
    # Perform hierarchical clustering using Ward's method
    linkage_matrix = linkage(focus_df, method="ward")

    # Plot dendrogram
    plt.figure(figsize=(12, 8))
    dendrogram(linkage_matrix, labels=focus_df.index, leaf_rotation=90, leaf_font_size=10)
    plt.title("Hierarchical Clustering of Species Based on Focus Categories")
    plt.xlabel("Species")
    plt.ylabel("Distance")
    
    # Save and show the plot
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

# to get the number of clusters to use as n_clusters in kmeans_clustering
def elbow_method(matrix, max_clusters=10, output_file=os.path.join(output_folder,"elbow_plot.png")):
    """
    Determines the optimal number of clusters using the Elbow Method.

    Args:
        matrix (pd.DataFrame): The dataset for clustering.
        max_clusters (int): The maximum number of clusters to test.
        output_file (str): File to save the elbow plot.

    Returns:
        None
    """
    inertias = []
    for k in range(1, max_clusters + 1):
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(matrix)
        inertias.append(kmeans.inertia_)

    # Plot the elbow curve
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, max_clusters + 1), inertias, marker='o', linestyle='--')
    plt.title('Elbow Method for Optimal Clusters', fontsize=16)
    plt.xlabel('Number of Clusters (k)', fontsize=12)
    plt.ylabel('Inertia (Sum of Squared Distances)', fontsize=12)
    plt.grid(True)
    plt.savefig(output_file)
    plt.show()

def kmeans_clustering(focus_df, n_clusters=3, output_file=os.path.join(output_folder,"kmeans_clusters.csv")):
    # Perform k-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    focus_df["Cluster"] = kmeans.fit_predict(focus_df)
    
    # Print cluster summary
    print(focus_df.groupby("Cluster").mean())  # Mean focus category values for each cluster
    
    # Save cluster assignments
    focus_df.to_csv(output_file)
    print(f"Cluster assignments saved to {output_file}")
    
    return focus_df

def plot_cluster_heatmap(focus_df, output_file=os.path.join(output_folder,"cluster_heatmap.png")):
    # Sort species by their cluster assignments
    sorted_df = focus_df.sort_values(by="Cluster")
    
    # Plot heatmap
    plt.figure(figsize=(14, 10))
    sns.heatmap(
        sorted_df.drop(columns=["Cluster"]),
        cmap="viridis",
        cbar=True,
        annot=True,
        linewidths=0.5,
        linecolor="gray",
        annot_kws={"fontsize": 8}  # Adjust annotation font size
    )
    plt.title("Heatmap of Focus Categories Grouped by Clusters", fontsize=16, weight="bold")
    plt.xlabel("Focus Categories", fontsize=12)
    plt.ylabel("Species", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=10)  # Rotate focus category names
    plt.yticks(fontsize=10)
    plt.tight_layout()  # Avoid clipping
    
    # Save and show the plot
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()



# ====================================RUN======================================

# Set the folder containing topGO files
matrix = create_go_term_matrix(input_folder, output_file=os.path.join(output_folder, "go_term_binary_matrix.csv"))


# Define focus categories
focus_categories = {
    "Photosynthesis": ["GO:0015979", "GO:0009765", "GO:0019684", "GO:0009654"],
    "Carbon Fixation": ["GO:0015977", "GO:0018313", "GO:0048488"],
    "Silicon Metabolism": ["GO:0043949", "GO:0009542", "GO:0071941"],
    "Nutrient Uptake": ["GO:0015698", "GO:0015706", "GO:0006820", "GO:0015695"],
    "Osmoregulation": ["GO:0006811", "GO:0009651", "GO:0009507"],
    "Lipid Metabolism": ["GO:0008610", "GO:0044241", "GO:0006633"],
    "Oxidative Stress Response": ["GO:0006979", "GO:0000302", "GO:0055114"],
    "Nitrogen Metabolism": ["GO:0042128", "GO:0019740", "GO:0006807"],
    "Iron Uptake": ["GO:0006826", "GO:0015682", "GO:0005506"],
    "Diatom-Specific Traits": ["GO:0042303", "GO:0009826", "GO:0016044"]
}

# Analyze focus categories
focus_df = analyze_focus_categories(matrix, focus_categories)
if focus_df.empty:
    print("Error: No matching GO terms found in focus categories. Please check your focus categories or input matrix.")
    exit()
focus_df.to_csv(os.path.join(output_folder, "focus_categories_analysis.csv"))


# Plot PCA of focus categories | did not help at all
# pca_df = plot_pca(matrix)


# Perform k-means clustering with elbow method | did not help at all
elbow_method(matrix, max_clusters=10) # The Elbow Method evaluates the sum of squared distances (inertia) between points and their assigned cluster centers. 
                                      # It identifies the point where adding more clusters results in diminishing returns in reducing the inertia.

# Perform hierarchical clustering
plot_hierarchical_clustering(focus_df) # with ward's linkage method

# Perform k-means clustering with 4 clusters
focus_df_kmeans = kmeans_clustering(focus_df, n_clusters=3)
plot_cluster_heatmap(focus_df_kmeans)
