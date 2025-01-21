# Additional-Scripts

## 1. `conservedProteins.py`

This script processes orthogroup data for a set of protein sequences and selects conserved orthogroups present in at least 50% of species.

1. Mapping Protein IDs to Species:
   - Reads protein FASTA files for multiple species.
   - Maps protein IDs to species using a predefined dictionary of species names.

2. Processing Orthogroup Data:
   - Reads the `Orthogroups.txt` file containing protein orthogroups.
   - Assigns species information to each orthogroup based on the protein-to-species mapping.

3. Filtering Orthogroups:
   - Calculates the total number of species and determines a minimum threshold (50% of total species).
   - Selects orthogroups that are present in at least the threshold number of species.

4. Saving Results:
   - Writes the selected orthogroups and their protein information to a text file (`Conserved_Orthogroups.txt`).

5. Preparing Protein Sequences:
   - Extracts protein sequences for the selected orthogroups.
   - Saves them to a combined FASTA file (`Conserved_Proteins.faa`) with headers modified to include species names, making the file ready for tools like BRAKER2/3.
  
**Input Files**

- Protein FASTA Files: contains protein sequences for different species. (`e.g., ncbi_dataset/data/GCA_000372725.1/GCA_000372725.1.faa`)  
- Orthogroups.txt File:  contains orthogroup data (with each line formatted as:  `Orthogroup_ID: Protein_ID1 Protein_ID2 ... Protein_IDN`; e.g., `Orthogroups.txt` )  

**Output Files**

- Conserved Orthogroups File: text file listing orthogroups that are present in at least 50% of species, along with their associated proteins (e.g., `Conserved_Orthogroups.txt`)
- Combined Protein FASTA File: single FASTA file containing protein sequences for the selected orthogroups, headers are modified to include species names (e.g., `Conserved_Proteins.faa`)

**Applications**

- Identification of conserved orthogroups for evolutionary or comparative genomic studies.  
- Preparing input files for tools like BRAKER2/3 for gene prediction or annotation refinement.  
- Functional and taxonomic analysis of conserved protein families across species.  


## 2. `categorizeGOterms.py`

This script processes Gene Ontology (GO) annotations for transcripts to categorize GO terms into their primary categories and analyze their distribution across a proteome.

1. GO Term Parsing:  
   - Downloads and parses the `go-basic.obo` file to extract GO term categories:  
     - Biological Process (BP): Functional processes like mitosis or photosynthesis.  
     - Molecular Function (MF): Activities like enzymatic functions or binding.  
     - Cellular Component (CC): Locations like the nucleus or plasma membrane.  

2. Data Categorization:  
   - Processes a GO annotations file (e.g., `Asterionella_formosa_topgo.txt`) containing transcript IDs and associated GO terms.  
   - Categorizes GO terms for each transcript into BP, MF, and CC based on the parsed data.  

3. Distribution Analysis:  
   - Identifies the distribution of GO annotations for transcripts across single, double, or all three categories (BP, MF, CC).  
   - Counts transcripts in different category combinations.  

4. Visualization:  
   - Creates a bar plot to visualize the distribution of transcripts across GO categories.  

**Input Files**

- GO annotations file: Contains transcript IDs and comma-separated GO terms (e.g., `Asterionella_formosa_topgo.txt`).  
- GO term-to-category mapping file: Generated from `go-basic.obo` (e.g., `go_categories.csv`).  

**Outputs**

- Category Distribution Summary:  the number of transcripts annotated in one, two, or all three GO categories.
- Enhanced GO Data: detailed dataset with transcript IDs, GO terms, and their categorization into BP, MF, and CC.  
- Visualization: bar plot showing the distribution of transcripts across GO category combinations.  

**Applications**

- Comparative proteomics to assess GO annotation distributions in single-exon genes or across species.  
- Functional annotation analysis for transcriptomic studies. 


## 3. `gettopGO.py`

**Input Files**

- Gene-to-GO Mapping File:  contains mappings of gene IDs to Gene Ontology (GO) terms (`input/Asterionella_formosa_topgo.txt`)
  Format:  
  ```
  Gene_ID<TAB>GO:ID1,GO:ID2,...
  ```  
- Gene Ontology DAG File (go-basic.obo):  provides the hierarchical structure of GO terms (Downloaded automatically from: `http://purl.obolibrary.org/obo/go/go-basic.obo`, saved to: `input/go-basic.obo`)  

**Output Files**
- Human-Readable Gene-to-GO Mapping File:  converts GO IDs into descriptive terms for each gene (`Asterionella_formosa_GOreadable.txt`).  
  Format:  
  ```
  Gene<TAB>Description
  ```    
- Pie Chart of Top GO Terms: a visual representation of the most frequent GO terms in the dataset (`Asterionella_formosa_GOterms_piechart.png`)

**Applications**

- Analyze the distribution of Gene Ontology (GO) terms associated with genes for functional annotation.  
- Visualize the most frequently occurring GO terms in a dataset to identify key biological processes, molecular functions, and cellular components.  
- Generate human-readable gene-to-GO mappings for downstream analysis or reporting.  


## 4. `fantasiaResultsProcessing.py`

**Input Files**

- Gene-to-GO Mapping Files:  each file maps species to GO terms (`input/*_topGo.txt`).
  Format:  
  ```
  Gene_ID<TAB>GO:ID1,GO:ID2,...
  ```  

**Output Files**

- Binary GO Term Matrix: csv with species as rows and GO terms as columns, indicating presence (`1`) or absence (`0`) of terms (`output/FantasiaAnalysis/go_term_binary_matrix.csv`).  
- PCA Plot: plot to visualize species based on their GO term profiles (`output/FantasiaAnalysis/pca_plot.png`)
- Focus Categories Analysis: csv summarizing the count of GO terms in each category for each species (`output/FantasiaAnalysis/focus_categories_analysis.csv`).     
- Hierarchical Clustering Dendrogram: dendrogram of species based on focus categories (`output/FantasiaAnalysis/hierarchical_clustering_dendrogram.png`) 
- Elbow Plot: this plot showing the optimal number of clusters for K-means clustering (`output/FantasiaAnalysis/elbow_plot.png`)
- K-means Clusters:  cSV with species and their K-means cluster assignments (`output/FantasiaAnalysis/kmeans_clusters.csv` ).  
- Cluster Heatmap: heatmap of focus categories across K-means clusters (`output/FantasiaAnalysis/cluster_heatmap.png`).  


**Applications**

This script performs the following steps:

1. GO Term Matrix Creation: generates a binary matrix from gene-to-GO mapping files.

2. PCA: reduces dimensionality and visualizes species based on their GO term profiles.

3. Focus Categories Analysis: counts the number of GO terms in predefined categories for each species.

4. Hierarchical Clustering: clusters species based on focus categories and generates a dendrogram.

5. Elbow Method: determines the optimal number of clusters for K-means by analyzing inertia.

6. K-means Clustering: performs K-means clustering on species and saves the results.

7. Cluster Heatmap: visualizes the distribution of focus categories across K-means clusters.

This analysis helps identify functional similarities among species based on their GO term profiles.
