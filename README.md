# DOCUMENTATION of the Codes for Interactome
This document illustrates the use of the codes to implement the methods described in the Interactome article
## Table of contents
1. [Pre requisites](#pre-requisities)
2. [Similarity Network Fusion(SNF)](#similarity-network-fusion)
    1. [Single Biome Clustering](#single-biome-clustering)
    2. [Dual Biome Clustering](#dual-biome-clustering)
    3. [Merged Biome Clustering](#merged-biome-clustering)
3. [Weighted SNF](#weighted-snf)
4. [Co-occurence analysis](#co-occurence-analysis)
## Pre requisites
You will need the following softwares and packages to run the codes flawlesssly.
For Linux systems, Run the following 
```
sudo apt-get update && apt-get install -y python2.7 \
	python-pip \
	libcurl4-openssl-dev \
	libopenblas-base \
	gdebi-core \
	wget \
	default-jdk\
	gfortran\
	r-base
```
1. Python 2.7
    - [pandas](https://pandas.pydata.org/)
    - [numpy](https://numpy.org/)
    - [sklearn](https://scikit-learn.org/stable/index.html)
2. R 3.5.1
    - [SNFtool](https://cran.r-project.org/web/packages/SNFtool/index.html)
    - [vegan](https://cran.r-project.org/web/packages/vegan/index.html)
    - [dunn.test](https://cran.r-project.org/web/packages/dunn.test/index.html)
    - [reticulate](https://cran.r-project.org/web/packages/reticulate/index.html)
    - [mboost](https://cran.r-project.org/web/packages/mboost/index.html)
    - [boot](https://cran.r-project.org/web/packages/boot/boot.pdf)
    - [gMCP](https://cran.r-project.org/web/packages/gMCP/index.html)
    - [minet](https://www.bioconductor.org/packages/release/bioc/html/minet.html)

**Note:**
For windows and other operating system try running [this docker image](https://hub.docker.com/repository/docker/jayanthkumar/co-occurance_analysis), to implement co-occurence analysis.
## Similarity Network Fusion
Navigate to `./SNF_Analysis` 
``` bash
SNF_Analysis
   ├── Data
   ├── Dual_biome_clustering
   │   ├── results
   │   ├── tables
   │   └── Tuning_k
   │       ├── b+f
   │       ├── b+v
   │       └── f+v
   ├── Merged_biome_clustering
   │   ├── results
   │   └── Tuning_k
   └── Single_biome_clustering
       ├── results
       └── Tables
```
The datasets used in the article are available in the `./Data` directory.
### Single biome clustering
___
To run single biome clustering, execute the below command
``` bash 
Rscript ./Single_biome_clustering/snf.R
```
This code uses spectral clustering to cluster the individual biomes and charcterize the clusters using appropriate statistical tests. The results of the analysis, ie the cluster labels and similarity matrices  are written into the `results` directory as `./results/*_labels.csv` and `./results/*_matrix.csv` where `*` represents bacteria, fungi or virus datasets.   
The cluster charaterization results are written as `./Single_biome_clustering/*_Dunn-test.txt` and `./Single_biome_clustering/* table.csv` where `*` is bacteria, fungi or virus.

**Calculating the silhouette scores of the clusters**
To calculate the silhouette scores, we use the function `silhouette_score` defined in `./results/sil.py`. 
The following python code can be used to calculate the silhouette values of the clusters
```python
import sil
arr=pandas.read_csv("filename")
labels=pandas.read_csv("labels")
silhouette_score(arr,labels)
```
### Dual biome clustering
___
To regenerate clustering results of dual biome merging, i.e "bacteria+fungi", "fungi+virus" and "virus+bacteria", execute the below command
```bash
Rscript ./Dual_biome_clustering/snf.R
```
This code uses Similarity Network Fusion(SNF) to merge 2 microbiomes and clusters the integrated biome using spectral clustering. Further, appropriate statistical tests are used to charcterize the clusters. The results of the analysis, i.e. the cluster labels and similarity matrices  are written into the `results` directory as `./results/*_labels.csv` and `./results/*_matrix.csv` where `*` represents "bacteria+fungi", "fungi+virus"..etc datasets.   
The cluster charaterization results are written out as `./Dual_biome_clustering/*_Dunn-test.txt` and `./Dual_biome_clustering/* table.csv` where `*` "bacteria+fungi", "fungi+virus"..etc.

**Silhouette scores of the clusters can be computed as described above**
### Merged biome clustering
___
To regenerate clustering results of integrated microbiome , i.e "bacteria+fungi+virus", execute the below command
```bash
Rscript ./Merged_biome_clustering/snf.R
```
This code uses Similarity Network Fusion(SNF) to integrate all three microbiomes and cluster them using spectral clustering. Further, appropriate statistical tests are used to charcterize the clusters. The results of the analysis, i.e. the cluster labels and similarity matrices are written into the `results` directory as `./results/labels.csv` and `./results/matrix.csv` 
The cluster charaterization results are written out as `./Dual_biome_clustering/Dunn-test.txt` and `./Dual_biome_clustering/table.csv`.
___
## Weighted SNF
The weighted SNF analysis codes are to executed similar to SNF analysis codes as described above
Below is the directory structure of `./wSNF_Analysis`
```
wSNF_Analysis
	├── All_biomes
	│   ├── function_snf.R
	│   ├── results
	│   │   ├── labels.csv
	│   │   └── matrix.csv
	│   ├── sil.py
	│   ├── snf.R
	│   └── table.csv
```
`SNF_weighted_iter` function of `function_snf.R` is the modified similarity network fusion that accounts for weights.
**Note:** weighted SNF needs atleast 3 similarity matrices 
___
## Co-occurence analysis
This section illustrates the implementation of co-occurence analysis as described in the article.
```
Co-occurance_*
	├── Bray-Curtis
	├── Cytoscape_Visuvalization
	├── GBLM
	├── Merge_p-val_scores
	├── merge.py
	├── MI
	├── Microbes.csv
	├── Pearsons
	├── run.sh
	└── Spearman
```
Execute `python merge.py` from the `./Co-occurance_filter` directory to run filtered version of co-occurence analysis and for unfiltered version of co-occurence analysis execute `python merge.py` from the `./Co-occurance_no_filter`

`python merge.py` outputs `Microbes.csv` the base data that the program uses

**Co-occurence networks from different measures**
Execute the following shell script
```bash
chmod +x runs.sh
./run.sh
```
Running this would populate the results of different measures into their respective directories.
To merge these co-occurence networks
```bash
cd ./Merge_p-val_scores/
#Set the variable "cluster" to "cluster1" in Merge.R
Rscript Merge.R 
#Set the variable "cluster" to "cluster2" in Merge.R
Rscript Merge.R 
python plot_fun.py
```
The above code would result in outputs `cluster1_Adj_cyto.csv` and `cluster2_Adj_cyto.csv` which are then imported into cytoscape using an adjacency matrix reader for further downstream processing.
