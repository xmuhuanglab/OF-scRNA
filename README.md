# OF-scRNA

<div align=center>
<img src="https://user-images.githubusercontent.com/95668602/173569851-05147ce0-347a-4c2d-98f8-c3e753ac05c9.png" width="400px"/>
</div>

**A framework to identify and characterize oncofetal cells of tumor microenvironment based on single-cell RNA sequencing: a case study for hepatocellular carcinoma**

## Overview
OF-scRNA is a framework designed to answer the following questions:
- Are there onco-fetal cells in the tumor?
- What are the similarities between onco-fetal cells and fetal cells? Why do onco-fetal cells reacquire fetal similarities? Which cells are transformed into onco-fetal cells?
- Is there clinical translation value of fetal-like features?

## Installation
We wrapped the OF-scRNA framework as a command line program: **OF-scRNA.py**<br />
You can download it via:
<div align=center>
<img src="https://user-images.githubusercontent.com/95668602/173772517-da6ce32d-f0db-4f65-877e-f62578e9d2df.png" width="400px"/>
</div>
and unzip the file

```bash
unzip OF-scRNA-main.zip

# enter the folder
cd OF-scRNA-main
```

**Data** in https://drive.google.com/file/d/1k26hI6hjvg7WQpZf4Ps7IcTMAp1rnzYl/view?usp=sharing.
You have to download it in OF-scRNA-main folder and unzip it via:
```bash
# You may need to modify file permissions
# chmod 750 data.tgz
tar -zvxf data.tgz
```

We also wrapped the docker image for the program to run:
```bash
docker pull huaqianghuang/of-scrna:0.1
```

## Usage

### 1. Build the container to run OF-scRNA

> You must modify the contents of [XXXXXXX]

```bash
# -v host_folder:container_folder  | Mount the host folder
docker run --name=of -dit -p 8001:8001 -h of --restart unless-stopped -v [/the_absolute_path_you_unzip_the_.zip]:/cluster/huanglab/hhuang/project/Cancer_dev/Liver/OF-scRNA/ huaqianghuang/of-scrna:0.1

# Now you enter the container
docker exec -it of bash

# change to the work directory
cd /cluster/huanglab/hhuang/project/Cancer_dev/Liver/OF-scRNA/
```

### 2. Workflow of OF-scRNA

#### Usage of OF-scRNA.py

```bash
python OF-scRNA.py -h

python OF-scRNA.py [modulename] -h

# such as
python OF-scRNA.py findoncofetal -p ./parameters/parameters.txt
```

#### 2.1 Identification of onco-fetal cells

Cells that **specifically enrich in tumors** and **resemble fetal cells** were defined as candidate **onco-fetal cells**

So, we use two score to identify onco-fetal cells

**relative enrichment score**: specifically enrich in tumor

**connectivity score**: resemble fetal cells

```bash
# I have run the following steps so you don't need to run
# 1. First extract celltype yout want interest, such as fibroblast
# python OF-scRNA.py extract-cells -p ./parameters/parameters.txt
# 2. Then find the highly variable genes, do PCA
# python OF-scRNA.py scanpy-workflow -p ./parameters/parameters.txt
# 3. Use cluster similarity spectrum algorithm to correct batch effect, then clustering
# python OF-scRNA.py css -p ./parameters/parameters.txt
# -------------- After subtype annotation, we get the data [Fibro_css_type.h5ad] ------
```

**Input:**

Anndata object, have log normalize expression matrix, high variable genes and PCA, CSS (X_css dimension), batch label, subtype label

```bash
# 4. To examine whether there are onco-fetal cells
python OF-scRNA.py findoncofetal -p ./parameters/parameters.txt
```

**output:** /1.Identification/

**FindOncofetal.pdf**: Heatmap display relative enrichment score and connectivity score
<div align=center>
<img src="https://s2.loli.net/2022/05/12/2Pjs9JlFTxz3CQI.png"/>
</div>

**relative_enrichment.csv**: The data to plot heatmap of relative enrichment score

**paga_subset.csv**: The data to plot heatmap of connectivity score

**Fibro_final.rds**: The Seurat object after calculate connectivity score

**Fibro_final.h5ad**: The anndata object after calculate connectivity score

#### 2.2 Characterization of onco-fetal cells through transcriptional program

To evaluate the similarity of transcriptional program between mCAF and Fetal mFibro and find shared genes

**Input:**

Anndata object, have log normalize expression matrix, high variable genes, PCA, CSS (X_css dimension) and subtype label

```bash
# 5. find shared genes
python OF-scRNA.py findgene -p ./parameters/parameters.txt
```

**output:** /2.Gene/

**top50-markers-jaccard.pdf**: Jaccard similarity (top) of top50 markers from fetal enriched subtypes (y axis) with top50 markers of tumor or adjacent normal enriched subtypes (x axis)

<div align=center>
<img src="https://s2.loli.net/2022/05/12/IdGU5zZpjqcfAF7.png"/>
</div>

**tn-top50-markers.csv**: top50 differentially express genes between adult tissue subtypes, only gene name

**fetal-top50-markers.csv**: top50 differentially express genes between fetal liver tissue subtypes, only gene name

**full-tn-top50-markers.csv**: top50 differentially express genes between adult tissue subtypes, contains gene name, pvalue, logFC, score (criteria to choose top50)

**full-fetal-top50-markers.csv**: top50 differentially express genes between fetal tissue subtypes, contains gene name, pvalue, logFC, score (criteria to choose top50)

**overlap_genes.csv**: Overlapping genes of differentially express genes from onco-fetal cells and its corresponding fetal cells

**oncofetal.csv**: The onco-fetal signature

**non_oncofetal.csv**: The non onco-fetal signature

#### 2.3 Characterization of onco-fetal cells through gene regulatory networks

To better understand the onset of embryonic reprogramming of tumor cells and find the shared TFs

**Input:**

1. Seurat object on celltype level, have log normalize expression matrix, subtype label
2. Cistarget database (https://resources.aertslab.org/cistarget/)

```bash
# 6. find shared TFs
python OF-scRNA.py findgrn -p ./parameters/parameters.txt
```

**output:** /3.GRN

**correlation-heatmap.pdf**: Correlation of binary regulon activity percentage between each subtype. 
<div align="center">
<img src="https://s2.loli.net/2022/05/12/jwgoNthKI13Hdn2.png">
</div>

**RSS-filtering.pdf**: Dot plot shows regulons that onco-fetal cells or its corresponding fetal cells turn on after filtering regulator specific score
<div align="center">
<img src="https://s2.loli.net/2022/05/12/gHSkCODsQYT28oG.png">
</div>

**rss-table.csv**: Dataframe of the regulon specific score of each regulon in each subtype

#### 2.4 Characterization of onco-fetal cells through Cellular communication

To investigate how fetal signals reacquired by onco-fetal cells affect TME

**Input:**

tumor_data, normal_data, fetal_data: seurat object, have log normalize expression matrix, type label (subtypes + other celltypes)

```bash
# 7. find shared LR-pairs
python OF-scRNA.py findccc -p ./parameters/parameters.txt
```

**output:** /4.CCC/

**tumor.cellchat.orig.rds, normal.cellchat.orig.rds, fetal.cellchat.orig.rds**: The CellChat object created at the beginning

**tumor.cellchat.final.rds, normal.cellchat.final.rds, fetal.cellchat.final.rds**: The CellChat object after run cellchat workflow

**outgoing-compare.pdf**: Number of significant ligand-receptor pairs from onco-fetal_cells/normal_cells/corresponding_fetal_cells to other cells in different sample origin. The edge width is proportional to the indicated number of ligand-receptor pairs.

<div align="center">
<img src="https://s2.loli.net/2022/05/12/SaXUVg24xl3RrQw.png">
</div>

**LR-bubble.pdf**: Dot plot of shared ligand-receptors of onco-fetal cells and its corresponding fetal cells. 
<div align="center">
<img src="https://user-images.githubusercontent.com/95668602/173576366-44d5f4ba-27c0-487f-b6c2-2f3f22dde12d.png" width="500px"><img src="https://user-images.githubusercontent.com/95668602/173576558-b5aabb5c-6d1b-41d6-8cb9-1edd13826182.png" width="500px">
</div>

#### 2.5 Trajectory inference

To investigate the potetion origin of onco-fetal cells

**Input:**

spliced.csv: matrix of spliced mRNA of tumor and normal cells

unspliced.csv: matrix of spliced mRNA of tumor and normal cells

cell_embeddings.csv: Before UMAP

metadata.csv: cellinfo, such as subtype label

```bash
# 9. find the potential origin of onco-fetal cells
python OF-scRNA.py trajectory -p ./parameters/parameters.txt
```

**output:** /5.Trajectory/

**scvelo_velocity.png**: Steady state RNA velocity of tumor and adjacent normal fibroblasts
<div align="center">
<img src="https://s2.loli.net/2022/05/12/mEUkFpicjJQzZfC.png">
</div>

**scvelo_diff_genes.csv**: The differentially regulated genes of each subtpye

#### 2.6 Clinical relevance of onco-fetal cells

To inverstigate the clinical relevance of onco-fetal cells

**Input:**

1. Bulk rna data (gene in row, patient in column)
2. Clinical information (have overall survival (time), state )
3. Onco-fetal signature

```bash
# 10. evaluate the clinical relevance of onco-fetal signature
python OF-scRNA.py clinical -p ./parameters/parameters.txt
```

**output:** /6.Clinical/

**oncofetal-gsva-TCGA-KM-plot.pdf**: Kaplan-Meier curves for overall survival of HCC patients in the TCGA ( cohort stratified by GSVA score of mCAF onco-fetal signature. The p value was calculated by the log-rank test.
<div align="center">
<img src="https://s2.loli.net/2022/05/12/a9s1lzxqguVkMUO.png">
</div>
