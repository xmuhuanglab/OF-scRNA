# OF-scRNA

<div align=center><img src="https://user-images.githubusercontent.com/95668602/173528440-d558a8cb-ee50-41ec-a1aa-7b993e65b032.png"/></div>

**A framework to identify and characterize oncofetal cells of tumor microenvironment based on single-cell RNA sequencing: a case study for hepatocellular carcinoma**

## Overview
OF-scRNA is a framework designed to answer the following questions
- Are there onco-fetal cells in the tumor?
- What are the similarities between onco-fetal cells and fetal cells? Why do onco-fetal cells reacquire fetal similarities? Where the onco-fetal cells come from?
- 


## Usage

### 1. Copy the demo data and script (.tgz) 

> You must modify the contents of [XXXXXXX]

```bash
cp /cluster/huanglab/hhuang/project/Cancer_dev/Liver/OF-scRNA/tmp/demo.tgz [your_dir_to_save]

# eg. cp /cluster/huanglab/hhuang/project/Cancer_dev/Liver/OF-scRNA/tmp/demo.tgz ./tmp/
```

### 2. Unzip the file

```bash
tar -zxvf demo.tgz
```

**output:** /OF-scRNA/

**OF-scRNA.py**: the command line program to run

**of-scrna.tar**: the docker image

**code**: Rscript of the OF-scRNA framework

**Other directories**: some data

### 3. Build the container to run OF-scRNA

> You must modify the contents of [XXXXXXX]

```bash
cd OF-scRNA/

# You need to run docker in node2!!! Something is wrong with the node1 file system, can't write h5ad files
# You must have access to docker
# su docker-lab
# tugU6WScPq
docker load -i of-scrna.tar

# -v host_folder:container_folder  | Mount the host folder
docker run --name=of -dit -p 8001:8001 -h of --restart unless-stopped -v [/the_absolute_path_you_unzip_the_.tgz]:/cluster/huanglab/hhuang/project/Cancer_dev/Liver/OF-scRNA/ of-scrna:0.1

# Now you enter the container
docker exec -it of bash

# change to the work directory
cd /cluster/huanglab/hhuang/project/Cancer_dev/Liver/OF-scRNA/
```

### 4. Workflow of OF-scRNA

#### Usage of OF-scRNA.py

```bash
python OF-scRNA.py -h

python OF-scRNA.py [modulename] -h

# such as
python OF-scRNA.py findoncofetal -p ./parameters/parameters.txt
```



#### 4.1 Identification of oncofetal cells

Cells that **specifically enrich in tumors** and **resemble fetal cells** were defined as candidate **oncofetal cells**

So, we use two score to identify oncofetal cells

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
# 4. To examine whether there are oncofetal cells
python OF-scRNA.py findoncofetal -p ./parameters/parameters.txt
```

**output:** /1.Identification/

**FindOncofetal.pdf**: Heatmap display relative enrichment score and connectivity score

![image-20220512145929938](https://s2.loli.net/2022/05/12/2Pjs9JlFTxz3CQI.png)

**relative_enrichment.csv**: The data to plot heatmap of relative enrichment score

**paga_subset.csv**: The data to plot heatmap of connectivity score

**Fibro_final.rds**: The Seurat object after calculate connectivity score

**Fibro_final.h5ad**: The anndata object after calculate connectivity score



#### 4.2 Characterization of oncofetal cells through transcriptional program

To evaluate the similarity of transcriptional program between mCAF and Fetal mFibro and find shared genes

**Input:**

Anndata object, have log normalize expression matrix, high variable genes, PCA, CSS (X_css dimension) and subtype label

```bash
# 5. find shared genes
python OF-scRNA.py findgene -p ./parameters/parameters.txt
```

**output:** /2.Gene/

**top50-markers-jaccard.pdf**: Jaccard similarity (top) of top50 markers from fetal enriched subtypes (y axis) with top50 markers of tumor or adjacent normal enriched subtypes (x axis)

![image-20220512150441344](https://s2.loli.net/2022/05/12/IdGU5zZpjqcfAF7.png)

**tn-top50-markers.csv**: top50 differentially express genes between adult tissue subtypes, only gene name

**fetal-top50-markers.csv**: top50 differentially express genes between fetal liver tissue subtypes, only gene name

**full-tn-top50-markers.csv**: top50 differentially express genes between adult tissue subtypes, contains gene name, pvalue, logFC, score (criteria to choose top50)

**full-fetal-top50-markers.csv**: top50 differentially express genes between fetal tissue subtypes, contains gene name, pvalue, logFC, score (criteria to choose top50)

**overlap_genes.csv**: Overlapping genes of differentially express genes from oncofetal cells and its corresponding fetal cells

**oncofetal.csv**: The oncofetal signature

**non_oncofetal.csv**: The non oncofetal signature



#### 4.3 Characterization of oncofetal cells through gene regulatory networks

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

![image-20220512200937380](https://s2.loli.net/2022/05/12/jwgoNthKI13Hdn2.png)

**RSS-filtering.pdf**: Dot plot shows regulons that oncofetal cells or its corresponding fetal cells turn on after filtering regulator specific score

![image-20220512200957658](https://s2.loli.net/2022/05/12/gHSkCODsQYT28oG.png)

**rss-table.csv**: Dataframe of the regulon specific score of each regulon in each subtype



#### 4.4 Characterization of oncofetal cells through Cellular communication

To investigate how fetal signals reacquired by oncofetal cells affect TME

**Input:**

tumor_data, normal_data, fetal_data: seurat object, have log normalize expression matrix, type label (subtypes + other celltypes)

```bash
# 7. find shared LR-pairs
python OF-scRNA.py findccc -p ./parameters/parameters.txt
```

**output:** /4.CCC/

**tumor.cellchat.orig.rds, normal.cellchat.orig.rds, fetal.cellchat.orig.rds**: The CellChat object created at the beginning

**tumor.cellchat.final.rds, normal.cellchat.final.rds, fetal.cellchat.final.rds**: The CellChat object after run cellchat workflow

**outgoing-compare.pdf**: Number of significant ligand-receptor pairs from oncofetal_cells/normal_cells/corresponding_fetal_cells to other cells in different sample origin. The edge width is proportional to the indicated number of ligand-receptor pairs.

![image-20220512201058956](https://s2.loli.net/2022/05/12/SaXUVg24xl3RrQw.png)

**LR-bubble.pdf**: Dot plot of shared ligand-receptors of oncofetal cells and its corresponding fetal cells. 



#### 4.5 Regulation

To investigate whether hepatocytes can induce reprogramming of oncofetal cells

**Input:**

Seurat object, total cells, have log normalize expression matrix, celltype label, batch label

```bash
# 8. find which celltype drive the oncofetal reprogramming
python OF-scRNA.py drive -p ./parameters/parameters.txt
```

**output:** /4.CCC/

**Hepa_Fibro_ligand_activity_target.pdf**: Heatmap shows the regulatory potential of hepatocyte's gene regulate fibroblast gene

![image-20220512201159447](https://s2.loli.net/2022/05/12/X2Vw3ShGKbDnyB4.png)



#### 4.6 Trajectory inference

To investigate the potetion origin of oncofetal cells

**Input:**

spliced.csv: matrix of spliced mRNA of tumor and normal cells

unspliced.csv: matrix of spliced mRNA of tumor and normal cells

cell_embeddings.csv: Before UMAP

metadata.csv: cellinfo, such as subtype label

```bash
# 9. find the potential origin of oncofetal cells
python OF-scRNA.py trajectory -p ./parameters/parameters.txt
```

**output:** /5.Trajectory/

**scvelo_velocity.png**: Steady state RNA velocity of tumor and adjacent normal fibroblasts

![image-20220512201252187](https://s2.loli.net/2022/05/12/mEUkFpicjJQzZfC.png)

**scvelo_diff_genes.csv**: The differentially regulated genes of each subtpye



#### 4.7 Clinical relevance of oncofetal cells

To inverstigate the clinical relevance of oncofetal cells

**Input:**

1. Bulk rna data (gene in row, patient in column)
2. Clinical information (have overall survival (time), state )
3. Oncofetal signature

```bash
# 10. evaluate the clinical relevance of oncofetal signature
python OF-scRNA.py clinical -p ./parameters/parameters.txt
```

**output:** /6.Clinical/

**oncofetal-gsva-TCGA-KM-plot.pdf**: Kaplan-Meier curves for overall survival of HCC patients in the TCGA ( cohort stratified by GSVA score of mCAF onco-fetal signature. The p value was calculated by the log-rank test.

![image-20220512201319399](https://s2.loli.net/2022/05/12/a9s1lzxqguVkMUO.png)
</font>
