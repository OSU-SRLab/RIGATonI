# RIGATONI
## 1.0: Installation
To install RIGATONI please use the code below.
```{r}
devtools::install_github('OSU-SRLab/RIGATONI')
```
For detailed information regarding the functions and steps described below, please refer to the RIGATonI Manuscript at https://www.biorxiv.org/content/10.1101/2024.03.02.583103v1

R version 4.1 or greater is required.

## 2.0: Run RIGATONI step by step
### 2.1: Make Gene List
RIGATONI helps you assess whether there is a functional change in a gene based on expression data. To get started, you should create two gene lists to assess the function of your gene of interest. You can either make them manually or use the function makeGeneList() below.
```{r}
makeGeneList('TP53')
```
makeGeneList takes an argument of a gene name in the form of a string (HGNC symbols) and uses the STRING database as default to create two gene lists. We recommend using the makeGeneList unless you have detailed and extensive knowledge of the gene of the interest or would like to assess the impact of the gene of interest on non-canonical pathways or genes.
If you would perfer to make your own upstream and downstream gene lists, you can using the example code below as a template. All genes should be written in capital letters.
Upstream genes refer to genes which act on the gene of interest and affect the expression of the gene of interest directly.
Downstream genes refer to genes on which the gene of interest acts in any capacity.
```{r}
gene_list = list()
gene_list['upstream'] <- c('MY', 'UPSTREAM', 'GENES')
gene_list['downstream'] <- c('MY', 'DOWNSTREAM', 'GENES')
```
Make sure your genes are in the same format as your gene expression file.
### 2.2: Predicting Function
To predict function, you first have to build regression models using WT or control samples. Next, you should apply those regression models to the experimental or mutant samples. Use the example code below as a template.
```{r}
Regression <- buildRegression(ControlRNA, gene_of_interest, gene_list)
Function <- getMutantFunction(Regression, MutantRNA, gene_of_interest)
```
ControlRNA and MutantRNA should be dataframes with rownames as gene names and colunm names as sample IDs. The gene_list can be made either through makeGeneList() or manually (see 'Make Gene List' for more information).
### 2.3: Predicting Immune phenotype
Immune phenotype predictions are provided as a probability that a given sample is "hot" or inflammed. To run this function, all that is needed is a dataframe of gene counts where rownames are gene names (HGNC symbols) and column names are sample IDs. Please follow the template shown below.
```{r}
getImmuneProb(MutantRNA)
```
## 3.0 Run RIGATONI with default settings all at once
You can also run RIGATONI with a single function and default settings using the template below.
```{r}
runRIGATONI(gene_of_interest, ControlRNA, MutantRNA)
```
The gene of interest is an HGNC symbol in the form of a character string. ControlRNA and MutantRNA should be dataframes with gene names for rownames and sample IDs for column names. We expect the gene names to be the same for both dataframes, and we expect the gene names to be in the HGNC symbol format.
