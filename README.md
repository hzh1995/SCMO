# SCMO
## Introduction
SCMO is a single-cell data-enhanced survival prediction model for colorectal cancer (CRC) that integrates multi-omics data. By leveraging CRC single-cell data as a reference, we identified CRC-specific tumor microenvironment (TME) characteristics. Subsequently, using these TME features as a reference, we reconstructed the TME profiles for each patient from bulk RNA-seq data in the TCGA-CRC cohort. We further developed a survival prediction model using a Siamese Neural Network (SNN), incorporating TME features alongside RNA, DNA, and microbiome characteristics. This integrative approach enables a comprehensive assessment of CRC biology and improves prognostic accuracy through multi-omics feature fusion.

In this work, we integrated the Ecotyper algorithm with the Siamese Neural Network (SNN) model from PORPOISE.



## Step1 Discovery of Cell States and Ecotypes in scRNA-seq CRC Data
In this step of the analysis, single-cell data derived from Colorectal Cancer (CRC) was utilized as a reference. This allowed us to discover ecotypes and cell states that are specifically associated with the tumor microenvironment in CRC.

The script that does cell type and ecotype discovery is:
```
Rscript EcoTyper_discovery_bulk.R -c config_discovery_scRNA.yml
```
```
## usage: EcoTyper_discovery_scRNA.R [-c <PATH>] [-h]
## 
## Arguments:
##   -c <PATH>, --config <PATH>
##                         Path to the config files [required].
##   -h, --help            Print help message.
```
This script takes as input file a configuration file in YAML format. The configuration file for this step can be referenced  in :config_discovery_scRNA.yml
```
default :
  Input :    
    Discovery dataset name : "DiscoveryOutput_CRC_sample_detaild_annotation"
    Expression matrix : "/SC_CRC/CRC_sample/scRNA_object_sample_new.txt"    
    Annotation file : "/SC_CRC/CRC_sample/scRNA_object_sample_new_Annotation.txt" 
    Annotation file column to scale by : NULL
    Annotation file column(s) to plot : []
    
  Output :
    Output folder : "DiscoveryOutput_scRNA"

  Pipeline settings :
    #Pipeline steps:
    #   step 1 (extract cell type specific genes)
    #   step 2 (cell state discovery on correrlation matrices)
    #   step 3 (choosing the number of cell states)
    #   step 4 (extracting cell state information)
    #   step 5 (cell state re-discovery in expression matrices)
    #   step 6 (extracting information for re-discovered cell states)
    #   step 7 (cell state QC filter)
    #   step 8 (ecotype discovery)
    Pipeline steps to skip : [] 
    # Accepted values: 
    # "cell type specific" - select genes overexpressed in a cell type
    # <integer> - e.g. 1000, select top <integer> genes with highest variance in a cell type
    # "no filter" - use all genes
    Filter genes : "cell type specific"
    Number of threads : 140
    Number of NMF restarts : 50
    Maximum number of states per cell type : 50
    Cophenetic coefficient cutoff : 0.95
    #The p-value cutoff used for filtering non-significant overlaps in the jaccard matrix used for discovering ecotypes in step 8. Default: 1 (no filtering).
    Jaccard matrix p-value cutoff : 0.001
    Minimum number of states in ecotypes : 3
```
Due to file size limitations, the two files scRNA_object_sample_new.txt and scRNA_object_sample_new_Annotation.txt could not be uploaded.

**Input**

**Discovery dataset name**
The name of the discovery dataset is user-defined. If users want to recover ecotypes and  cell states, they can utilize this dataset by specifying it with the -d/--discovery parameter. This enables the incorporation of user-discovered ecotypes and cell states as references for the analysis.
```
Discovery dataset name : "DiscoveryOutput_CRC_sample_detaild_annotation"
```

**Expression matrix**
```
Expression matrix : "/SC_CRC/CRC_sample/scRNA_object_sample_new.txt"
```
Expression matrix field should contain the path to a tab-delimited file containing the expression data, with genes as rows and cells as columns. The expression matrix should be in the TPM, CPM or other suitable normalized space. The specific format is as follows:
```
data = read.delim("/SC_CRC/CRC_sample/scRNA_object_sample_new.txt", nrow = 5)
head(data[,1:5])
```
```
    Gene SMC01.T_AACTTTCGTTCGGGCT SMC01.T_ACAGCCGTCCCTAATT SMC01.T_ACGGAGAAGGATGGAA
1   A1BG                        0                        0                        0
2   A1CF                        0                        0                        0
3    A2M                        0                        0                        0
4  A2ML1                        0                        0                        0
5 A4GALT                        0                        0                        0
  SMC01.T_ACTGAGTGTTGACGTT
1                        0
2                        0
3                        0
4                        0
5                        0
```
**Annotation file**
```
Annotation file : "/SC_CRC/CRC_sample/scRNA_object_sample_new_Annotation.txt" 
```
The Annotation file must contain the Celltyoe column, other information can be added as required.The specific format is as follows:
```
annotation = read.delim("/SC_CRC/CRC_sample/scRNA_object_sample_new_Annotation.txt", nrow = 5)
head(annotation)
```
```
                        ID  Sample         CellType
1 SMC01.T_AACTTTCGTTCGGGCT SMC01.T Epithelial.cells
2 SMC01.T_ACAGCCGTCCCTAATT SMC01.T Epithelial.cells
3 SMC01.T_ACGGAGAAGGATGGAA SMC01.T Epithelial.cells
4 SMC01.T_ACTGAGTGTTGACGTT SMC01.T Epithelial.cells
5 SMC01.T_ACTTTCATCCTTTCTC SMC01.T Epithelial.cells
```
**Annotation file column to scale by**
```
Annotation file column to scale by : NULL
```
This parameter can be set to NULL if the user only focuses on one type of cancer analysis.
In our study, this field is not used and therefore set to NULL.

**The output section**
The Output section contains a single field, Output folder, which specifies the path where the final output will be saved. This folder will be created if it does not exist.
```
DiscoveryOutput_scRNA
```
**Pipeline settings**

**Pipeline steps to skip**
```
#Pipeline steps:
    #   step 1 (extract cell type specific genes)
    #   step 2 (cell state discovery on correlation matrices)
    #   step 3 (choosing the number of cell states)
    #   step 4 (extracting cell state information)
    #   step 5 (cell state re-discovery in expression matrices)
    #   step 6 (extracting information for re-discovered cell states)
    #   step 7 (cell state QC filter)
    #   step 8 (ecotype discovery)
    Pipeline steps to skip : [] 
```
If the user needs to adjust the parameters, then the user can select the analysis steps to be skipped.For example, if the Cophenetic coefficient cutoff used in step 3 needs adjusting, the user might want to skip steps 1-2 and re-run from step 3 onwards.

**Filter non cell type specific genes**
```
# Accepted values: 
    # "cell type specific" - select genes overexpressed in a cell type
    # <integer> - e.g. 1000, select top <integer> genes with highest variance in a cell type
    # "no filter" - use all genes
    Filter genes : "cell type specific"
```
The option "cell type specific" will use the genes that are specific to each cell type. The option <integer> will select the top number of genes with the highest variance within each cell type. The option "no filter" will use all genes without applying any filtering.

**Number of threads**
```
Number of threads : 120
```
Specify the number of threads to be used for analysis, depending on the server's performance as decided by the user.

**Number of NMF restarts**
```
Number of NMF restarts : 50
```
To obtain a stable solution, NMF is generally run multiple times with different seeds, and the solution that best explains the discovery data is chosen. Additionally, the variation of NMF solutions across restarts with different seeds is quantified using Cophenetic coefficients and used in step 4 of EcoTyper for selecting the number of states. The parameter Number of NMF restarts specifies how many restarts with different seed should EcoTyper perform for each rank selection, in each cell type. Since this is a very time consuming process, in this example we only use 5. However, for publication-quality results, we recommend at least 50 restarts.

**Maximum number of states per cell type**
```
Maximum number of states per cell type : 20
```
Maximum number of states per cell type specifies the upper end of the range for the number of states possible for each cell type. The lower end is 2.

**Cophenetic coefficient cutoff**
```
Cophenetic coefficient cutoff : 0.95
```
This field indicates the Cophenetic coefficient cutoff, in the range [0, 1], used for automatically determining the number of states in step 4. Lower values generally lead to more clusters being identified. In this particular example, we set it to 0.95.

**Jaccard matrix p-value cutoff**
```
Jaccard matrix p-value cutoff : 0.001
```
Ecotype identification on step 8 is performed by clustering a jaccard matrix that quantifies the sample overlap between each pair of states. Prior to performing ecotype identification, the jaccard matrix values corresponding to pairs of states for which the sample overlap is not significant are set to 0, in order to mitigate the noise introduced by spurious overlaps. The value provided in this field specifies the p-value cutoff above which the overlaps are considered non-significant. When the number of samples in the scRNA-seq dataset is small, such as in the current example, we recommend this filter is disabled (p-value cutoff = 1), to avoid over-filtering the jaccard matrix. However, we encourage users to set this cutoff to lower values (e.g. 0.05), if the discovery scRNA-seq dataset contains a number of samples large enough to reliably evaluate the significance of overlaps.

**Minimum number of states in ecotypes**
```
Minimum number of states in ecotypes : 3
```
The ecotypes with less cell states than indicated in this field will be filtered out.
## Step2 Recovery of Cell States and Ecotypes in TCGA-CRC Bulk Data
For this step, we used  the TCGA bulk samples from CRC.
The script used to perform recovery in bulk data is :
```
Rscript EcoTyper_recovery_bulk.R -h
```
```
## usage: EcoTyper_recovery_bulk.R [-d <character>] [-m <PATH>] [-a <PATH>]
##                                 [-c <character>] [-t <integer>] [-o <PATH>]
##                                 [-h]
## 
## Arguments:
##   -d <character>, --discovery <character>
##                         The name of the discovery dataset used to define cell
##                         states and ecotypes. Accepted values: 'Carcinoma' will
##                         recover the cell states and ecotypes defined across
##                         carcinomas, as described in the EcoTyper carcinoma
##                         paper, 'Lymphoma' will recover the cell states and
##                         ecotypes defined in diffuse large B cell lymphoma
##                         (DLBCL), as described in the EcoTyper lymphoma paper,
##                         '<MyDiscovery>' the value used in the field 'Discovery
##                         dataset name' of the config file used for running
##                         EcoTyper discovery ('EcoTyper_discovery.R') script.
##                         [default: 'Carcinoma']
##   -m <PATH>, --matrix <PATH>
##                         Path to a tab-delimited file containing the input bulk
##                         tissue expression matrix, with gene names on the first
##                         column and sample ids as column names [required].
##   -a <PATH>, --annotation <PATH>
##                         Path to a tab-delimited annotation file containing the
##                         annotation of samples in the input matrix. This file
##                         has to contain in column 'ID' the same ids used as
##                         column names in the input matrix, and any number of
##                         additional columns. The additional columns can be
##                         plotted as color bars in the output heatmaps.
##                         [default: 'NULL']
##   -c <character>, --columns <character>
##                         A comma-spearated list of column names from the
##                         annotation file to be plotted as color bar in the
##                         output heatmaps. [default: 'NULL']
##   -t <integer>, --threads <integer>
##                         Number of threads. [default: '10']
##   -o <PATH>, --output <PATH>
##                         Output directory path. [default: 'RecoveryOutput']
##   -h, --help            Print help message.
```
The script takes the following arguments:

-d/–discovery: The name of the discovery dataset used for defining cell states. The name of the discovery dataset is the value provided in the Discovery dataset name field of the configuration file used for running cell state discovery. In our tutorial, the name of the discovery dataset is DiscoveryOutput_CRC_sample_detaild_annotation.

-m/–matrix: Path to the input expression matrix. The expression matrix should be in the TPM or FPKM space for bulk RNA-seq and non-logarithmic (exponential) space for microarrays. It should have gene symbols on the first column and gene counts for each sample on the next columns. Column (sample) names should be unique. Also, we recommend that the column names do not contain special characters that are modified by the R function make.names, e.g. having digits at the beginning of the name or containing characters such as space, tab or -. The CRC cancer bulk data used in this step looks as follows:
```
data = read.delim("bulk_CRC_data.txt", nrow = 5)
head(data[,1:5])
```
```
      gene TCGA.3L.AA1B.01A TCGA.4N.A93T.01A TCGA.4T.AA8H.01A TCGA.5M.AAT4.01A
1   TSPAN6        5.2200268        5.1540641         4.498977         5.528208
2     TNMD        0.4035851        0.9348164         1.511457         1.107263
3     DPM1        5.3095168        5.6079264         4.616169         5.481517
4    SCYL3        1.9458830        1.9691536         1.485403         1.788333
5 C1orf112        1.4670772        0.9829993         1.138072         1.915354
```
-a/–annotation: Path to a tab-delimited annotation file (not required). If provided, this file should contain a column called ID with the same values as the columns of the expression matrix. Additionally, this file can contain any number of columns, that can be used for plotting as color bars in the output heatmaps (see argument -c/–columns).
```
annotation = read.delim("bulk_CRC_annotation.txt", nrow = 5)
head(annotation)
```
```
                ID OS OS.time
1 TCGA.3L.AA1B.01A  0   15.83
2 TCGA.4N.A93T.01A  0    4.87
3 TCGA.4T.AA8H.01A  0   12.83
4 TCGA.5M.AAT4.01A  1    1.63
5 TCGA.5M.AAT6.01A  1    9.67
```
-c/–columns: A comma-separated list of column names from the annotation file (see argument -a/–annotation) to be plotted as color bars in the output heatmaps. By default, the output heatmaps contain as color bar the cell state label each cell is assigned to. The column names indicated by this argument will be added to that color bar.

-t/–threads: Number of threads. Default: 10.

-o/–output: Output folder. The output folder will be created if it does not exist.

The command line for recovering the carcinoma cell states and ecotypes in our study bulk data is:
```
Rscript EcoTyper_recovery_bulk.R -d DiscoveryOutput_CRC_sample_detaild_annotation -m bulk_CRC_data.txt -a bulk_CRC_annotation.txt -o RecoveryOutput_CRC_detailed_p0.001cut0.95_idx0.2
```

## Step3 SCMO SNN model
**data**

The training and testing dataset used in this study can be found at the following paths.
```
PORPOISE-SNN/CRC_newtest/train/tcga_coadread_all_clean.csv
```
```
PORPOISE-SNN/CRC_newtest/test/tcga_coadread_all_clean.csv
```
**splits**

For evaluating the algorithm's performance, we randomly partitioned each dataset using 5-fold cross-validation. Splits for CRC are found at the following paths:
```
PORPOISE-SNN/splits/5fold_crc_train/tcga_coadread
```
```
PORPOISE-SNN/splits/5fold_crc_test/tcga_coadread
```
**Running command-line**

To run experiments using the SNN, AMIL, and MMF networks defined in this repository, experiments can be run using the following generic command-line:
```
CUDA_VISIBLE_DEVICES=<DEVICE ID> python main.py --which_splits <SPLIT FOLDER PATH> --split_dir <SPLITS FOR CANCER TYPE> --mode <WHICH MODALITY> --model_type <WHICH MODEL> --result_dir <RESULT PATH>
```
Specific instructions for using the parameters can be found in the file ```main.py```.

**Training set SNN**
```
CUDA_VISIBLE_DEVICES=0,1,2,3 python main.py --which_splits 5fold_crc_train --split_dir tcga_coadread --mode omic --reg_type omic --model_type snn --results_dir ./CRC_newtest/train/CRC_RNA_DNA_MIC_ECO
```
**Test set SNN**
```
CUDA_VISIBLE_DEVICES=0,1,2,3, python main.py --which_splits 5fold_crc_test --split_dir tcga_coadread --mode omic --reg_type omic --model_type snn --results_dir ./CRC_newtest/test/CRC_RNA_DNA_MIC_ECO
```