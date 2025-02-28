# SCMO
## Introduction
SCMO is a single-cell data-enhanced survival prediction model for colorectal cancer (CRC) that integrates multi-omics data. By leveraging CRC single-cell data as a reference, we identified CRC-specific tumor microenvironment (TME) characteristics. Subsequently, using these TME features as a reference, we reconstructed the TME profiles for each patient from bulk RNA-seq data in the TCGA-CRC cohort. We further developed a survival prediction model using a Siamese Neural Network (SNN), incorporating TME features alongside RNA, DNA, and microbiome characteristics. This integrative approach enables a comprehensive assessment of CRC biology and improves prognostic accuracy through multi-omics feature fusion.

In this work, we integrated the Ecotyper algorithm with the Siamese Neural Network (SNN) model from PORPOISE.

**Ecotyper**(https://github.com/digitalcytometry/ecotyper)

**PORPOISE-SNN**(https://github.com/mahmoodlab/PORPOISE)

## Step1 Discovery of Cell States and Ecotypes in scRNA-seq CRC Data
The configuration file we used is referenced at
/Ecotyper/config_used.yml
```
Rscript EcoTyper_discovery_bulk.R -c config_discovery_scRNA.yml
```

## Step2 Recovery of Cell States and Ecotypes in TCGA-CRC Bulk Data

```
Rscript EcoTyper_recovery_bulk.R -d DiscoveryOutput_CRC_sample_detaild_annotation -m bulk_CRC_data.txt -a bulk_CRC_annotation.txt -o RecoveryOutput_CRC_detailed_p0.001cut0.95_idx0.2
```

## Step3 SCMO SNN model

```
python main.py --which_splits 5fold_crc_train --split_dir tcga_coadread --mode omic --reg_type omic --model_type snn --results_dir ./CRC_newtest/train/CRC_RNA_DNA_MIC_ECO
```
```
python main.py --which_splits 5fold_crc_test --split_dir tcga_coadread --mode omic --reg_type omic --model_type snn --results_dir ./CRC_newtest/test/CRC_RNA_DNA_MIC_ECO
```
