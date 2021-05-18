# RAW

Here all raw data exists. This directory should be considered as read only - just leave what we got as it is.

All description must expend how and why you create this also cite if it have

## Folder tree should have structure

Same as **Folder tree should have structure** at test folder.

## Truth structure

```tree
raw
├── G1K_chr22.vcf.gz
├── README.md
├── hg19.fa.gz
├── hg38.fa.gz
├── infinium-omni2-5-8v1-5-a1-manifest-file-csv.gz
└── VN_chr22.vcf.gz
```

## Description

**G1K_chr22.vcf.gz** from [Data project 1000 genomes](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). Download file **ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz**[1] or phase3 vcf file about chr22 if file [1] not exits and change it name to **G1K_chr22.vcf.gz**. Or run this script:

```script
wget <link download file> -O G1K_chr22.vcf.gz
```

**infinium-omni2-5-8v1-5-a1-manifest-file-csv.gz** from [Omni kit](https://support.illumina.com/array/array_kits/humanomni2_5-8_beadchip_kit/downloads.html). Download omni file csv format from **Infinium Omni2.5-8 v1.5 Product Files** which description is "Manifest, cluster, and LIMS product descriptor files for the Infinium Omni2.5-8 v1.5 kit" and **change that compressed file to .gz type**.

**hg19.fa.gz** from [reference genome hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and **hg38.fa.gz** **hg19.fa.gz** from [reference genome hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

**VN_chr22.vcf.gz** privated data, have this or not still okie!
