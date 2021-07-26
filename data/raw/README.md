# RAW

Here all raw data exists. This directory should be considered as read only - just leave what we got as it is.

All description must expend how and why you create this also cite if it have

## Folder tree should have structure

Same as **Folder tree should have structure** at test folder.

## Truth structure

```tree
raw
├── G1K_chr22_hg38.vcf.gz
├── G1K_chr22_hg38.vcf.idx
├── README.md
├── VN_chr22.vcf.gz
├── hg38.fa.gz
└── infiniumomni2-5-8-v1-3-a2.csv.gz
```

## Description

Change dir to this folder then run script

**infiniumomni2-5-8-v1-3-a2.csv.gz** from [Omni kit](https://support.illumina.com/array/array_kits/humanomni2_5-8_beadchip_kit/downloads.html) have [version 1.3](ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/humanomni25/v1-3/infiniumomni2-5-8-v1-3-a2-manifest-file-csv.zip) support **hg38** (last seen 5/21/2021).

**Note**: When download successfull, **change that compressed file to .gz type**.

For the impatient:

```script
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/humanomni25/v1-3/infiniumomni2-5-8-v1-3-a2-manifest-file-csv.zip
unzip -p infiniumomni2-5-8-v1-3-a2-manifest-file-csv.zip | gzip -c > infiniumomni2-5-8-v1-3-a2.csv.gz
rm infiniumomni2-5-8-v1-3-a2-manifest-file-csv.zip
```

**hg38.fa.gz** from [reference genome hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

```script
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

**G1K_chr22_hg38.vcf.gz** from  [data project 1000 genomes](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/). It is vcf about chr22, download it and change it name to **G1K_chr22_hg38.vcf.gz**

```script
wget <link download file> -O G1K_chr22_hg38.vcf.gz
```

For the impatient:

```script
## This for chrom 22
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -O G1K_chr22_hg38.vcf.gz
```

```script
## This for chrom 20
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -O G1K_chr20_hg38.vcf.gz
```

**Note**: you can also down it from origin source but phease check **Sub data** below to see more.

**VN_chr22.vcf.gz** privated data, have this or not still okie!
**VN_chr20.vcf.gz** privated data, have this or not still okie!

## Sub data

**hg19.fa.gz** from [reference genome hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)

```script
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O ./data/raw/hg19.fa.gz
```

**G1K_chr22_hg38.vcf.gz** can liftover from **G1K_chr22_hs37d5.vcf.gz**.

**G1K_chr22_hs37d5.vcf.gz** from [Data project 1000 genomes](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). Download file **ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz**[1] or phase3 vcf file about chr22 if file [1] not exits and change it name to **G1K_chr22_hs37d5.vcf.gz**. Or run this script:

```script
wget <link download file> -O G1K_chr22_hs37d5.vcf.gz
bcftools annotate G1K_chr22_hs37d5.vcf.gz --rename-chrs chrs_name_map_file.txt -o G1K_chr22_hs37d5_v2.vcf.gz -O z
mv -f G1K_chr22_hs37d5_v2.vcf.gz G1K_chr22_hs37d5.vcf.gz
```

For the impatient:

```script
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -O G1K_chr22_hs37d5.vcf.gz
bcftools annotate G1K_chr22_hs37d5.vcf.gz --rename-chrs chrs_name_map_file.txt -o G1K_chr22_hs37d5_v2.vcf.gz -O z
mv -f G1K_chr22_hs37d5_v2.vcf.gz G1K_chr22_hs37d5.vcf.gz
```

**Note**: bcftools script use for change chrom name from 22 to chr22

**infiniumomni2-5-8v1-5-a1.csv.gz** from [Omni kit](https://support.illumina.com/array/array_kits/humanomni2_5-8_beadchip_kit/downloads.html) [version 1.5](https://webdata.illumina.com/downloads/productfiles/humanomni25/v1-5/infinium-omni2-5-8v1-5-a1-manifest-file-csv.zip)

For the impatient:

```script
wget https://webdata.illumina.com/downloads/productfiles/humanomni25/v1-5/infinium-omni2-5-8v1-5-a1-manifest-file-csv.zip
unzip -p infinium-omni2-5-8v1-5-a1-manifest-file-csv.zip | gzip -c > infiniumomni2-5-8v1-5-a1.csv.gz
rm infinium-omni2-5-8v1-5-a1-manifest-file-csv.zip
```

To view number variants at file vcf:

```script
bcftools view -H <vcf file> | wc -l
```

## No more avaliable

This data have reference genome is **hs37d5** based on NCBI **GRCh37**. But we study up to date use latest reference genome **hg38** so we need **liftover** from hs37d5 to hg38.

**Note**: **picard** from gatk should installed. Unless, run this script:

```script
conda install -c bioconda picard
```

To convert we should down chain file, we can find it [here](https://hgdownload.soe.ucsc.edu/downloads.html). Each reference genome (exclude hs37d5) have chain file to other ref. Bcuz it ref base on GRCh37/hg19 so we can use [chain file h19 to hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz) or [b37 to hg38](https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/funcotator/data_sources/gnomAD/b37ToHg38.over.chain)

Liftover by run this script. [Chain file h19 to hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz) was included on scipt at firt line.

```script
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
picard CreateSequenceDictionary --REFERENCE ./hg38.fa.gz
picard -Xmx6g LiftoverVcf -CHAIN hg19ToHg38.over.chain -INPUT G1K_chr22_hs37d5.vcf.gz -OUTPUT G1K_chr22_hg38.vcf.gz -REFERENCE_SEQUENCE hg38.fa.gz -REJECT G1K_chr22_hs37d5Tohg38_refect.vcf.gz
```
