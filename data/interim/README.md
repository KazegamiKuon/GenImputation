# INTERIM

Data on intermediate format between raw and processed. Not raw and also not ready yet.

All description must expend how and why you create this also cite if it have

## Folder tree should have structure

Same as **Folder tree should have structure** at test folder.

## Truth structure

```tree
interim
├── G1K_chr20_hg38.hap.gz
├── G1K_chr20_hg38.legend.gz
├── G1K_chr20_hg38_true.hap.gz
├── G1K_chr20_hg38_true.legend.gz
├── G1K_chr20_hg38_true.sample
├── process.ipynb
├── G1K_chr22_biallelic.vcf.gz
├── README.md
└── VN_chr22_biallelic.vcf.gz
```

## Description

**Declare variable**

```
source_file="./data/raw/G1K_chr20_hg38.vcf.gz"
output_file="./data/interim/G1K_chr20_biallelic.vcf.gz"
```

**G1K_chr20_biallelic.vcf.gz** was created by remove multiallelic and duplicate position (keep first position). File was created by this script at root project:

```script
bcftools norm -m -both $source_file -O z -o $output_file
bcftools index $output_file
```

```
bcftools norm -m -both $source_file -O z -o $output_file
```

**G1K_chr20_hg38** file is data by run cell at **process.ipynb**
