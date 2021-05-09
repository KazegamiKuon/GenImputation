# INTERIM

Data on intermediate format between raw and processed. Not raw and also not ready yet.

All description must expend how and why you create this also cite if it have

## Folder tree should have structure

Same as **Folder tree should have structure** at test folder.

## Truth structure

```tree
interim
├── G1K_chr22_biallelic.vcf.gz
└── README.md
```

## Description

**G1K_chr22_biallelic.vcf.gz** was created by remove multiallelic and duplicate position (keep first position). **learning_vcf** environment must activate. File was created by this script at root project:

```script
bcftools annotate --collapse all ./data/raw/G1K_chr22.vcf.gz | bcftools view -m2 -M2 -v snps -o ./data/interim/G1K_chr22_biallelic.vcf.gz -O z
```
