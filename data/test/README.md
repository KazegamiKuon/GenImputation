# TEST

Data use to test/evaluation

All description must expend how and why you create this also cite if it have

## Folder tree should have structure

```tree
test                        - this folder, any test data must be here
│   README.md               - description about this folder test
|
└───dataname                - this folder's name is your data name
│   │   data.txt            - data
|   |   README.md           - description about this data
│   │
│   └───subfolder           - subfolder when data div to batch/patch
│       │   data001.txt
│       │   ...
│   
└───...
```

## Truth structure

```tree
test
├── G1K_chr22_biallelic_gtruth.log              - ground truth for imputing. See below.
├── G1K_chr22_biallelic_gtruth.recode.vcf.gz
├── G1K_chr22_biallelic_input.log               - input for model. See below.
├── G1K_chr22_biallelic_input.recode.vcf.gz
├── G1K_chr22_biallelic_predict.dose.vcf.gz     - predicted by minimac4. See at description below.
├── G1K_chr22_biallelic_predict.info            - created when predict
├── G1K_chr22_biallelic_predict.logfile         - created when predict
└── README.md
```

## Description

**3 script below** create input and gtruth for minimac 4 impute

**G1K_chr22_biallelic_gtruth.recode.vcf.gz** file was created by this script at root project:

```script
vcftools --gzvcf ./data/interim/G1K_chr22_biallelic.vcf.gz --keep ./data/external/test_100_samples.txt --exclude-positions ./data/external/omni_chr22_position.csv --out ./data/test/G1K_chr22_biallelic_gtruth --recode --recode-INFO-all
gzip ./data/test/G1K_chr22_biallelic_gtruth.recode.vcf

```

**G1K_chr22_biallelic_input.recode.vcf.gz** file was created by this script at root project:

```script
vcftools --gzvcf ./data/interim/G1K_chr22_biallelic.vcf.gz --keep ./data/external/test_100_samples.txt --positions ./data/external/omni_chr22_position.csv --out ./data/test/G1K_chr22_biallelic_input --recode --recode-INFO-all
gzip ./data/test/G1K_chr22_biallelic_input.recode.vcf
```

**G1K_chr22_biallelic_predict.dose.vcf.gz** file was created by this script at this folder:

```script
minimac4 --refHaps ../train/G1K_chr22_biallelic_train.recode.m3vcf.gz --haps G1K_chr22_biallelic_input.recode.vcf.gz --allTypedSites --log --cpus 4 --prefix G1K_chr22_biallelic_predict
```