# TRAIN

Data use t train

All description must expend how and why you create this also cite if it have

## Folder tree should have structur

Same as **Folder tree should have structure** at test folder.

## Truth structure

```tree
train
├── G1K_chr22_biallelic_train.log                   - process log file.
├── G1K_chr22_biallelic_train.recode.erate          - was created when run script create "m3vcf" file.
├── G1K_chr22_biallelic_train.recode.m3vcf.gz       - was created by minimac3. See at description below.
├── G1K_chr22_biallelic_train.recode.rec            - was created when run script create "m3vcf" file.
├── G1K_chr22_biallelic_train.recode.vcf.gz         - data for training. See at description below.
└── README.md
```

## Description

**NOTE**: see **README.md** at data root folder to use minimac.

**To impute by minimac3** data should be biallelic. Then see below to prepare data.

**m3vcf.gz**:

```script
minimac3 --refHaps <data path to train> --processReference --prefix <ouput name prefix>
```

## For the impatient

**G1K_chr22_biallelic_train.recode.m3vcf.gz** file was created by this script at this folder:

```script
minimac3 --refHaps G1K_chr22_biallelic_train.recode.vcf.gz --processReference --prefix G1K_chr22_biallelic_train.recode
```
