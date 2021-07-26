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

Data have chr22 or 22 mean #CHROM property is format chr**chrom name** or **chrom name**.

**Declare Variable chr22 hs37d5**

```script
data_source="./data/interim/to_biallelic/G1K_chr22_hs37d5_biallelic.vcf.gz"
train_prefix="./data/train/G1K_chr22_hs37d5_biallelic_train.vcf.gz"
test_sample="./data/external/test_100_samples.txt"
m3vcf_prefix="./data/train/G1K_chr22_hs37d5_biallelic_train"
```

```
chrm_map="data/raw/chrs_name_map_file_chr_to_num.txt"
train_prefix_num="./data/train/G1K_22_hs37d5_biallelic_train.vcf.gz"
```

**<data_name>_train.vcf.gz** file was created by this script at this folder:

```script
bcftools view $data_source -o $train_prefix -O z -S ^$test_sample
bcftools index $train_prefix
```

Connvert from chr**number** to **number**. 

```script
bcftools annotate $train_prefix --rename-chrs $chrm_map -o $train_prefix_num -O z
m3vcf_prefix="./data/train/minimac/G1K_22_hs37d5_biallelic_train"
train_prefix=$train_prefix_num
bcftools index $train_prefix
```

**<data_name>_train.m3vcf.gz** file was created by this script at this folder:

```script
minimac3 --refHaps $train_prefix --processReference --prefix $m3vcf_prefix
```
