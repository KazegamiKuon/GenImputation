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

**Declare variable:**

```script
source="./data/interim/G1K_chr20_biallelic.vcf.gz"
gtruth="./data/test/G1K_chr20_biallelic_gtruth.vcf.gz"
input="./data/test/G1K_chr20_biallelic_input.vcf.gz"
dir_input="./data/test/G1K_chr20_biallelic_input"
dir_input_file="./data/test/G1K_chr20_biallelic_input/0000.vcf.gz"
genotype="This is path of genotype vcf file to exstract that position from gtruth"
sample="./data/external/test_100_samples.txt"
train="./data/train/G1K_chr20_biallelic_train.m3vcf.gz"
predict_prefix="./data/test/G1K_chr22_biallelic_predict"
```

**G1K_chr20_biallelic_gtruth.vcf.gz** file was created by this script at root project:

```script
bcftools view $source -o $gtruth -O z -S $sample
```

**G1K_chr20_biallelic_input.vcf.gz** file was created by this script at root project:

```script
bcftools index $gtruth
bcftools index $genotype
bcftools isec $gtruth $genotype -p $dir_input -n =2 -w 1 -O z
cp $dir_input_file $input
```

**G1K_chr20_biallelic_predict.dose.vcf.gz** file was created by this script at this folder by run at dockercontainer:

```script
minimac4 --refHaps $train --haps $input --allTypedSites --log --cpus 4 --prefix $predict_prefix
```