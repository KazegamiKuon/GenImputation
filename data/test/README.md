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
├── README.md
├── chr22_gtruth.log                - process log.
├── chr22_gtruth.recode.vcf         - ground truth for imputing. See below.
├── chr22_input.log                 - process log.
├── chr22_input.recode.vcf          - input for model. See below.
├── chr22_predict.dose.vcf.gz       - predicted by minimac4. See at description below.
├── chr22_predict.dose.vcf.gz.csi   - created by bcftools. See at description below.
├── chr22_predict.info              - created when predict "chr22_predict.dose.vcf.gz".
└── chr22_predict.logfile           - created when predict "chr22_predict.dose.vcf.gz".
```

## Description

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

**chr22_predict.dose.vcf.gz** file was created by this script at this folder:

```script
minimac4 --refHaps ../train/chr22.recode.m3vcf.gz --haps chr22_input.recode.vcf --allTypedSites --log --cpus 4 --prefix chr22_predict
```

**chr22_predict.dose.vcf.gz.csi** file was created by this script at this folder:

```script
bcftools index chr22_predict.dose.vcf.gz
```
