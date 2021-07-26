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

Script should run at root project

**Declare variable for chrom 20**

```
source_file="./data/raw/G1K_chr20_hg38.vcf.gz"
output_file="./data/interim/to_biallelic/G1K_chr20_hg38_biallelic.vcf.gz"
```

**Declare variable for chrom 20**

```
source_file="./data/raw/G1K_chr22_hs37d5.vcf.gz"
output_file="./data/interim/to_biallelic/G1K_chr22_hs37d5_biallelic.vcf.gz"
target_file="./data/interim/omni_bcftools_targets_22.targets.tsv.gz"
chrm_map="data/raw/chrs_name_map_file_chr_to_num.txt"
output_file_num="./data/interim/to_biallelic/G1K_22_hs37d5_biallelic.vcf.gz"
```

**G1K_<chrom>_biallelic.vcf.gz** was created by remove multiallelic and duplicate position (keep first position). File was created by this script:

```script
bcftools norm -m -both $source_file -O z -o $output_file
bcftools index $output_file
```

**omni_bcftools_targets_22.targets.gz** file is data by run cell at **process.ipynb** and run some script:

```script
bcftools annotate $output_file --rename-chrs $chrm_map -o $output_file_num -O z
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $output_file_num | bgzip -c > $target_file

### Run cell to subset line that data in omni chip

tabix -s1 -b2 -e2 $target_file
```

**G1K_chr20_hg38** file is data by run cell at **process.ipynb**
