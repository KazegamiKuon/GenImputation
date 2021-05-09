# DATA

This description about some file at this folder.

Ex:
process.ipynb will call lib or func at untils or lib/datasets folder to view, test, process ... data

## Folder tree

Use this cmd to show Folder tree.

```script
tree -I zarr
```

```tree
data                        - this folder, any data must be here
├── README.md           - description about this folder.
├── external            - any data we consider as being external.
├── interim             - data on intermediate format between raw and processed. Not raw and also not ready yet.
├── process.ipynb       - test/process/ ... we can call lib to process data at here.
├── raw                 - here all raw data exists. This directory should be considered as read only - just 
|                       leave what we got as it is.
├── test                - data use to test/evaluation.
├── train               - data use to train.
└── val                 - data use to validation/evaluation.
```

## Script/ cmd could usefully

Work with vcf file, [learning_vcf](https://github.com/davetang/learning_vcf_file) environment would helpfull. You should install it.

Preprocess data to biallelic (only one ALT each variant) data exclude duplicate position.

```script
bcftools view -m2 -M2 -v snps --collapse all -o <vcf_path to save> -O z <vcf_path to process>
```

**Description:** [common options](http://samtools.github.io/bcftools/bcftools.html#common_options) to see more option

```desciption
-m2 -M2 -v snps         option to only view biallelic SNPs.
-o                      output name.
-O z                    output type, z mean output type is vcf will compressed.
--collapse              controls how to treat records with duplicate positions 
                        and defines compatible records across multiple input file.
                        In the case of records with the same position, only the first will be considered and appear on output.
```

## Docker

Get minimac was ready below:

```script
docker pull thehung92phuyen/imputation:v5.2
```

Mount container to this folder data

```script
docker run --name minimac -v C://Users//cuongdq4//Documents//VBDI//GenImputation//:/home/GenImputation -it thehung92phuyen/imputation:v5.2

Detail:
docker run --name <container name> -v <absolute path without container>:<absolute path within container> -it thehung92phuyen/imputation:v5.2
```

To start and attach cotainer to terminal.

```script
docker start minimac
docker container attach minimac

Details:
docker start <container name>
docker container attach <container name>
```
