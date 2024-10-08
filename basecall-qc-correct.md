## 测序数据的质控 | quality control of sequencing data

这个列表收集了短读长、长度长的质控工具。

## 数据模拟 | reads simulation
https://github.com/lh3/wgsim

### basecall

### qc
- [[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)] - [Java, 0.12.0, 2023.03] - [A quality control tool for high throughput sequence data.]
- [[fastp](https://github.com/OpenGene/fastp)] - [C++/C, v0.23.4, 2023.05] - [2018.09, _Bioinformatics_] - [[fastp: an ultra-fast all-in-one FASTQ preprocessor.](https://doi.org/10.1093/bioinformatics/bty560)]
- [[fastp](https://github.com/OpenGene/fastp)] - [C++/C, v0.23.4, 2023.05] - [2023.05, _iMeta_] -[[Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp.](https://doi.org/10.1002/imt2.107)]
- [[Trimmomatic](https://github.com/usadellab/Trimmomatic)] - [Java, v0.39, 2021.01] - [2014,04, _Bioinformatics_] - [[Trimmomatic: a flexible trimmer for Illumina sequence data.](https://doi.org/10.1093/bioinformatics/btu170)]
- [[Cutadapt](https://github.com/marcelm/cutadapt)] - [Python, v4.9, 2024.06] - [2010.na, _EMBnet J_] - [[Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads.](https://doi.org/10.14806/ej.17.1.200)]
- [[FastUniq](https://sourceforge.net/projects/fastuniq/files/)] - [C++, 1.1, 2012.05] - [2012.12, _Plos ONE_] - [[FastUniq: A Fast De Novo Duplicates Removal Tool for Paired Short Reads.](https://doi.org/10.1371/journal.pone.0052249)]
- [[PoreChop][https://github.com/rrwick/Porechop]] - [C++, v0.2.4, 2018.10] - [adapter trimmer for Oxford Nanopore reads. **Unsupported since Oct 2018**]
- [[Sickle](https://github.com/najoshi/sickle)] - [C, v1.33, 2014.7] - [Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files]
- [[PEAT](https://github.com/jhhung/PEAT)] - [C++, v1.2.4, 2016.2, No longer maintained, recommond EARRINGS] - [2015.1, _BMC Bioinformatics_] - [PEAT: an intelligent and efficient paired-end sequencing adapter trimming algorithm.](https://doi.org/10.1186/1471-2105-16-S1-S2)
- [[EARRINGS](https://github.com/jhhung/EARRINGS)] - [C++, v1.1.0, 2024.6] - [2021.7, _Bioinformatics_] - [[EARRINGS: an efficient and accurate adapter trimmer entails no a priori adapter sequences](https://doi.org/10.1093/bioinformatics/btab025)]
- [[NanoPack](https://github.com/wdecoster/nanopack)] - [Python, No release, ] - [2018.3, _Bioinformatics_] - [NanoPack: visualizing and processing long-read sequencing data.](https://doi.org/10.1093/bioinformatics/bty149)
- [NanoPlot](https://github.com/wdecoster/NanoPlot) - [Python, v1.42.0, 2023.10] - [2023.5, _Bioinformatics_] - [NanoPack2: population-scale evaluation of long-read sequencing data](https://doi.org/10.1093/bioinformatics/btad311)
- [chopper](https://github.com/wdecoster/chopper) - [Rust, v0.9.0, 2024.8] - [2023.5, _Bioinformatics_] - [[NanoPack2: population-scale evaluation of long-read sequencing data.](https://doi.org/10.1093/bioinformatics/btad311)]
- [nanoQC](https://github.com/wdecoster/nanoQC) - [Python, no release, na] - [2018.3, _Bioinformatics_] - [[NanoPack: visualizing and processing long-read sequencing data](https://doi.org/10.1093/bioinformatics/bty149)]

### correct

### 长读长一致序列 | consensus sequences from long-reads
- [[medaka](https://github.com/nanoporetech/medaka)] - [Python, v2.0.0, 2024.9] - [a tool to create consensus sequences and variant calls from nanopore sequencing data.]
