
## 1. Check the quality of the raw sequencing data (QC report) using:
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)
- [MultiQC](https://multiqc.info/docs/)

## 2. Pre-processing
1. Remove reads containing adapters using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html)
2. Remove reads containing N > 10% (N represents base that could not be determined)
3. Remove reads where the Qscore (Quality value) of over 90% bases of the read is <= 10

## 3. Align reads to a reference using [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) 


3. Count the number of reads assigned to each contig/gene
4. Extract counts and store in a matrix
5. Create column metadata (sampleinfo) table
6. Analyze count data using DESEQ2



