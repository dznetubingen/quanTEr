# quanTEr
Nextflow pipeline for QC, mapping and quantificaton of bulk RNAseq on genic and transposon levels.
Prerequisits:
* nextflow
* docker


To get the docker image with necessary tools installed run:

```
docker pull savytskanatalia/quantefication2
```
Examples of commands to run:

```
nextflow run fastp.nf -c fastp.config --reads '/absolute/path/to/files/*_R{1,2}.fastq.gz' -with-docker savytskanatalia/quantefication2
nextflow run quanTEr.nf -c quanTEr.config --reads "/absolute/path/to/files/*_R{1,2}.fastq.gz" ... -with-docker savytskanatalia/quantefication2
```




The pipeline consists of two scrpts. 

fastp.nf - performs QC of sequencing data using fastp, aggregating reports with MultiQC and uses pecheck to confirm synchronicity of reads in paired-end data.

quanTEr.nf or quanTEr2.nf - performs alignment of reads to reference genome using STAR, then uses TElocal to quantify transposons and genes from the data. STAR index can be built, if not supplied directly. For TElocal index file download one of the prebuilt indices distributed by Hammel's lab and change file extension to ".ind" for compatibility with older TElocal version that this pipeline uses (e.g. with command 'mv downloaded.uncompressed.file file.ind'). Note, in the latter case you need to build STAR index using same version of genome assembly as was used for index creation.


merge_count_tables.R - set of R functions to merge separate sample-individual TElocal count tables into three common tables, TE loci, TE subfamilies and genes.

For DZNE members: there are prebuilt custom TElocal index and a number of pre-built STAR-indices for human genome (hg38) for reads of varying length.





Further information on TElocal can be found on TElocal`s authors git:
https://github.com/mhammell-laboratory/TElocal

Further information on fastp...:
https://github.com/OpenGene/fastp



Savytska Natalia, DZNE Tuebingen, 2021
