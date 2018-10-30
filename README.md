#   Pilon-eval

## 0. Introduction

A bash wrapper script for Pilon.

## 1. Prerequisites

This script depends on:
	Java runtime 1.7 or later
	seqtk		(https://github.com/lh3/seqtk)
	Bowtie2		ver. >2.1	(http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
	Samtools	ver. >1.4	(https://github.com/samtools/samtools)

## 2. Installation

No need to compile since this is a bash wrapper script. Put this where you want to install. 
In the directory, change its permission to be executable.

```
git clone https://github.com/Yukihirokinjo/Pilon-eval.git
cd Pilon-eval

$ chmod u+x Pilon-eval.bash
```

Thereafter, add the path to your PATH.  
For example..
```
$ echo 'export PATH=/path/to/Pilon-eval_dir/:$PATH' >> ~/.bashrc
$ source ~/.bashrc
```

To check the installation, you can type
```
$ java -Xmx16G  -jar  $(which pilon-1.22.jar)
```
If you can see pilon usage, the instllation was successful.


## 3. Running 

### 3.1 Input data

Input contigs/scaffolds for this script must be in fasta format. 
As input read files for this script, Illumina paired-end reads are assumed.

### 3.2 Quick start

```
$ Pilon-eval.bash  -i input_scaffolds.fasta -1 reads_1.fastq -2 reads_2.fastq -o output_dir
```

### 3.3 Command line options
--------------------------------------------------------------------------------

	Mandatory:
	-i		<FILE>	Input contigs/scaffolds file.

	-1		<FILE>	Input read file (forward).

	-2		<FILE>	Input read file (reverse).


	Optional:
	-o		<str>	Output directory (default: "Pilon_eval_<current time>").

	-c		<int>	Number of threads to be used for computation (default: 1).

	-k		<int>	Size of K-mer for internal assembler in Pilon (default: 47). 

--------------------------------------------------------------------------------


### 3.4 Output directories/files

--------------------------------------------------------------------------------
	./EstInsSize		:This directory contains initial mapping results to be used for insert size estimation.

	./pilon.fasta		:Output contigs file generated by Pilon.

--------------------------------------------------------------------------------

