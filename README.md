[![Build Status](https://travis-ci.org/yhoogstrate/flaimapper.svg?branch=master)](https://travis-ci.org/yhoogstrate/flaimapper)

# Fragment Location Annotation Identification Mapper #
---
	 ___    _      _____  _         _____  ___    ___    ___    ___         ___
	(  _`\ ( )    (  _  )(_)/'\_/`\(  _  )(  _`\ (  _`\ (  _`\ |  _`\     /'_  )
	| (_(_)| |    | (_) || ||     || (_) || |_) )| |_) )| (_(_)| (_) )   (_)_) |
	|  _)  | |  _ |  _  || || (_) ||  _  || ,__/'| ,__/'|  _)_ | ,  /     _(_ <
	| |    | |_( )| | | || || | | || | | || |    | |    | (_( )| |\ \    ( )_) |
	(_)    (____/'(_) (_)(_)(_) (_)(_) (_)(_)    (_)    (____/'(_) (_)   `\____)
---
![FlaiMapper](https://github.com/yhoogstrate/flaimapper/raw/master/share/flaimapper.png)
---

- [Download & Installation](#download--installation)
    - [From source (native)](#from-source-native)
    - [Install via bioconda](#install-via-bioconda)
    - [Read the GPL3 free software license](#read-the-gpl3-free-software-license)
    - [Check if it works!](#check-if-it-works)
    - [Uninstall FlaiMapper](#uninstall-flaimapper)
    - [Galaxy](#galaxy)
- [Alignment](#alignment)
    - [Obtain reference](#obtain-reference)
         - [Full genome alignment](#full-genome-alignment)
         - [ncRNAdb09 alignment](#ncrnadb09-alignment)
         - [Combination](#combination)
    - [Trim adapters](#trim-adapters)
         - [Quality trimming](#quality-trimming)
    - [Choose aligner](#choose-aligner)
         - [Full genome alignment](#full-genome-alignment-1)
         - [ncRNAdb09 alignment](#ncrnadb09-alignment-1)
    - [Read collapsing](#read-collapsing)
    - [Multi-mapping](#multi-mapping)
    - [Alignment indexing](#alignment-indexing)
- [Run FlaiMapper-3](#run-flaimapper-3)
    - [Usage](#usage)
    - [Input: BAM](#input-bam)
         - [The "\-\-parameters"-argument](#the---parameters-argument)
         - [The "\-\-fasta"-argument](#the---fasta-argument)
    - [Output: formats](#output-formats)
- [Reproduce article data](#reproduce-article-data)
- [Authors & Citing](#authors--citing)

## Download & Installation
### From source (native) ####
To get the latest version of FlaiMapper-3 please download with Git (make sure it is installed) from GitHub using the following terminal command:

	git clone https://github.com/yhoogstrate/flaimapper.git

You can install FlaiMapper-3 with corresponding dependencies via pip with the following terminal command:

	cd flaimapper
	cd src
	
	sudo pip install .

If you are not admin or would rather like have a contained version of FlaiMapper-3 you can use a virtual environment as well:

	cd flaimapper
	cd src
	
	virtualenv -p python3 .
	source activate .venv/activate/bin
	
	python setup.py build
	python setup.py install

### Install via bioconda

For any installed instance of (ana/mini/)conda, make sure you have the bioconda repository added:

	conda config --add channels bioconda
	conda config --add channels conda-forge
	
	# you can confirm this usually by reading the following file:
	cat ~/.condarc
	
	# setup an environment with the latest version of flaimapper
	conda create -n flaimapper flaimapper
	source activate conda

### Read the GPL3 free software license

If you have downloaded FlaiMapper-3, you automatically agree with the GNU General Public License v3.0. Please read the license to make sure you understand what you are allowed to do with it, and what free software means. It also states that the software is distributed in the hope to be useful, but without any warranty. The license can be accessed directly at the following url:

*	[https://github.com/yhoogstrate/flaimapper/raw/master/LICENSE](https://github.com/yhoogstrate/flaimapper/raw/master/LICENSE)

### Check if it works!

Last but not least, please check if FlaiMapper-3 runs. You can do this by running either a full analysis or just run the following terminal command: 

	flaimapper --help
	
	# If you have installed nose, you can automatically run the same test as Travis does
	# Make sure you are in ./flaimapper/src and run:
	nose

**If FlaiMapper-3 givers errors**, warnings or does not install, please do not hazitate and either **submit the bug** or send a fix to the GitHub repository at the following url: <A HREF='https://github.com/yhoogstrate/flaimapper' TARGET='_new'>https://github.com/yhoogstrate/flaimapper</A>.

If you did not experience an error: congratulations, you have just installed FlaiMapper-3!

### Uninstall FlaiMapper

To remove FlaiMapper-3 from your system directories (without removing your data files), proceed with the following terminal command:

	# Depending on whether you used su or sudo to install flaimapper:
	sudo pip uninstall flaimapper
	pip uninstall flaimapper
	
	# If you installed via virtualenv:
	cd ~/flaimapper/src
	rm -rf .venv
	
	# via conda
	conda env remove -n flaimapper

The downloaded source files can only be removed manually. Be aware that references files for *ncRNAdb09* and the data of the analysis as demonstrated in the corresponding article are also located within the FlaiMapper-3 root directory.

### Galaxy

If you think installing, configuring FlaiMapper-3, creating references and understanding the commandline interface takes too long, or you just want to have a visual interface, you can make use of Galaxy. We have a public instance available at the following url:

*	[http://bioinf-galaxian.erasmusmc.nl/galaxy/](http://bioinf-galaxian.erasmusmc.nl/galaxy/)

If you have a personal Galaxy instance, you can install FlaiMapper-3 via the toolshed at the following url:

*	[https://toolshed.g2.bx.psu.edu/view/yhoogstrate/flaimapper](https://toolshed.g2.bx.psu.edu/view/yhoogstrate/flaimapper)

You can find the Galaxy wrapper for FlaiMapper-3 at the following url:

*	[https://github.com/galaxyproject/tools-iuc/tree/master/tools/flaimapper/](https://github.com/galaxyproject/tools-iuc/tree/master/tools/flaimapper/)

## Alignment

A presentation on this topic available over here: [https://humgenprojects.lumc.nl/trac/humgenprojects/raw-attachment/wiki/RNA-seq-course/small_RNA_BioSB_RNA-Seq_2016.pdf](https://humgenprojects.lumc.nl/trac/humgenprojects/raw-attachment/wiki/RNA-seq-course/small_RNA_BioSB_RNA-Seq_2016.pdf)

### Obtain reference

There are multiple strategies for alignment in small RNA-Seq analysis, leaving it up to each researcher to decide which strategy fits best.

#### Full genome alignment

In this strategy reads are aligned to an enitre reference genome (e.g hg19).

*	Pro's: because you take into account the entire reference genome you will also align to previously unannotated ncRNAs (those that are not present in e.g. ncRNAdb09).
*	Con's: because small RNA-seq reads are relatively small, the chance to finding a read somewhere else in the genome by chance is relatively large. In such a case, reads may also align to places in the genome where they are not derived from. Some ncRNAs undergo maturation (splicing or ligation). This requires more complex alignment strategies.

#### ncRNAdb09 alignment

In this strategy reads are aligned to a list of annotated small ncRNAs with HUGO nomenclature (http://www.genenames.org/rna).
This list contains pre-miRNAs, snoRNAs, scRNAs, snRNAs, SNARs, vaultRNAs, miscRNAs and Y-RNAs, all with 10bp genomic extensions on both sides.
Also, the reference contains mature tRNAs from the genomic tRNA database (http://gtrnadb.ucsc.edu/) without introns and with CCA suffix.
The reference does not include Piwi- and rRNAs.
All identical sequences have been merged into single entries, but sequences may still partially overlap with others.

If you use this reference instead of a full reference genome, you have to be aware of the following:

*	Pro's: alignment is simpler (and faster) because it should not be aware of splicing (if mature ncRNAs are included).
*	Con's: you restrict yourself to a limited part of the genome and therefore you will miss any ncRNA that is not within this database.

#### Combination

In the first phase, align to targeted regions withing the genome (the ncRNAdb09 for example). In the next phase, align all unmapped reads as a classical RNA-Seq experiment to the reference genome.

*	Pro's: you solve the issues addressed above.
*	Con's: the alignment will take more time the and methodology is more complex which will probably require advanced scripting.

Some aligners have a related feature implemented: an indexed transcriptome besides the indexed reference genome. You might think of tophat's "<CODE>--transcriptome-index</CODE>" and RNA-STAR's "<CODE>sjdbGTFfile</CODE>".
However, idealy you want to have a three phase alignment: (1) ncRNAs, (2) transcriptome, (3) entire genome.
We are not aware of software that is capable of doing this directly, but with `samtools view` you must be able to get pretty close.

### Trim adapters

Small RNA-Seq protocols often (always?) provide reads that are contaminated with so called adapter sequences.
These sequences are manually added to your small ncRNAs and have to be removed in order to align properly.
These adapter sequences are often specific per protocol/sequencer, so please ensure you have access to the correct adapters.
There are several tools that can remove them.
Here a list of some of the tools that can be used to remove adapters:

*	CLC Bio (<FONT COLOR='red'>commercial</FONT>):	[http://www.clcbio.com/](http://www.clcbio.com/)
*	Cutadapt:	[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/), [https://github.com/marcelm/cutadapt/](https://github.com/marcelm/cutadapt/)

We are aware of samples that also have fixed-size random primers, besides their adapters (with a specific sequence). Be aware that these also have to be trimmed. This usually requires multiple trimming steps.

It is also convenient to remove trimmed reads that are really small (e.g. <= 15bp).

#### Quality trimming

Although it is common in RNA-Seq to trim low quality bases from your reads, we strongly recommend NOT to do this prior to running FlaiMapper-3. This will introduce a bias; fragments will most likely become shortened. If you really doubt the quality of a read, we advise you to neglect/discard the read.

### Choose aligner

The main complexity in RNA-Seq is splicing. There are several widely used free alignment programs for RNA-Seq. The same principle holds for small RNA-Seq. However, we are (at the moment) not aware of splicing events in ncRNAs other than tRNAs. The splice junctions in tRNAs are small (~10bp) but might be too long to be aligned using classicial alignment without adjusted settings. Therefore, if you want to align reads to pre-tRNAs, you want your aligner to understand splicing. Otherwise, if you want to use a 'classical' (non-splicing-aware) aligner, you want your introns to be removed prior to alignment. If you are not focussing on tRNAs at all, you probably also do not need your aligner to be aware of splicing and you can continue with classical aligners.

#### Full genome alignment

In case you want to align to the entire genome, you are not making use of mature ncRNAs. Also, if you want to include results from immature ncRNAs that undergo splicing (e.g. pre-tRNAs) you should proceed with a splicing-aware aligner. The following aligners are popular and should have sufficient documentation in order to run an alignment:

*	TopHat: [http://ccb.jhu.edu/software/tophat/index.shtml](http://ccb.jhu.edu/software/tophat/index.shtml)
*	RNA STAR: [http://code.google.com/p/rna-star/](http://code.google.com/p/rna-star/)
*	GSNAP: [http://research-pub.gene.com/gmap/](http://research-pub.gene.com/gmap/)
*	SubRead: [http://subread.sourceforge.net/](http://subread.sourceforge.net/)
*	HISAT: [http://ccb.jhu.edu/software/hisat/index.shtml](http://ccb.jhu.edu/software/hisat/index.shtml)

#### ncRNAdb09 alignment

If your reference consists of mature ncRNAs or you are sure you do not take results of ncRNAs that undergo alternative splicing (tRNAs) into account, you can use the following 'classical' aligners:

*	BWA: [http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)
*	bowtie: [http://bowtie-bio.sourceforge.net/bowtie2/index.shtml](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
*	SubRead: [http://subread.sourceforge.net/](http://subread.sourceforge.net/)
*	NovoAlign (<FONT COLOR='red'>commercial</FONT>): [http://www.novocraft.com/](http://www.novocraft.com/)
*	CLC Bio (small RNA-Seq module; <FONT COLOR='red'>commercial</FONT>):	 [http://www.clcbio.com/](http://www.clcbio.com/)
	*	Although we have used CLC for our analysis, we **do not** recommend using it prior to FlaiMapper-3. Exporting to SAM/BAM file aggregates all reads with an identical sequence and exporting the tables does not provide the coordinates of the aligned reads. The SAM/BAM aggregation affects the peak-detection of FlaiMapper-3. We have solved this issue by doing a first alignment round in CLC, to link the reads to their corresponding pre-cursor ncRNAs. We then apply a second alignment using MUSCLE (http://nar.oxfordjournals.org/content/32/5/1792.long), wrapped by a program called SSLM (http://www.gatcplatform.nl/), to find the exact coordinates of the reads linked to their precursor. We wrote a program to converts the SSLM format into BAM to ensure compatibility with other tools. To convert MUSCLE's output as wrapped by SSLM into BAM proceed with the following command(s):

		sslm2sam -o alignment.sam sslm_directory
		samtools view -h -bS alignment.sam > alignment.unsorted.bam
		
		# Depending on version of samtools - 1.+:
		samtools sort -o alignment.sorted.bam alignment.unsorted.bam
		
		# 0.18/0.19:
		samtools sort alignment.unsorted.bam > alignment.sorted.bam
		
		# Cleanup
		rm -r osslm_directory; rm alignment.unsorted.bam ; rm alignment.sam

Most aligners require tool-specific reference files (indexed versions of the reference genome). You can find the ncRNAdb09 reference for the following aligners:

| *Aligner* | *Reference file(s)* |
|:----------|:--------------------|
| TopHat & Bowtie | [ncrnadb09.ebwt.tar.gz](https://github.com/yhoogstrate/flaimapper/raw/master/share/annotations/ncRNA_annotation/ncrnadb09.ebwt.tar.gz) |
| TopHat2 & Bowtie2 | [ncrnadb09.bt2.tar.gz](https://github.com/yhoogstrate/flaimapper/raw/master/share/annotations/ncRNA_annotation/ncrnadb09.bt2.tar.gz) |

### Read collapsing

In small RNA-Seq it is not uncommon to apply read collapsing. This technique simply merges all reads with an identical sequence into one sequence, which sharply reduces the number of aligned reads. The consequence is that the peak detection will be applied upon much lower numbers, lowering the resolution of the experiment. This will most likely affect your results in a negative way. Therefore we advice you **not to use read collapsing** before FlaiMapper-3, unless you have a clear reason to believe it will give a better answer to your biological question(s) or improve your outcome.

### Multi-mapping

Some small ncRNAs have multiple genomic copies or share identical regions of their sequence with others. Reads that align to such a location(s) are called multi-map reads, since they have multiple candidate genomic origins. There are several strategies to deal with a multi-map read. Imagine we detected a read 6 times and it aligns 100% correctly to 3 different ncRNAs (multi-map regions). There are several strategies to deal with this situation:

1.	**All reads are assigned to only one of the multi-map locations.** In this case you would have 1 location with 6 reads, and two locations with 0 reads. Disadvantage: you are (most likely) under-representing the other locations.
2.	**All reads are mapped to all multi-map locations.** In this case you would end up with 6 aligned reads mapped to each of the 3 locations. This means that the total number of aligned reads has been multiplied with the number of multi-maps (3*6=18). Disadvantage: you are most likely over-representing most/all of the locations.
3.	**All reads are proportionally or randomly distributed over all multi-map locations.** If you do this proportionally, you would end up with 2 aligned reads for each of the 3 locations. Disadvantage: for all locations you are (most likely) either over- or under- representing them since.

Each of these strategies are based upon different assumptions and have their own disadvantages. It is difficult to state which strategy is the best, also because this may be dependent on where you want to use the outcome for. In case you want to enlist all possible fragments in your sample(s), we propose to use strategy 2. A single multi-map read, as measured by the sequencer, can not be unambiguously indicate its true genomic origin. Therefore, we believe that by retaining as much as possible multi-map reads, flaimapper will probably report as many as possible scenario's and annotate all possibilities. Since the intention is to enlist all possibilities, we would rather find a fragment more often due to its homology than missing some because of underrepresentation. Surely, down-stream research can still validate a fragment annotation in a multi-map region if necessary.

### Alignment indexing

The pysam library requires (position) sorted bam-files with a corresponding index file. Although the described alignment methods should by default produce (position) sorted bam-files, you can accomplish this by running the following command (notice that the second parameter has **no** *.bam* suffix):

	samtools sort output_unsorted.bam output

Only if a bam-file is position sorted it can be indexed. This is done using the following command:

	samtools index output.bam

## Run FlaiMapper-3

After you have aligned your reads, you can proceed with FlaiFapper. Here we will explain which commands you need to type in your terminal to analyse your data using FlaiMapper-3.

### Usage

The usage of FlaiMapper-3 (using BAM formatted files as input) is as follows:

	usage: flaimapper [-h] [-V] [-v | -q] [-p PARAMETERS] [-o OUTPUT] [-f {1,2}]
	                  [-r FASTA] [--offset5p OFFSET5P] [--offset3p OFFSET3P]
	                  alignment_file
	
	positional arguments:
	  alignment_file        indexed SAM or BAM file
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -V, --version         show program's version number and exit
	  -v, --verbose
	  -q, --quiet
	  -p PARAMETERS, --parameters PARAMETERS
	                        File containing the filtering parameters, using
	                        default if none is provided
	  -o OUTPUT, --output OUTPUT
	                        output filename; '-' for stdout
	  -f {1,2}, --format {1,2}
	                        file format of the output: [1: table; per fragment],
	                        [2: GTF (default)]
	  -r FASTA, --fasta FASTA
	                        Single reference FASTA file (+faid index) containing
	                        all genomic reference sequences
	  --offset5p OFFSET5P   Offset in bp added to the exon-type annotations in the
	                        GTF file. This offset is used in tools estimating the
	                        expression levels (default=4)
	  --offset3p OFFSET3P   Offset in bp added to the exon-type annotations in the
	                        GTF file. This offset is used in tools estimating the
	                        expression levels (default=4)

The '<CODE>\-\-verbose</CODE>' and '<CODE>\-\-quiet</CODE>' arguments change the level of verbosity. If '<CODE>\-\-verbose</CODE>' is enabled, FlaiMapper-3 will give more details about progress.

### Input: BAM

The FlaiMapper-3 binary that uses BAM formatted input files is called '*flaimapper*' and can be executed with the following terminal command:

For alignment **to reference genomes (e.g. hg19):**

	flaimapper \
	    -r hg19_full.fasta \
	    -o results_flaimapper.tabular.txt \
	    alignment_01.bam

For alignment **to ncRNAdb09:**

	flaimapper \
	    -r ncrnadb09.fasta \
	    -o results_flaimapper.tabular.txt \
	    alignment_02.bam

Remark that the backslashes are used to make the command continue at the next line and they can be removed when the command is written on a single line.

#### The "<CODE>\-\-parameters</CODE>"-argument

The filter function uses a set of parameters, which after installation can be found by running the following python code:

	python
	
	>>> from pkg_resources import resource_filename
	>>> print resource_filename("flaimapper","data/parameters.default.txt")

This is a tabular file with in the first column the position relative to the peak and in the second column the percentage to duck the other peaks intensity.
So, if at `-5` bases from your peak is another peak with an intensity of 200, and it will be filtered with 24.9%, the value after ducking will be 49,8.

#### The "<CODE>\-\-fasta</CODE>"-argument

The tabular output formatsis able to provide the fragments sequences of a FASTA file is provided.
If the file is not provided the program will not terminate, but will indicate wheterh the file was provided  in the header column and empty sequences shall be reported.
This file needs to be given separately because the BAM format does not store actual referenec sequences.
To ensure FlaiMapper-3 has access to the reference sequence(s), you can provide the reference genome as a single indexed FASTA file using the '<CODE>\-\-fasta</CODE>' argument as follows:

	flaimapper -r ncrnadb09.fa [...]

	flaimapper -r hg19_full.fa [...]

For ncRNAdb09, the FASTA file (and corresponding index) are available at the following url:

[ncrnadb09.fa](https://raw.github.com/yhoogstrate/flaimapper/master/share/annotations/ncRNA_annotation/ncrnadb09.fa)

[ncrnadb09.fa.fai](https://raw.github.com/yhoogstrate/flaimapper/master/share/annotations/ncRNA_annotation/ncrnadb09.fa.fai)

You need to have a <U>corresponding index file</U> under the name '*<prefix>.fa.fai*' besides your FASTA reference file. From FlaiMapper-3 and onwards this file will be automatically generated.

*Yet we do not provide any other reference genome than the Human ncRNAdb09.*

### Input: SSLM
The SSLM input has been terminated. There is, however, a converter to SAM available. To convert any SSLM directory to BAM, proceed with:

	sslm2sam -o alignment.sam sslm_directory
	samtools view -h -bS alignment.sam > alignment.unsorted.bam
	
	# Depending on version of samtools - 1.+:
	samtools sort -o alignment.sorted.bam alignment.unsorted.bam
	
	# 0.18/0.19:
	samtools sort alignment.unsorted.bam > alignment.sorted.bam
	
	# Cleanup
	rm -r osslm_directory; rm alignment.unsorted.bam ; rm alignment.sam

### Input: multiple alignments

FlaiMapper-3 is not able to deal with multiple input files anymore. Use `samtools merge prior` to running FlaiMapper-3 instead.

### Output: formats

FlaiMapper-3 can export results into the following formats:

- Tabular #1, per fragment
  * The start- and end-positions are 0-based.
- GTF/GFF
  * Every entry contains 2 lines: one of type *sncdRNA* (small non-coding derived RNA), which is the actual prediction and one of type *exon*, used by default by expression estimation tools.

The output format can be chosen with the '<CODE>\-f</CODE>' or the '<CODE>\-\-format</CODE>' argument, where the following argument have the following meaning:

## Reproduce article data

The raw figures used for the publication could be (re-)generated by running the scripts in the '*[scripts](https://github.com/yhoogstrate/flaimapper/tree/master/scripts/)*' directory.
The master script '*[scripts/analysis.sh](https://github.com/yhoogstrate/flaimapper/blob/master/scripts/analysis.sh)*' should run all analysis sequentially and make the directory structure within FlaiMapper-3's '*[output](https://github.com/yhoogstrate/flaimapper/tree/master/output)*' directory.

These files were generated using FlaiMapper-1. FlaiMapper-3 simplifies the interface but is backwards incompatible. You can still checkout older revisions of FlaiMapper in the git archive.

## Authors & Citing

The people who have contributed to this project are:

 - Youri Hoogstrate
 - Guido Jenster
 - Elena S. Martens-Uzunova

You can find the corresponding article at:

http://dx.doi.org/10.1093/bioinformatics/btu696

Youri Hoogstrate, Guido Jenster, and Elena S. Martens-Uzunova
FlaiMapper: computational annotation of small ncRNA derived fragments using RNA-seq high throughput data
Bioinformatics first published online October 22, 2014 doi:10.1093/bioinformatics/btu696
