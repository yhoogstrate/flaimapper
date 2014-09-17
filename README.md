# Fragment Location Annotation Identification Mapper #
---
	 ___    _      _____  _         _____  ___    ___    ___    ___   
	(  _`\ ( )    (  _  )(_)/'\_/`\(  _  )(  _`\ (  _`\ (  _`\ |  _`\ 
	| (_(_)| |    | (_) || ||     || (_) || |_) )| |_) )| (_(_)| (_) )
	|  _)  | |  _ |  _  || || (_) ||  _  || ,__/'| ,__/'|  _)_ | ,  / 
	| |    | |_( )| | | || || | | || | | || |    | |    | (_( )| |\ \ 
	(_)    (____/'(_) (_)(_)(_) (_)(_) (_)(_)    (_)    (____/'(_) (_)

---

- [Download & Installation](#download--installation)
    - [Install dependencies](#install-dependencies)
    - [Download FlaiMapper](#download-flaimapper)
        - [Latest version from GitHub](#latest-version-from-github)
        - [Read & agree with GPL3 free software license](#read--agree-with-gpl3-free-software-license)
    - [Compile & Install](#compile--install)
    - [Check if it works!](#check-if-it-works)
    - [Uninstall FlaiMapper](#uninstall-flaimapper)
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
- [Run FlaiMapper](#run-flaimapper)
    - [Usage](#usage)
    - [Input: BAM](#input-bam)
         - [The "\-\-mask"-argument](#the---mask-argument)
         - [The "\-\-fasta"-argument](#the---fasta-argument)
    - [Input: SSLM](#input-sslm)
    - [Input: multiple alignments](#input-multiple-alignments)
    - [Output: formats](#output-formats)
- [Reproduce article data](#reproduce article data)
- [Authors & Citing](#authors--citing)

## Download & Installation
### Install dependencies

Make sure you have *python2*, *pip* and the *pysam* library installed. Pysam is not by default installed on many systems. You can find installation details on at the following urls:

*	[https://github.com/pysam-developers/pysam](https://github.com/pysam-developers/pysam)
*	[https://github.com/pysam-developers/pysam/master/INSTALL](https://github.com/pysam-developers/pysam/master/INSTALL)

The easiest way to install the latest version of pysam is via pip using the following terminal command:

	sudo apt-get install python-dev python-pip

	sudo pip install --upgrade pyrex

	sudo pip install --upgrade pysam


### Download FlaiMapper
#### Latest version from GitHub

To get the latest version of FlaiMapper please download with Git (make sure it is installed) from GitHub using the following terminal command:

	sudo apt-get install git

	git clone --recursive https://github.com/yhoogstrate/flaimapper.git

Or you can download the latest source as ZIP package from the following url: [here (right mouse click; save as)](https://github.com/yhoogstrate/flaimapper/archive/master.zip). Extract this package to a directory called "flaimapper".

#### Read & agree with GPL3 free software license

If you have downloaded FlaiMapper, you should read the GNU General Public License v3.0 to make sure you understand what you are allowed to do with the source files and what free software means. It also states that the software is distributed in the hope to be useful but without any warranty. The license can be accessed directly at the following url:

*	[https://github.com/yhoogstrate/flaimapper/raw/master/LICENSE](https://github.com/yhoogstrate/flaimapper/raw/master/LICENSE)

### Compile & Install
Browse into the src directory you just have created: 

	cd flaimapper/src

The installation procedure of FlaiMapper first converts the python code into byte code, followed by a mechanism that installs flaimapper in the system directories. You can achieve this by running the following two commands: 

	python setup.py build

	sudo python setup.py install

### Check if it works!

Last but not least, please check if FlaiMapper runs. You can do this by running either a full analysis or just run the following terminal command: 

	flaimapper --help

If FlaiMapper givers errors, warnings or does not install, please do not hazitate and either submit the bug or send a fix to the GitHub repository at the following url: <A HREF='https://github.com/yhoogstrate/flaimapper' TARGET='_new'>https://github.com/yhoogstrate/flaimapper</A>

Otherwise: congratulations, you have just installed FlaiMapper!

### Uninstall FlaiMapper

To remove FlaiMapper automatically from your system directories (without removing your data files), proceed with the following terminal command:

	sudo pip uninstall flaimapper

The downloaded source files can only be removed manually. Be aware that references files for *ncRNAdb09* and the data of the analysis as demonstrated in the corresponding article are also located within this directory.

## Alignment

### Obtain reference

We propose multiple reference strategies for alignment, but we let it up to the researcher to decide which strategy is preferred.

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

Some aligners have a related feature implemented: an indexed transcriptome besides the indexed reference genome. You might think of tophat's "<CODE>--transcriptome-index</CODE>" and RNA-STAR's "<CODE>sjdbGTFfile</CODE>". However, idealy you want to have a three phase alignment: (1) ncRNAs, (2) transcriptome, (3) entire genome. We are not aware of software that is capable of doing this (yet?).

### Trim adapters

Small RNA-Seq protocols often provide reads that are contaminated with so called adapter sequences.
These sequences are manually added to your small ncRNAs and have to be removed in order to align properly.
These adapter sequences are often specific per protocol/sequencer, so please ensure you have access to the correct adapters.
There are several tools that can remove them.
Here a list of some of the tools that can be used to remove adapters:

*	CLC Bio (<FONT COLOR='red'>commercial</FONT>):	[http://www.clcbio.com/](http://www.clcbio.com/)
*	FASTX toolkit:	[http://hannonlab.cshl.edu/fastx_toolkit/](http://hannonlab.cshl.edu/fastx_toolkit/)
*	Scythe: [https://github.com/vsbuffalo/scythe](https://github.com/vsbuffalo/scythe)

We are aware of samples that also have fixed-size random primers, besides their adapters (with a specific sequence). Be aware that these also have to be trimmed. This usually requires multiple trimming steps.

#### Quality trimming

Although it is common in RNA-Seq to trim low quality bases from your reads, we strongly recommend you NOT to do this prior to FlaiMapper. This will introduce a bias; fragments will most likely become shortened. If you really doubt the quality of a read, we advise you to neglect and discard the entire read!

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
	*	Although we have used CLC for our analysis, we **do not** recommend using it prior to FlaiMapper. Exporting to SAM/BAM file aggregates all reads with an identical sequence and exporting the tables does not provide the coordinates of the aligned reads. The SAM/BAM aggregation affects the peak-detection of FlaiMapper. We have solved this issue by doing a first alignment round in CLC, to link the reads to their corresponding pre-cursor ncRNAs. We then apply a second alignment using MUSCLE (http://nar.oxfordjournals.org/content/32/5/1792.long), wrapped by a program called SSLM (http://www.gatcplatform.nl/), to find the exact coordinates of the reads linked to their precursor. We wrote a program to converts the SSLM format into BAM to ensure compatibility with other tools. To convert MUSCLE's output as wrapped by SSLM into BAM proceed with the following command(s):

			sslm2sam "output_sslm" -m ncrnadb09.gtf -o output_unsorted.sam
		
			samtools view -h -bS output_unsorted.sam > output_unsorted.bam
		
			samtools sort output_unsorted.bam output
		
			rm -r output_sslm; rm output_unsorted.bam

Most aligners require tool-specific reference files (indexed versions of the reference genome). You can find the ncRNAdb09 reference for the following aligners:

| *Aligner* | *Reference file(s)* |
|:----------|:--------------------|
| TopHat & Bowtie | [ncrnadb09.ebwt.tar.gz](https://github.com/yhoogstrate/flaimapper/raw/master/share/annotations/ncRNA_annotation/ncrnadb09.ebwt.tar.gz) |
| TopHat2 & Bowtie2 | [ncrnadb09.bt2.tar.gz](https://github.com/yhoogstrate/flaimapper/raw/master/share/annotations/ncRNA_annotation/ncrnadb09.bt2.tar.gz) |

If you think installing, configuring and creating references takes too long, you can make use of Galaxy. We have a public instance with ncRNAdb09 installed available at the following url:

*	[http://galaxy-sandbox.trait-ctmm.cloudlet.sara.nl/](http://galaxy-sandbox.trait-ctmm.cloudlet.sara.nl/)

### Read collapsing

In small RNA-Seq it is not uncommon to apply read collapsing. This technique simply merges all reads with an identical sequence into one sequence, which sharply reduces the number of aligned reads. The consequence is that the peak detection will be applied upon much lower numbers, lowering the resolution of the experiment. This will most likely affect your results in a negative way. Therefore we advice you **not to use read collapsing** before FlaiMapper unless you have a clear reason to believe it will give a better answer to your outcome or biological question(s).

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

## Run FlaiMapper

After you have aligned your reads, you can proceed with FlaiFapper. Here we will explain which commands you need to type in your terminal to analyse your data using FlaiMapper.

### Usage

The usage of FlaiMapper (using BAM formatted files as input) is as follows:

	usage: flaimapper [-h] [-V] [-v | -q] [-o OUTPUT] [-f FORMAT] -m MASK
	                  [-r FASTA]
	                  alignment_files [alignment_files ...]
	
	positional arguments:
	  alignment_files       indexed SAM or BAM files compatible with pysam
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -V, --version         show program's version number and exit
	  -v, --verbose
	  -q, --quiet
	  -o OUTPUT, --output OUTPUT
	                        output filename; '-' for stdout
	  -f FORMAT, --format FORMAT
	                        file format of the output: [1: table; per fragment],
	                        [2: table; per ncRNA], [3: genbank]
	  -m MASK, --mask MASK  GTF/GFF3 mask file (precursors)
	  -r FASTA, --fasta FASTA
	                        Single reference FASTA file (+faid index) containing
	                        all genomic reference sequences

The usage of FlaiMapper (using SSLM formatted data as input) is as follows:


positional arguments:
  alignment_directories
                        SSLM formatted output directories

	usage: flaimapper-sslm [-h] [-V] [-v | -q] [-o OUTPUT] [-f FORMAT]
	                       alignment_directories [alignment_directories ...]
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -V, --version         show program's version number and exit
	  -v, --verbose
	  -q, --quiet
	  -o OUTPUT, --output OUTPUT
	                        output filename; '-' for stdout
	  -f FORMAT, --format FORMAT
	                        file format of the output: [1: table; per fragment],
	                        [2: table; per ncRNA], [3: genbank]
	  -m MASK, --mask MASK  GTF/GFF3 mask file (precursors)
	  -r FASTA, --fasta FASTA
	                        Single reference FASTA file (+faid index) containing
	                        all genomic reference sequences

From this follows that you can find the version of your installed flaimapper with the following commands:

	flaimapper --version

	flaimapper-sslm --version

The '<CODE>\-\-verbose</CODE>' and '<CODE>\-\-quiet</CODE>' arguments change the level of verbosity. If '<CODE>\-\-verbose</CODE>' is enabled, FlaiMapper will give more details about progress.

### Input: BAM

The FlaiMapper binary that uses BAM formatted input files is called '*flaimapper*' and can be executed with the following terminal command:

For alignment **to reference genomes (e.g. hg19):**

	flaimapper \
	    -m ncrnadb09_hg19.gtf \
	    -r hg19_full.fasta \
	    -o results_flaimapper.tabular.txt \
	    alignment_01.bam

For alignment **to ncRNAdb09:**

	flaimapper \
	    -m ncrnadb09.gtf \
	    -r ncrnadb09.fasta \
	    -o results_flaimapper.tabular.txt \
	    alignment_02.bam

Remark that the backslashes are used to make the command continue at the next line and they can be removed when the command is written on a single line.

#### The "<CODE>\-\-mask</CODE>"-argument

The BAM (and SAM) alignment formats are tabular file formats that store a reads absolute start position, and a formal description of how the alignment proceeds. This also includes the reads sequence and can be provided with a quality score.
The BAM/SAM format describes the alignment of a read in relation to a so called *Reference sequence*.
It is very important to understand this, because the reference sequences may have a different meaning in the two proposed types of experiments.

 - In case you <U>align to a reference genome</U> (e.g. hg19), each *Reference sequence* represents a chromsome (or a contig). Within this chromosome, multiple ncRNAs can be located.
 - In case you <U>align to ncRNAdb09</U>, each *Reference sequence* represents exactly one mature ncRNA.

The consequence of this is that FlaiMapper must know where in which Reference sequence the ncRNAs are located. In FlaiMapper these so called MASK regions are provided as GTF/GFF files with the '<CODE>\-m</CODE>' or '<CODE>\-\-mask</CODE>' argument. For alignment to reference genome hg19, you have to provide the following argument:

	flaimapper -m ncrnadb09_hg19.gtf [...]

In ncRNAdb09, each Reference sequence represent exactly single ncRNA.
Therefore, the provided MASK for ncRNAdb09 describes only one ncRNA per (entire) reference sequence.

For alignment to ncRNAdb09, you have to provide the following argument:

	flaimapper -m ncrnadb09.gtf [...]

Currently we serve the ncRNAdb09 MASK as a GTF/GFF file for the following reference genomes:

| **Ref. Genome** | **Ref. Genome ID** | **GTF file** | **GTF index** |
|:----------------|:--------------------|:--------------|:----------------|
| ncRNAdb09 (Human) | ncrnadb09 | [ncrnadb09.gtf](https://github.com/yhoogstrate/flaimapper/raw/master/share/annotations/ncRNA_annotation/ncrnadb09.gtf) | [ncrnadb09.gtf.tbi](https://github.com/yhoogstrate/flaimapper/raw/master/share/annotations/ncRNA_annotation/ncrnadb09.gtf.tbi) |
| Human Feb. 2009 \(GRCh37/hg19\) | hg19 | [ncrnadb09_hg19.gtf](https://github.com/yhoogstrate/flaimapper/raw/master/share/annotations/ncRNA_annotation/ncrnadb09_hg19.gtf) | [ncrnadb09_hg19.gtf.tbi](https://github.com/yhoogstrate/flaimapper/raw/master/share/annotations/ncRNA_annotation/ncrnadb09_hg19.gtf.tbi) |

#### The "<CODE>\-\-fasta</CODE>"-argument

In contrast to formats that only contrain genomic coordines, like BED and GTF, the tabular output formats and GenBank also provide the fragments sequences.
It is important to understand is that within the BAM/SAM format no sequences of the reference genome are stored. Therefore it is not possible (feasible) to extract the sequence of a fragment from a BAM file.
To ensure FlaiMapper has access to the reference sequence(s), you can provide the reference genome as a single indexed FASTA file using the '<CODE>\-\-fasta</CODE>' argument as follows:

	flaimapper -r ncrnadb09.fa [...]

	flaimapper -r hg19_full.fa [...]

For ncRNAdb09, the FASTA file (and corresponding index) are available at the following url:

[ncRNAdb09\_with\_tRNAs\_and\_Pseudogenes\_\_21\_oct\_2011\_\_hg19.fasta](https://raw.github.com/yhoogstrate/flaimapper/master/share/annotations/ncRNA\_annotation/ncRNA\_annotation/ncRNAdb09\_with\_tRNAs\_and\_Pseudogenes\_\_21\_oct\_2011\_\_hg19.fasta)

[ncRNAdb09\_with\_tRNAs\_and\_Pseudogenes\_\_21\_oct\_2011\_\_hg19.fasta.fai](https://raw.github.com/yhoogstrate/flaimapper/master/share/annotations/ncRNA\_annotation/ncRNA\_annotation/ncRNAdb09\_with\_tRNAs\_and\_Pseudogenes\_\_21\_oct\_2011\_\_hg19.fasta.fai)

You need to have a <U>corresponding index file</U> under the name '*<prefix>.fa.fai*' besides your FASTA reference file.

*Yet we do not provide any other reference genome than the Human ncRNAdb09.*

### Input: SSLM

**<FONT COLOR="red">The SSLM format is deprecated. If possible, proceed with BAM/SAM files.</FONT>**

The SSLM format is the output format of Short Sequence Location Mapper:

[http://www.gatcplatform.nl/SSLM/index.html](http://www.gatcplatform.nl/SSLM/index.html)

Earlier analysis made use of the ncRNA analysis program SSLM which wraps the aligner MUSCLE ([http://dx.doi.org/10.1093/nar/gkh340](http://dx.doi.org/10.1093/nar/gkh340)). We initially proceeded with SSLM's corresponding data format, because it provided all neccesairy information (including reference sequence, in contrast to BAM). Furher analysis made FlaiMapper evolve into a tool, still specific for SSLM data. Because the BAM format has become the de facto standards for alignment data we have changed the input format to BAM. Because we still want to keep our software backwards compatible with the initial analysis (in order to reproduce the results), we also want to keep the deprecated SSLM format supported. The directory structure of an SSLM experiment is as follows:

	.
	├── idreadable.txt
	└── validated
	    ├── file100.fa
	    ├── file101.fa
	    ├── file102.fa
	    ├── file103.fa
	    ├── file104.fa
	    ├── file105.fa
	    └── file106.fa

The _idreable.txt_ is a tab-delimted file and links the name of a ncRNA to its alignment (as a 'validated/*.fa' file). The _idreadable.txt_ has a syntax comparable to:

	sequence  »  filename
	>NAME=TRNAValAAC&LOCI=[chr1:180184276-180184348:strand=-]&SOURCE=UCSC&SOURCE-ACCESSION=chr1.tRNA63-ValAAC&GENOME=hg19  »  file100
	>HGNC=38236&HUGO-Symbol=MIR3198-1&HUGO-Name=microRNA_3198-1&LOCI=[chr22:18246936-18247035:strand=-]&SOURCE=RefSeq&SOURCE-ACCESSION=NR_036168&GENOME=hg19  »  file101

An individual alignment file in the FASTA (_*.fa_) format corresponds to an alignment of one ncRNA sequence, is located in a directory valled '**validated**' and has a sytax comparable to:

	>NAME=TRNAVALAAC&LOCI=[CHR1:180184276-180184348:STRAND=-]&SOURCE=UCSC&SOURCE-ACCESSION=CHR1.TRNA63-VALAAC&GENOME=HG19
	GTTTCCATAGTGTACTGGTTATCACATTCACCTAACACGCGAAAGGTCCTTGGTTTGAAACCAGGCAGAAACACCA
	>ReadID#1_GTTTCCATAGTGTAGTGGTTATC_hits8
	GTTTCCATAGTGTAGTGGTTATC-----------------------------------------------------
	>ReadID#2_GTTTCCATAGTGTAGTGG_HITS7
	GTTTCCATAGTGTAGTGG----------------------------------------------------------
	>ReadID#3_GTTTCCATAGTGTAGTGGTTAT_x6
	GTTTCCATAGTGTAGTGGTTAT------------------------------------------------------

The '<CODE>_hits...</CODE>' and '<CODE>_x...</CODE>' suffixes are indicators for the number of indentical copies of the read.

To run flaimapper using SSLM data, we provide the paths of the directories that contain the SSLM data.
In the article experiment '*[SRP006788](https://github.com/yhoogstrate/flaimapper/tree/master/share/small_RNA-seq_alignments/SRP006788)*' was analysed.
To analyse 1 of its 6 SSLM experiments '*[SRR207111_HeLa18-30](https://github.com/yhoogstrate/flaimapper/tree/master/share/small_RNA-seq_alignments/SRP006788/SRR207111_HeLa18-30)*', we run FlaiMapper (SSLM) as follows:

	flaimapper-sslm \
	    -m share/annotations/ncRNA_annotation/ncrnadb09.gtf
	    -r share/annotations/ncRNA_annotation/ncrnadb09.fa
	    -o output/FlaiMapper/SRP006788/01.a_output_flaimapper.txt \
	    share/small_RNA-seq_alignments/SRP006788/SRR207111_HeLa18-30

FlaiMapper will export the results as table into '*[output/FlaiMapper/SRP006788/01.a_output_flaimapper.txt](https://github.com/yhoogstrate/flaimapper/raw/master/output/FlaiMapper/SRP006788/01.a_output_flaimapper.txt)*'.

*The SSLM version of FlaiMapper does not support entire reference genome alignment.*

### Input: multiple alignments

FlaiMapper is able to deal with multiple input files. In certain situations you want to enhance your resulotion by combining datasets, multiple runs for example. You can simply enhance your resolution by stacking the alignment. Instead of creating merged BAM or SSLM files, you can simply tell FlaiMapper to use multiple alignments and FlaiMapper will simply reads through these alignments as if they were one alignment. So, if you provide multiple input files you will get only one output file based on the concatenated data. If you want to get **individual outputs for any of your samples**, you have to **run FlaiMapper separately on each sample!**

The last argument of FlaiMapper is simply a 1-or-multiple. You can run FlaiMapper on multiple files by separating all desired files with a space:

	flaimapper \
	    -f 1 \
	    -o output/FlaiMapper/SRP002175/01_output_flaimapper.txt \
	    -m share/annotations/ncRNA_annotation/ncrnadb09.gtf \
	    -r share/annotations/ncRNA_annotation/ncrnadb09.fa \
	        share/small_RNA-seq_alignments/SRP002175/SRR038852.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038853.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038854.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038855.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038856.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038857.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038858.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038859.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038860.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038861.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038862.bam \
	        share/small_RNA-seq_alignments/SRP002175/SRR038863.bam

	flaimapper-sslm \
	    -f 1 \
	    -m share/annotations/ncRNA_annotation/ncrnadb09.gtf \
	    -r share/annotations/ncRNA_annotation/ncrnadb09.fa \
	    -o output/FlaiMapper/SRP002175/01_output_flaimapper.txt \
	        share/small_RNA-seq_alignments/SRP002175/SRR038852 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038853 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038854 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038855 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038856 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038857 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038858 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038859 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038860 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038861 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038862 \
	        share/small_RNA-seq_alignments/SRP002175/SRR038863

Remark that the backslashes are used to make the command continue at the next line and they can be removed when the command is written on a single line.

### Output: formats

FlaiMapper can export results into the following formats:

- Tabular #1, per fragment
- Tabular #2, per ncRNA
- GenBank

The implementation of the following formats is under development:

- BED
- GTF/GFF

The output format can be chosen with the '<CODE>\-f</CODE>' or the '<CODE>\-\-format</CODE>' argument, where the following argument have the following meaning:

- <CODE>\-f 1</CODE>&nbsp; &nbsp; &nbsp; &nbsp; Tabular #1: <CODE>Fragment&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Precursor&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Fragment-start&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Fragment-stop&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Sequence&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Corresponding-reads</CODE>
- <CODE>\-f 2</CODE>&nbsp; &nbsp; &nbsp; &nbsp; Tabular #2: <CODE>Precursor&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Curated&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Fragment-1-start&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Fragment-1-stop&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Fragment-1-sequence&nbsp; <FONT COLOR="gray">&#187;</FONT>&nbsp; Fragment-2-...</CODE>
- <CODE>\-f 3</CODE>&nbsp; &nbsp; &nbsp; &nbsp; GenBank

The location of the output is defined with the '<CODE>\-o</CODE>' or '<CODE>\-\-output</CODE>' argument. If the argument is left empty or equal to '<CODE>\-</CODE>', FlaiMapper will write directly to stdout.

## Reproduce article data

All documents used for the publication can be (re-)generated by running the scripts in the '*[scripts](https://github.com/yhoogstrate/flaimapper/tree/master/scripts/)*' directory.
The master script '*[scripts/analysis.sh](https://github.com/yhoogstrate/flaimapper/blob/master/scripts/analysis.sh)*' should run all analysis sequentially and make the directory structure within FlaiMapper's '*[output](https://github.com/yhoogstrate/flaimapper/tree/master/output)*' directory.

## Authors & Citing

The people who have contributed to this project are:

 - Youri Hoogstrate
 - Elena S. Martens-Uzunova
 - Guido Jenster
