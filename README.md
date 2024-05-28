![ProActiveV2](https://github.com/jlmaier12/ProActive/assets/45083046/e298f3c0-c8b2-477f-a37c-993b90bdef6e)

Detect elevations and gaps in read coverage on metagenome contigs or on assembled genomes. Regions of elevated read coverage may be associated with mobile genetic elements (MGE) that are actively replicating and/or are highly abundant, like prophages and transposases and translocases. Gaps in read coverage may be associated with genetic elements that are not homogenously integrated/present in the bacterial 'host' population. ProActive is best used as a screening method to identify genetic regions for further investigation.

### Decription
`ProActive` detects elevations and gaps in read coverage using read coverage pattern-matching. Elevations or gaps in read coverage can be caused by differential abundance of specific genetic elements. For example, an elevation in read coverage may be caused by prophage activation. When a prophage activates and enters the lytic cycle, its genome begins replicating and the ratio of phage:bacterial genomes in the cell begins to increase. Because there are more phage genomes than bacterial genomes, during sequencing more phage reads are generated than bacterial. When these reads are mapped back to their associated contig, the read coverage of the prophage region will be elevated in comparison to the read coverage of the bacterial genome on either side of the prophage. This same principle applies to temperate phage who are highly abundant in the environment as well as other mobile genetic elements that are freely present in the environment at a higher ratio than the originiating or 'host' genome. 

Figure 1 from Kieft and Anantharaman (2022) is an excellent visualization of this phenomenon with active prophage:

![msystems 00084-22-f001](https://github.com/jlmaier12/ProActive/assets/45083046/7f1d4e54-8ae8-406e-940e-5da311718dba)

**Reference** Kieft K, Anantharaman K. Deciphering Active Prophages from Metagenomes. mSystems. 2022 Apr 26;7(2):e0008422. doi: 10.1128/msystems.00084-22. Epub 2022 Mar 24. PMID: 35323045; PMCID: PMC9040807.

Conversely, a gap in read coverage may indicate genetic heterogeniety in the associated bacterial population. Genetic varients with and without specific genetic elements, like prophage or certain genes, will produce differential abundances of sequencing reads that may form read coverage gaps. The formation of read coverage gaps due to genetic varients is dependent on the assembly (i.e. if the assembler assembles the genetic varients as seperate entities or not) 

The increase or decrease in read coverage at MGE locations can be assessed and inferences can be made about associated activity/abundance, however, in order to assess elevation or gaps in read coverage, the genomic coordinates of the MGE(s) must be known. Locations of these genetic elements are usually idenfitied, at least in part, using gene annotations which means unannotated MGEs are missed. ProActive, however, is reference-independent and therefore bypasses this barrier in identification.

ProActive works by detecting the 'block-like' pattern generated by elevations or gaps in read coverages associated with integrated genetic elements. ProActive scans the mapped read coverages of a genome or metagenome assmebled contigs and search for elevations or gaps in read coverage. ProActive generates a summary table containing the sizes of all the detected read coverage elevations and gaps, the associated start and stop positions, and the associated ratio of elevated:non-elevated read coverage. 

ProActive's read coverage pattern-matching is only as good as the provided data. There are other, non-biological, situations that can create gaps or elevations in read coverage. For example, chimeric assemblies or contaminants can also create odd read coverage artifacts. Manual curation of ProActive's results is required, however, ProActive has several metrics with which a user can filter results to ease the manual curation process. In addition, ProActive's pattern-matching search can be modified with several user-defined input parameters. See the tutorial below for details. 


### Data input for `ProActive`

Map your sequencing read to the associated genome or metagenome contigs. Any read mapper can be used, however we recommend using [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/) with `ambiguous=random`, `qtrim=lr`, and `minid=0.97`. 

Create pileup files using `BBMap` `pileup.sh` with 100 bp window sizes:
```{bash}
$ pileup.sh in=YOUR_SORTED_READ_MAPPING.bam out=output.pileupcovstats bincov=pileup_bincov100.txt binsize=100 stdev=t
```
The pileup_bincov100.txt file produced with `pileup.sh` can be used directly as input for ProActive.

NOTE: ProActive will filter out contigs shorter than 30kbp. If your metagenome assembly consists of majority short contigs, ProActive may not be the right tool to use. Filtering out contigs less than 30kbp prior to read mapping will result in smaller bam and pileup files and subsequent faster processing by ProActive. 
Hint: filtering contigs by length can be done easily with BBMap's `reformat.sh` function (`reformat.sh in=contigs.fasta out=filtered_contigs.fasta minlength=30000`)

### Install and run `ProActive`

##### Install `ProActive`:

```R
library(devtools)
install_github("jlmaier12/ProActive", ref="master")
library(ProActive)
```

##### Using `ProActive`:

Import your pileup file:
```R
my_pileup <- read.delim("pileup_bincov100.txt",  header=FALSE, comment.char="#")
```

Run `ProActive()`:
```R
ProActiveResults <- ProActive(pileup=my_pileup, mode="metagenome", windowsize=2000, minsize=10000, maxsize=80000, chunksize=NA, nones=FALSE, cleanup=TRUE) 
```
Input parameters:
- pileup: Your pileup file with 100bp binsizes (aka windowsize)
- mode: Either "genome" or "metagenome". 
- windowsize: Either 200, 500, 1000 and 2000. The number of base pairs to average coverage values over. Larger window sizes improve processing time but the resolution of read coverage patterns may be lost. 1000bp windowsize is the default.
- minsize: The minimum size (in base pairs) of elevated read coverage you would like ProActive to search for. Default is 10000. Value should be divisible by 100. 
- maxsize: The maximum size (in base pairs) of elevated read coverage you would like ProActive to search for. Default is NA (i.e no upper limit set). Value should be divisible by 100.
- chunksize: If mode="genome", ProActive will 'chunk' the genome to search for a region of elevated read coverage on each chunk. Since ProActive can only identify one region of elevated read coverage per 'chunk', choosing smaller chunks will increase sensitivity AND processing time. Default is 100,000bp per chunk, with the exception of the last chunk which is the remainder number of base pairs once the genome has been evenly 'chunked'. Chunk size should be divisible by 100.
- nones: TRUE or FALSE. Default is FALSE. TRUE if you would like to include the contigs classified as 'none' in your final pattern lists. This can be useful if you are comparing ProActive results between sample conditions and would like to overlay the read coverage pattern-matching results, including 'negative' results, for a specific contig or genome.
- cleanup: Either TRUE or FALSE. ProActive will clean and reformat your input pileup_bincov100.txt file for proper compatibility with functions of ProActive if cleanup=TRUE. Default is TRUE. 


Output:
The output of `ProActive()` is a list. Assign list items of interest to their own variables for further use:
```R
summarytable <- ProActiveResults$SummaryTable
VC_PatternMatches <- ProActiveResults$VeryConfidentPatternMatches
C_PatternMatches <- ProActiveResults$ConfidentPatternMatches
NC_PatternMatches <- ProActiveResults$NotConfidentPatternMatches
FilterOutSummaryTable <- ProActiveResults$FilteredOut
```
- ProActiveResults$SummaryTable: A dataframe containing the results of ProActive's pattern-matching on each contig or genome chunk. The table contains the confidence of each pattern match and the associated match size, start and stop positions, and elevation ratio. 
- ProActiveResults$VeryConfidentPatternMatches: Pattern matches where all the coverage values within the pattern match region are the largest coverage values on the contig or genome chunk. 
- ProActiveResults$ConfidentPatternMatches: Pattern matches where the maximum coverage value for the contig or genome chunk falls within the pattern match region.
- ProActiveResults$NotConfidentPatternMatches: Pattern matches where the maximum coverage value for the contig or genome chunk does NOT fall within the pattern match region.
- ProActiveResults$FilteredOut: A dataframe containing the contigs or genome chunks that were filtered out due to being shorter than 30kbp or having low read coverage. 

Plot mapped read coverage plots of contigs or genome chunks detected with regions of elevated read coverage. The 'block-like' pattern that ProActive used to make its prediction for each contig/genome chunk will be overlayed on the barplot in black:
```R
#plot confident predictions (C_PatternMatches) that have an associated elevated:not-elevated read coverage ratio (elevation_filter) of 1.5 or greater
ProActivePredictionPlots(pileup=my_pileup, ProActiveResults=ProActiveResults, ProActive_shapelist=C_PatternMatches, elevation_filter=1.5, cleanup=TRUE)
```
- pileup: Your pileup file with 100bp binsizes (aka windowsize)
- ProActiveResults: The output of `ProActive()`
- ProActive_shapelist: Either the ProActiveResults$VeryConfidentPatternMatches, ProActiveResults$ConfidentPatternMatches, or ProActiveResults$NotConfidentPatternMatches from the ProActive output. 
- elevation_filter: Only plot contigs or genome chunks with a pattern match greater than a specified mean elevated:non-elevated read coverage ratio value. Default is NA. 
- cleanup: Either TRUE or FALSE. ProActive will clean and reformat your input pileup_bincov100.txt file for proper compatibility with functions of ProActive if cleanup=TRUE. Default is TRUE. 


##### Search for gene annotations of interest:

Generating your gff file:
Your gff file may be generated with any annotation tool. If you are annotating metagenome contigs, it is a good idea to only annotate contigs greater than 30kbp as `ProActive()` filters out contigs less than 30kbp.  

gff file format:
Depending on how your gff file was generated, it may contain additional information above and below the main table containing gene annotation information. To import your gff file to R, it must have this additional inforamtion removed. The easiest way to do this is to select only the necessary information from your gff and copy it to a new file. This can be done on the command line using grep. Using nano, open your gff file to see what information it contains. Scroll down until you reach the tab-seperated table containing the gene annotation information. The first column should contain the reference names for your contigs or genome. The reference name should begin with a common identifier which can be used to extract the associated lines and move them to a new file with grep. 
```bash
$ grep ^COMMON_IDENTIFIER geneannots.gff > cleaned_geneannots.gff 
```

Import your gff file to R:
```R
gene_annotations <- read.delim("cleaned_geneannots.gff", header=FALSE)
```
Your gff file should have columns in the following order: contig or genome reference (i.e "NODE_#"), tool (i.e "prodigal"), CDS, start position, stop position, '.', '+/-', '0', annotations.

Clean and reformat your gff file with ProActive's built-in function. This only needs to be done once:
```R
clean_gene_annotations <- gff_cleanup(gene_annotations)
```

Search for gene annotations matching of interest:
```R
GeneAnnotationSearch(ProActiveResults=ProActiveResults, ProActive_shapelist=C_PatternMatches, pileup=my_pileup, gene_annots=clean_gene_annotations, geneorproduct="product", keywords=c("phage","tail","spike","needle"), genelocation="specific", bprange = 0, cleanup=TRUE, specificcontig=NA) 
```
Input parameters:
- ProActiveResults: The output of `ProActive()`
- ProActive_shapelist: Either the ProActiveResults$VeryConfidentPatternMatches, ProActiveResults$ConfidentPatternMatches, or ProActiveResults$NotConfidentPatternMatches from the ProActive output. 
- pileup: Your pileup file with 100bp binsizes (aka windowsize)
- gene_annots: A gff file that has been 'cleaned' with `gff_cleanup()`
- geneorproduct: Either "gene" or "product".
- keywords: The gene or product annotation key-word(s) you would like to search for. Case independent. Key-word(s) must be in quotes, comma-separated, and surrounded by c() i.e( c("antibiotic", "resistance", "drug") )
- genelocation: "specific" or "nonspecific". Search for gene-annotations that match provided keywords only within the pattern-match region with "specific" or search the entire contig or genome chunk for matching annotations with nonspecific". Default is "nonspecific".
- bprange: If you are searching for gene annotations with genelocation="specific", you may specify the region (in basepairs) that should be searched to the left and right of the pattern match in addition to the pattern match region itself. Default is 0 (if you choose genelocation="specific" and don't set the bprange, ProActive will only search for matching gene annotations within the pattern-match region)
- cleanup: Either TRUE or FALSE. ProActive will clean and reformat your input pileup_bincov100.txt file for proper compatibility with functions of ProActive if cleanup=TRUE. Default is TRUE. 
- specificcontig: Search for gene annotations matching the provided keywords only on a specific contig or genome chunk. Provide contig or genome chunk reference name in quotes.


Output:
`GeneAnnotationSearch()` will print read coverage plots for contigs or genome chunks that contain gene annotations that match the provided keyword(s). The locations of the matching gene annotations will be displated on the graph with black vertical lines. The start and stop position of the associated ProActive pattern-match are displayed with red vertical lines. The graphs are plotted without re-averaging the windowsizes in order to provide better resolution of gene annotation locations. 
