# ProActive
Active prophage detection in metagenome sequencing read coverages.

### Decription
ProActive detects regions of elevated read coverage representative of active propgages in metagenomic sequencing reads mapped to contigs. When a prophage activates and enters the lytic cycle, its genome begins replicating and the ratio of phage:bacterial genomes in the cell begins to increase. Because there are more phage genomes than bacterial genomes, during sequencing more phage reads are generated than bacterial. When these reads are mapped back to their associated contig, the read coverage of the prophage region will be elevated in comparison to the read coverage of the bacterial genome on either side of the prophage. 

Figure 1 from Kieft and Anantharaman (2022) is an excellent visualization of this phenomenon:

![msystems 00084-22-f001](https://github.com/jlmaier12/ProActive/assets/45083046/7f1d4e54-8ae8-406e-940e-5da311718dba)

**Reference** Kieft K, Anantharaman K. Deciphering Active Prophages from Metagenomes. mSystems. 2022 Apr 26;7(2):e0008422. doi: 10.1128/msystems.00084-22. Epub 2022 Mar 24. PMID: 35323045; PMCID: PMC9040807.

The increase in read coverage at prophage locations can be detected and inferences can be made about the lytic activity of the phage. However, in order to assess difference in read coverage between a prophage and the neighboring bacterial genome, the genomic coordinates of the prophage must be known. Prophage locations are usually idenfitied, at least in part, using phage-related gene annotations. Relying on gene annotations for active prophage detection means that only annotated phage are detected while the vast majority are missed. ProActive, however, is reference-independent and therefore bypasses this barier to active prophage identification.

ProActive works by detecting the 'shape' created by plotting the mapped read coverages of an active prophage on its associated contig. When the read coverages for an active prophage are plotted as a bar graph, the elevated read coverage at the prophage location creates a 'block-like' shape (visualized in the figure above). ProActive will scan through all the contigs in a metagenome and search each one for a 'block-like' shape with a base pair size that is realistic for that of a phage genome. ProActive will output a summary table containing all the conitgs with a potentially active prophage, the sizes of the active prophages, and read coverage graphs visualizing the active prophage with the shape that ProActive used to make its identification. As an optional feature, users may input a PROKKA gff file containing gene annotations for their contigs which ProActive will use to determine which of the identified active prophages have phage-related gene annotations. 

To be clear, while ProActive labels the elevated read coverage it detcts as active prophage, there are other situations that can create read coverage patterns similar to that of active prophage. For example, chimeric assemblies can create contigs with regions of elevated read coverage that are not associated with active prophage, but ProActive will likely identify this type of event as an active prophage. It is up to the user to assess the active prophage identifications and check for false positives. Typically, there tend to be more false positive identifications made in short, low-quality contigs. 
Additionally, ProActive will not identify dormant prophages or prophages that are only active in a small subset of a bacterial population as the elevated read coverage pattern will likely be hidden by the reads of bacteria without active prophage. 


### Data input for ProActive
ProActive requires one input file containinng read coverage summary information. We recommend the following workflow to get your metagenomic sequence data into the correct format for use with ProActive. 

On the command line, assemble your metagenome into contigs with your favorite assembler, I like MetaSPAdes. Next, map your metagenomic sequencing reads to the assembled contigs using BBMap and samtools. 

```bash
user@server:~$ reformat.sh app=t in=reads.R1.fastq.gz out=./appended_rawreads_for_mapping/allreads.fastq.gz
user@server:~$ reformat.sh app=t in=reads.R2.fastq.gz out=./appended_rawreads_for_mapping/allreads.fastq.gz
 
user@server:~$ bbmap.sh -Xmx400g threads=80 ambiguous=random qtrim=lr minid=0.97 nodisk=t ref=./assemblies/WT1_meta.fasta in1=./appended_rawreads_for_mapping/allreads.fastq.gz outm=./read_mapping/readmapping.bam 

user@server:~$ samtools sort readmapping.bam -o ./sorted_bams/readmapping_sorted.bam
user@server:~$ samtools index ./sorted_bams/readmapping_sorted.bam
```

Finally, generate read mapping summaries using Pileup, a function within the BBmap suite of bioinformatics tools. It is very important that you use a 'binsize' of 100 to ensure resolution is not lost by averaging over larger window sizes. The .bincov100 file created with Pileup can be used directly as input for ProActive.

```bash
user@server:~$  pileup.sh in=readmapping_sorted.bam out=readcoverage_summary.pileupcovstats bincov=readcoverage_summary.bincov100 binsize=100 stdev=t
```

### Install and run ProActive

##### Install ProActive:

```R
library(devtools)
install_github("jlmaier12/ProActive", ref="master")
library(ProActive)
```

##### Run ProActive:

Import your readcoverage_summary.bincov100 file:
```R
read_coverages <- read.delim("/readcoverage_summary.bincov100",  header=FALSE, comment.char="#")
```

Make active prophage identifications:
```R
active_prophages <- ProActive(read_coverages) 
```
The list created by prophage_activity_finder_func contains three objects. The first object is a table containing the summary prediction information for all the contigs in your metagenome. The second object is a cleaned table with the summary prediction information for only the contigs predicted as containing active prophages. The third object is a list containing the 'shape-match' information that ProActive used to make its identifications. ProActive uses the information within the third list object to rebuild the shapes it used to identify active prophages for plotting purposes. 

Plot mapped read coverage plots of contigs predicted as containing an active prophage. The 'shape' ProActive used to make its prediction for each contig will be overlayed on the barplot.:
```R
ProActivePredictionPlots(read_coverages, active_prophages[[2]], active_prophages[[1]]) 
```

##### Notes:
The mapped read coverages in your input file are averaged over 100 base pair windows, however ProActive can re-size the window size to both improve processing speed and reduce noise in the data. The larger the window size, the faster ProActive runs, but with larger window sizes, read coverage resolution is lost. Users may choose to re-size the windows to 200, 500, 1000, or 2000 base pairs. 1000 base pairs is the default window size for ProActive. If a window size other than 1000 is choosen, the user must specify this number as the second and fourth argument in the prophage_activity_finder_func and the ReadCoverageGraphs_activeprophage_func, respectively.  

