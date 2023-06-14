# ProActive
Active prophage detection in metagenome sequencing read coverages

### Decription
ProActive detects regions of elevated read coverage representative of active propgages in sequencing reads mapped to metagenome contigs. When a prophage activates and enters the lytic cycle, its genome begins replicating and the ratio of bacterial:phage genome in the cell begins to decrease. Because there are more phage genomes than bacterial genomes, the phage genome is overrepresented during sequencing and more reads are generated for the phage than the bacteria. When reads are mapped back to the contigs, the phage reads will map to its prophage location within the bacterial genome and due to the overabundance of phage reads, the prophage region read coverage will be elevated in comparison to the read coverage of the bacterial genome on either side of the prophage. The increase in read coverage at prophage locations can be detected and inferences can be made about the lytic activity of the associated phage. However, in order to assess variations in read coverage between the phage and bacterial genomes, the genomic coordinates of the prophages must be known. Genomic locations of prophages are usually idenfitied, at least in part, using phage-related gene annotations. Relying on gene annotations for prophage detection means that only annotated phage are detcted while the vast majority are missed. ProActive, however, is reference-independent and therefore bypasses this barier to active prophage identification.

ProActive works by detecting the 'shape' created by the read coverages of an active prophage. When the read coverages for an active prophage are plotted as a bar graph, the elevated read coverage at the prophage location creates a 'block-like' shape. ProActive scans through all the contigs in a metagenome and searches each one for a 'block-like' shape with a base pair size that is realistic for that of a phage genome. ProActive will output a summary table containing all the conitgs with a potentially active prophage, the sizes of the active prophages, and bar graphs visualizing the active prophage read coverage with the shape that ProActive used to make its final identification. As an optional feature, users may input a PROKKA gff file containing gene annotations for their contigs which ProActive will use to determine which of the identified active prophages have phage-related gene annotations. 

To be clear, while ProActive labels the elevated read coverage it detcts as active prophage, there are other circumstances that can create read coverage patterns similar to that of active prophage. For example, chimeric assemblies can create contigs with regions of elevated read coverage that are not associated with active prophage, but ProActive will likely identify this type of event as an active prophage. It is up to the user to assess the active prophage identifications and check for false positives. Typically, there tendd to be more false positive identifications made in short, low-quality contigs. 


### Data input for ProActive
ProActive requires one input file containg read coverage summary information to run. We recommend the following workflow to get your metagenomic sequence data into the correct format for input to ProActive. 

On the command line, assemble your metagenome into contigs with your favorite assembler, like MetaSPAdes. Then map the metagenomic sequencing reads to the assembled contigs using BBMap and samtools. 

```bash
user@server:~$ reformat.sh app=t in=reads.R1.fastq.gz out=./appended_rawreads_for_mapping/allreads.fastq.gz
user@server:~$ reformat.sh app=t in=reads.R2.fastq.gz out=./appended_rawreads_for_mapping/allreads.fastq.gz
 
user@server:~$ bbmap.sh -Xmx400g threads=80 ambiguous=random qtrim=lr minid=0.97 nodisk=t ref=./assemblies/WT1_meta.fasta in1=./appended_rawreads_for_mapping/allreads.fastq.gz outm=./read_mapping/readmapping.bam 

user@server:~$ samtools sort readmapping.bam -o ./sorted_bams/readmapping_sorted.bam
user@server:~$ samtools index ./sorted_bams/readmapping_sorted.bam
```

Finally, generate read mapping summaries using Pileup, a function within the BB suite of bioinformatics tools. It is very important that you use a 'binsize' of 100 to ensure resolution is not lost by averaging over larger window sizes. The .bincov100 file created with Pileup can be used directly as input for ProActive.

```bash
user@server:~$  pileup.sh in=./sorted_bams/readmapping_sorted.bam out=./coverage_stats/readcoverage_summary.pileupcovstats bincov=./coverage_stats/readcoverage_summary.bincov100 binsize=100 stdev=t
```

### Install and run ProActive

##### Install ProActive:

```R
library(devtools)
install_github("jlmaier12/ProActive")
library(ProActive)
```

##### Run ProActive:

Import your readcoverage_summary.bincov100 file:
```R
read_coverages <- read.delim("/readcoverage_summary.bincov100",  header=FALSE, comment.char="#")
```

Make active prophage identifications:
```R
active_prophages <- prophage_activity_finder_func(read_coverages) 
```
The object created by prophage_activity_finder_func is a list containing three objects. The first object is a table contanining the summary information for all the contigs in your metagenome. The second object is a table with the summary information for only the contigs predicted as containing active prophages. The third object is a list containing the 'shape-match' information that ProActive used to make its identifications. ProActive uses this information to rebuild the shape it used to identify an active prophage for plotting purposes. 

Plot read coverage barplots of contigs predicted as containing an active prophage. The 'shape' ProActive used to make its identification for each contig will be overlayed on the barplot.:
```R
shape_matching_plots_WC(read_coverages, active_prophages[[3]], active_prophages[[1]]) 
```
