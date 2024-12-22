
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ProActive

<!-- badges: start -->
<!-- badges: end -->

`ProActive` detects elevations and gaps in sequencing read coverage
using a pattern-matching algorithm. ProActive is best used as a
screening method to identify genetic regions for further investigation.

Elevations or gaps in read coverage can be caused by differential
abundance of specific genetic elements. For example, an elevation in
read coverage may be caused by prophage activation. When a prophage
activates and enters the lytic cycle, its genome begins replicating and
the ratio of phage:bacterial genomes in the cell begins to increase.
Because there are more phage genomes than bacterial genomes, during
sequencing more phage reads are generated than bacterial. When these
reads are mapped back to their associated reference sequence, the read
coverage of the prophage region will be elevated in comparison to the
read coverage of the bacterial genome on either side of the prophage.
This same principle applies to temperate phage who are highly abundant
in the environment as well as other mobile genetic elements that are
freely present in the environment at a higher ratio than the originating
or ‘host’ genome.

Conversely, a gap in read coverage may indicate genetic heterogeneity in
the associated bacterial population. Genetic variants with and without
specific genetic elements, like prophage or certain genes, will produce
differential abundances of sequencing reads that may form read coverage
gaps. The formation of read coverage gaps due to genetic variants is
dependent on the assembly (i.e. if the assembler assembles the genetic
variants as separate entities or not). Gaps in read coverage may also
form at regions with high mutation rates.

### Input files

#### Pileup file

ProActive detects read coverage patterns using a pattern-matching
algorithm that operates on pileup files. A pileup file is a file format
where each row summarizes the ‘pileup’ of reads at specific genomic
locations. Pileup files can be used to generate a rolling mean of read
coverages and associated base pair positions which reduces data size
while preserving read coverage patterns. **ProActive requires that input
pileups files** **be generated using a 100 bp window/bin size.**

Pileup files can be generated by either mapping sequencing reads to
either a metagenome or genome fasta. **Read mapping should be performed
using a high** **minimum identity (0.97 or higher) and random mapping of
ambiguous reads.** The pileup files needed for ProActive are generated
using the .bam files produced during read mapping.

Some read mappers, like
[BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/),
allow for the generation of pileup files in the
[`bbmap.sh`](https://github.com/BioInfoTools/BBMap/blob/master/sh/bbmap.sh)
command with use of the `bincov` output with the `covbinsize=100`
parameter/argument. **Otherwise, BBMap’s**
**[`pileup.sh`](https://github.com/BioInfoTools/BBMap/blob/master/sh/pileup.sh)**
**can convert .bam files produced by any read mapper to pileup files**
**compatible with ProActive using the `bincov` output with
`binsize=100`.**

#### gffTSV

ProActive optionally accepts a .gff file as input. The .gff file must be
associated with the same metagenome or genome used to create your pileup
file. The .gff file should be a TSV and should follow the same general
format described
[here](https://en.wikipedia.org/wiki/General_feature_format#:~:text=In%20bioinformatics%2C%20the%20general%20feature,DNA%2C%20RNA%20and%20protein%20sequences.).

## Installation

You can install the development version of ProActive from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jlmaier12/ProActive")
```

## Quick start

Metagenome mode:

``` r
library(ProActive)

MetagenomeProActive <- ProActive(
  pileup = sampleMetagenomePileup,
  mode = "metagenome",
  gffTSV = sampleMetagenomegffTSV
)
#> Preparing input file for pattern-matching...
#> Starting pattern-matching...
#> A quarter of the way done with pattern-matching
#> Half of the way done with pattern-matching
#> Almost done with pattern-matching!
#> Summarizing pattern-matching results
#> Finding ORFs in elevated or gapped regions of read coverage...
#> Finalizing output
#> Execution time: 3.06secs
#> 0 contigs were filtered out based on low read coverage
#> 0 contigs were filtered out based on length (< minContigLength)
#> 
#> Elevation       Gap NoPattern 
#>         3         3         1

MetagenomePlots <- plotProActiveResults(pileup = sampleMetagenomePileup,
                                        ProActiveResults = MetagenomeProActive)
```

Genome mode:

``` r
GenomeProActive <- ProActive(
  pileup = sampleGenomePileup,
  mode = "genome",
  gffTSV = sampleGenomegffTSV
)
#> Preparing input file for pattern-matching...
#> Starting pattern-matching...
#> A quarter of the way done with pattern-matching
#> Half of the way done with pattern-matching
#> Almost done with pattern-matching!
#> Summarizing pattern-matching results
#> Finding ORFs in elevated or gapped regions of read coverage...
#> Finalizing output
#> Execution time: 48.77secs
#> 0 contigs were filtered out based on low read coverage
#> 0 contigs were filtered out based on length (< minContigLength)
#> 
#> Elevation       Gap NoPattern 
#>        25         3        21

GenomePlots <- plotProActiveResults(pileup = sampleGenomePileup,
                                    ProActiveResults = GenomeProActive)
```
