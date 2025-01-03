test_that("ProActive works on metagenomes with defualts", {
  ##metagenome
  ProActiveTestMetagenome <- ProActive(sampleMetagenomePileup, mode="metagenome", gffTSV = sampleMetagenomegffTSV)
  expect_equal(ProActiveTestMetagenome, sampleMetagenomeResults)
})

test_that("ProActive works on genomes with defualts", {
  ##genome
  ProActiveTestGenome <- ProActive(exampleGenomePileupSubset, mode="genome", gffTSV = exampleGenomegffTSVSubset)
  expect_equal(ProActiveTestGenome, exampleGenomeSubsetResults)
})
