test_that("ProActive works on metagenomes with defualts", {
  ##metagenome
  ProActiveTestMetagenome <- ProActiveDetect(sampleMetagenomePileup,
                                       mode="metagenome",
                                       gffTSV = sampleMetagenomegffTSV)
  expect_equal(ProActiveTestMetagenome, sampleMetagenomeResults)
})

test_that("ProActive works on genomes with defualts", {
  ##genome
  ProActiveTestGenome <- ProActiveDetect(sampleGenomePileup,
                                   mode="genome",
                                   gffTSV = sampleGenomegffTSV)
  expect_equal(ProActiveTestGenome, sampleGenomeResults)
})
