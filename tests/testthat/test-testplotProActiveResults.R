test_that("plotProActiveResults works on genomes with defaults", {
  ##genome
  TestPlotsGenome <- plotProActiveResults(exampleGenomePileupSubset, exampleGenomeSubsetResults)
  expect_equal(TestPlotsGenome, examplePlotsGenomeSubset)
  })

test_that("plotProActiveResults works on metagenomes with defaults", {
  ##metagenome
  TestPlotsMetagenome <- plotProActiveResults(sampleMetagenomePileup, sampleMetagenomeResults)
  expect_equal(TestPlotsMetagenome, samplePlotsMetagenome)
})
