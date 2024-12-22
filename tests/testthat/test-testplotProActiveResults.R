

test_that("plotProActiveResults works on metagenomes with defaults", {
  ##metagenome
  TestPlotsMetagenome <- plotProActiveResults(sampleMetagenomePileup, sampleMetagenomeResults)
  testPlot <- TestPlotsMetagenome$NODE_1884
  expect_doppelganger("Test NODE 1884 plotting", testPlot)
})

test_that("plotProActiveResults works on genomes with defaults", {
   ##genome
   TestPlotsGenome <- plotProActiveResults(sampleGenomePileup, sampleGenomeResults)
   testPlot2 <- TestPlotsGenome$NC_003197.2_chunk_10
   expect_doppelganger("Test genome chunk 10 for plotting", testPlot2)
 })
