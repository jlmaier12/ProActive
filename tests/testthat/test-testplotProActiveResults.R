test_that("plotProActiveResults works on metagenomes with defaults", {
  ##metagenome
  TestPlotsMetagenome <- plotProActiveResults(sampleMetagenomePileup, sampleMetagenomeResults)
  expect_equal(TestPlotsMetagenome$NODE_1884, "samplePlot_NODE_1884")
})

# test_that("plotProActiveResults works on genomes with defaults", {
#   ##genome
#   TestPlotsGenome <- plotProActiveResults(sampleGenomePileup, sampleGenomeResults)
#   expect_equal(TestPlotsGenome$NC_003197.2_chunk_10, "samplePlot_NC_003197.2_chunk_10")
# })
