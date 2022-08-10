test_that("getEnrichedPath works not correct", {
  data(annotatedPeak)
  
  enriched.PATH <- getEnrichedPATH(annotatedPeak[1:500],
                                   orgAnn="org.Hs.eg.db", 
                                   feature_id_type="ensembl_gene_id",
                                   pathAnn="reactome.db", maxP=0.01,
                                   minPATHterm=10,
                                   multiAdjMethod=NULL)
  ann <- addGeneIDs(enriched.PATH[,2],
                    feature_id_type = "entrez_id",
                    orgAnn = org.Hs.eg.db,
                    IDs2Add = "symbol")
  enriched.PATH <- merge(ann, enriched.PATH)
  expect_true(nrow(enriched.PATH)>0)
})