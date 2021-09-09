test_that("enrichmentPlot works not correct", {
    data(enrichedGO)
    ## test enrichmentPlot can handle empty elements.
    enrichedGO$cc <- enrichedGO$cc[enrichedGO$cc$go.id=="A", ]
    enrichmentPlot(enrichedGO)
})