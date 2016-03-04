test_that("addMetadata works not correct", {
    expect_error(addMetadata(GRanges()))
})