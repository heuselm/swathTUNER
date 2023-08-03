
context("Functions for building isolation windows")

test_that("FUNCTION buildFixedWidthWindows", 
{
    expect_is(buildFixedWidthWindows(2,100,300), "data.frame")
    expect_equal(buildFixedWidthWindows(2,100,300), data.frame(start=c(100,200), end=c(200,300)), tolerance = 0.01)
})


test_that("FUNCTION buildWindowsEqualizedPrecFreq", 
{
    in1 = read.table("../../inst/extdata/LC-MS_peaks.tsv", header = TRUE, sep = "\t")
    expout1 = read.table("../../inst/extdata/windows_equalizePrecFreq.tsv", header = TRUE, sep = "\t")
    out1 = buildWindowsEqualizedPrecFreq(in1,36,350,1250)
    expect_is(out1, "data.frame")
    expect_true(all.equal(addOverlap(out1,1), expout1, tolerance = 0.001))
})


test_that("FUNCTION buildWindowsEqualizedTIC", 
{
    in1 = read.table("../../inst/extdata/LC-MS_peaks.tsv", header = TRUE, sep = "\t")
    expout1 = read.table("../../inst/extdata/windows_equalizeTIC.tsv", header = TRUE, sep = "\t")
    out1 = buildWindowsEqualizedTIC(in1,36,350,1250)
    expect_is(out1, "data.frame")
    expect_true(all.equal(addOverlap(out1,1), expout1, tolerance = 0.001))
})


test_that("FUNCTION addOverlap", 
{
    expect_equal(addOverlap(buildFixedWidthWindows(2,100,300), 0), data.frame(start=c(100,200), end=c(200,300)), tolerance = 0)
    expect_equal(addOverlap(buildFixedWidthWindows(2,100,300), 1), data.frame(start=c(99.5,200.5), end=c(199.5,300.5)), tolerance = 0.01)
})
