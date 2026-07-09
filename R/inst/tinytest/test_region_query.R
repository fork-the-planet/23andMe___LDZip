library(LDZipMatrix)
library(tinytest)

source("../../../scripts/source.R")

# load original matrix
orig <- read_tabular_ld_matrix("../../tests/data/test_tabular")

# list of compressed outputs
bits_list <- c(8, 16, 32, 99)
min_list  <- c(0.00, 0.001, 0.01, 0.1)

for (bits in bits_list) {
  for (min in min_list) {
    tol <- get_threshold(bits, min)

    prefix <- sprintf("../../tests/data/compressed_bits_tabular_%d_min_%s", bits, min)
    cat("Testing:", prefix, "\n")

    # load compressed
    ld <- LDZipMatrix(prefix)
    suppressMessages(buildIndex(ld))

    # Test region queries with and without "chr" prefix
    rsids_expected <- paste0("rs", 1:100)
    val_rsid <- fetchLD(ld, rsids_expected, rsids_expected, types="UNPHASED_R", simplify=F)$UNPHASED_R
    val_orig <- orig[rsids_expected, rsids_expected]

    # Test without "chr" prefix
    region_no_chr <- "2:10000-10099"
    val_region_no_chr <- fetchLD(ld, region_no_chr, region_no_chr, types="UNPHASED_R", simplify=F)$UNPHASED_R

    expect_equal(
        as.vector(val_region_no_chr),
        as.vector(val_rsid),
        tolerance = 1e-10,
        scale = 1,
        info = sprintf("Region query without chr should match rsID query for %s", prefix)
    )

    expect_equal(
        as.vector(val_region_no_chr),
        as.vector(val_orig),
        tolerance = tol + 1e-7,
        scale = 1,
        info = sprintf("Region query without chr should match original matrix for %s", prefix)
    )

    # Test with "chr" prefix
    region_with_chr <- "chr2:10000-10099"
    val_region_with_chr <- fetchLD(ld, region_with_chr, region_with_chr, types="UNPHASED_R", simplify=F)$UNPHASED_R

    expect_equal(
        as.vector(val_region_with_chr),
        as.vector(val_rsid),
        tolerance = 1e-10,
        scale = 1,
        info = sprintf("Region query with chr should match rsID query for %s", prefix)
    )

    expect_equal(
        as.vector(val_region_with_chr),
        as.vector(val_orig),
        tolerance = tol + 1e-7,
        scale = 1,
        info = sprintf("Region query with chr should match original matrix for %s", prefix)
    )

    # Test smaller region (rs50-rs60)
    region_small <- "2:10049-10059"
    rsids_small <- paste0("rs", 50:60)

    val_region_small <- fetchLD(ld, region_small, region_small, types="UNPHASED_R", simplify=F)$UNPHASED_R
    val_rsid_small <- fetchLD(ld, rsids_small, rsids_small, types="UNPHASED_R", simplify=F)$UNPHASED_R

    expect_equal(
        as.vector(val_region_small),
        as.vector(val_rsid_small),
        tolerance = 1e-10,
        scale = 1,
        info = sprintf("Small region query should match rsID query for %s", prefix)
    )

    # Test asymmetric region query (different row vs col regions)
    region_row <- "2:10000-10050"
    region_col <- "2:10051-10099"
    rsids_row <- paste0("rs", 1:51)
    rsids_col <- paste0("rs", 52:100)

    val_region_asym <- fetchLD(ld, region_row, region_col, types="UNPHASED_R", simplify=F)$UNPHASED_R
    val_rsid_asym <- fetchLD(ld, rsids_row, rsids_col, types="UNPHASED_R", simplify=F)$UNPHASED_R

    expect_equal(
        as.vector(val_region_asym),
        as.vector(val_rsid_asym),
        tolerance = 1e-10,
        scale = 1,
        info = sprintf("Asymmetric region query should match rsID query for %s", prefix)
    )

    # Test out-of-range region (should error)
    expect_error(
        fetchLD(ld, "2:99000000-99999999", "2:99000000-99999999", types="UNPHASED_R"),
        "No matching variants found in specified region",
        info = sprintf("Out-of-range region should error for %s", prefix)
    )
  }
}
