#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test LDZip Tabular Concat (Non-Overlapping) \033[0m"

BITS_LIST=(8 16 32 99)
MIN_LIST=(0.0 0.01 0.1)

get_threshold() {
  local bits="$1"
  local min="$2"
  awk -v b="$bits" -v m="$min" '$1 == b && $2 == m {print $3}' assets/thresholds.tsv
}


BITS_LIST=(8 16 32 99)
MIN_LIST=(0.01 0.1)


for min in "${MIN_LIST[@]}"; do
    echo "➤ Generating random plink matrices ..."
    min_r2=$(echo "$min * $min" | bc -l)
    for chr in 20 21 22; do
        ${PLINK2} \
        --pfile ../../assets/g1k.chr${chr} \
        --ld-window-kb 100 \
        --ld-window-r2 $min_r2 \
        --r-unphased ref-based cols=id,ref,alt \
        --out $OUT_DIR/chr${chr} > /dev/null
    done

    for bits in "${BITS_LIST[@]}"; do
    echo
    echo -e "\033[1;33m>>> Testing bits=$bits, min=$min \033[0m"

    echo "➤ Compressing..."
    for chr in 20 21 22; do
        ../bin/ldzip compress plinkTabular \
        --ld_file $OUT_DIR/chr${chr}.vcor \
        --snp_file ../../assets/g1k.chr${chr}.pvar \
        --output_prefix $OUT_DIR/compressed.chr${chr} \
        --min $min \
        --bits $bits \
        --min_col UNPHASED_R > /dev/null
    done

    # Test naive mode
    echo "➤ Concating (naive mode - non-overlapping)..."
    ../bin/ldzip concat \
    --inputs $OUT_DIR/compressed.chr20 $OUT_DIR/compressed.chr21 $OUT_DIR/compressed.chr22 \
    --output_prefix $OUT_DIR/concat_naive \
    --naive > /dev/null

    echo "➤ Decompressing naive concat..."
    ../bin/ldzip decompress \
    --input_prefix $OUT_DIR/concat_naive \
    --output_prefix $OUT_DIR/concat_naive_full \
    --type tabular > /dev/null

    echo "➤ Comparing naive roundtrip..."
    Rscript ../../scripts/corr.R \
        -orig $OUT_DIR/chr20.vcor $OUT_DIR/chr21.vcor $OUT_DIR/chr22.vcor\
        -new $OUT_DIR/concat_naive_full.vcor > /dev/null

    max_diff=$(jq .max_diff "$OUT_DIR/concat_naive_full_stats.json")
    max_allowed=$(get_threshold "$bits" "$min")

    echo "➤ Check naive max_diff[$max_diff] <= threshold[$max_allowed]"
    if awk "BEGIN {exit !($max_diff <= $max_allowed)}"; then
        echo "✅ Naive mode within threshold"
    else
        echo "❌ Naive mode exceeds threshold"
        exit 1
    fi

    # Test default mode (decompress/merge/recompress)
    echo "➤ Concating (default mode - non-overlapping)..."
    ../bin/ldzip concat \
    --inputs $OUT_DIR/compressed.chr20 $OUT_DIR/compressed.chr21 $OUT_DIR/compressed.chr22 \
    --output_prefix $OUT_DIR/concat_default > /dev/null

    echo "➤ Decompressing default concat..."
    ../bin/ldzip decompress \
    --input_prefix $OUT_DIR/concat_default \
    --output_prefix $OUT_DIR/concat_default_full \
    --type tabular > /dev/null

    echo "➤ Comparing default roundtrip..."
    Rscript ../../scripts/corr.R \
        -orig $OUT_DIR/chr20.vcor $OUT_DIR/chr21.vcor $OUT_DIR/chr22.vcor\
        -new $OUT_DIR/concat_default_full.vcor > /dev/null

    max_diff=$(jq .max_diff "$OUT_DIR/concat_default_full_stats.json")

    echo "➤ Check default max_diff[$max_diff] <= threshold[$max_allowed]"
    if awk "BEGIN {exit !($max_diff <= $max_allowed)}"; then
        echo "✅ Default mode within threshold"
    else
        echo "❌ Default mode exceeds threshold"
        exit 1
    fi
  done
done
