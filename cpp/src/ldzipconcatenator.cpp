#include "ldzipconcatenator.hpp"
#include "metadata.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <zstd.h>

namespace ldzip {

// Helper: binary copy file with stream buffer
template<typename T>
void copy_file_with_offset(const std::string& src, T& dst) {
    std::ifstream in(src, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open: " + src);
    dst << in.rdbuf();
}

// Helper: adjust compressed index with byte and column offsets
void adjust_index(const std::string& idx_file,
                  std::vector<uint64_t>& byte_offsets,
                  std::vector<uint64_t>& column_counts,
                  uint64_t cumulative_bytes,
                  uint64_t cumulative_cols,
                  bool skip_first) {
    std::ifstream in(idx_file, std::ios::binary | std::ios::ate);
    if (!in) throw std::runtime_error("Cannot open: " + idx_file);

    size_t compressed_size = in.tellg();
    in.seekg(0);
    std::vector<uint8_t> compressed(compressed_size);
    in.read(reinterpret_cast<char*>(compressed.data()), compressed_size);

    unsigned long long decompressed_size = ZSTD_getFrameContentSize(compressed.data(), compressed_size);
    if (decompressed_size == ZSTD_CONTENTSIZE_ERROR || decompressed_size == ZSTD_CONTENTSIZE_UNKNOWN)
        throw std::runtime_error("Cannot determine index size");

    std::vector<uint8_t> decompressed(decompressed_size);
    size_t result = ZSTD_decompress(decompressed.data(), decompressed_size, compressed.data(), compressed_size);
    if (ZSTD_isError(result))
        throw std::runtime_error("Index decompression failed: " + std::string(ZSTD_getErrorName(result)));

    size_t n_boundaries = decompressed_size / sizeof(uint64_t) / 2;
    uint64_t* idx_ptr = reinterpret_cast<uint64_t*>(decompressed.data());

    size_t start = skip_first ? 1 : 0;
    for (size_t i = start; i < n_boundaries; i++) {
        byte_offsets.push_back(idx_ptr[i] + cumulative_bytes);
        column_counts.push_back(idx_ptr[n_boundaries + i] + cumulative_cols);
    }
}

// Helper: compress and write index
void write_index(const std::string& path,
                 const std::vector<uint64_t>& byte_offsets,
                 const std::vector<uint64_t>& column_counts) {
    size_t n = byte_offsets.size();
    size_t raw_size = 2 * n * sizeof(uint64_t);
    std::vector<uint8_t> raw(raw_size);
    uint64_t* ptr = reinterpret_cast<uint64_t*>(raw.data());
    std::memcpy(ptr, byte_offsets.data(), n * sizeof(uint64_t));
    std::memcpy(ptr + n, column_counts.data(), n * sizeof(uint64_t));

    size_t bound = ZSTD_compressBound(raw_size);
    std::vector<uint8_t> compressed(bound);
    size_t compressed_size = ZSTD_compress(compressed.data(), bound, raw.data(), raw_size, 6);
    if (ZSTD_isError(compressed_size))
        throw std::runtime_error("Index compression failed");

    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(compressed.data()), compressed_size);
}

// Helper: validate and load metadata for naive concat
std::vector<MetaInfo> validate_and_load_metadata(const std::vector<std::string>& prefixes) {
    std::vector<MetaInfo> metas;
    size_t chunk_size = 0;
    for (size_t i = 0; i < prefixes.size(); i++) {
        MetaInfo meta = read_metadata_json(prefixes[i] + ".meta.json");
        if (meta.version != "3.0")
            throw std::runtime_error("Naive concat requires v3.0 format. Chunk " + std::to_string(i) + " is version " + meta.version);
        if (i == 0) chunk_size = meta.chunk_size;
        else if (meta.chunk_size != chunk_size)
            throw std::runtime_error("Chunk size mismatch");
        metas.push_back(meta);
    }
    return metas;
}

// Helper: copy binary file and adjust index
template<typename OutStreamT>
uint64_t copy_and_adjust_index(const std::string& bin_file,
                                const std::string& idx_file,
                                OutStreamT& out_stream,
                                std::vector<uint64_t>& byte_offsets,
                                std::vector<uint64_t>& column_counts,
                                uint64_t cumulative_bytes,
                                uint64_t cumulative_columns,
                                bool skip_first) {
    std::ifstream in(bin_file, std::ios::binary);
    in.seekg(0, std::ios::end);
    uint64_t bin_size = in.tellg();
    in.seekg(0);
    out_stream << in.rdbuf();
    adjust_index(idx_file, byte_offsets, column_counts, cumulative_bytes, cumulative_columns, skip_first);
    return bin_size;
}

// State for naive concat
struct NaiveConcatState {
    std::ofstream p_out, i_out;
    std::map<Stat, std::ofstream> x_outs;
    std::vector<uint64_t> i_byte_offsets, i_column_counts;
    std::map<Stat, std::vector<uint64_t>> x_byte_offsets, x_column_counts;
};

// Helper: initialize output streams and index arrays
NaiveConcatState initialize_naive_concat(const std::string& out_prefix, const MetaInfo& meta) {
    NaiveConcatState state;

    state.p_out.open(out_prefix + ".p.bin", std::ios::binary);
    state.i_out.open(out_prefix + ".i.bin", std::ios::binary);

    for (Stat s : All_Stats()) {
        if (meta.has_stat[s]) {
            state.x_outs[s].open(out_prefix + ".x." + stat_to_string(s) + ".bin", std::ios::binary);
            state.x_byte_offsets[s] = std::vector<uint64_t>();
            state.x_column_counts[s] = std::vector<uint64_t>();
        }
    }

    return state;
}

LDZipConcatenator::LDZipConcatenator(size_t nrows,
                                     size_t ncols,
                                     MatrixFormat format,
                                     const std::vector<Stat>& stats,
                                     Bits bits,
                                     const std::string& prefix,
                                     size_t chunk_size)
    : overlap_merger_(nrows, ncols, format, stats, bits, prefix, LDZipCompressor::Mode::ColumnStream, chunk_size) {

    // Set up large buffer for p_stream
    constexpr size_t buffer_size = 8 * 1024 * 1024;  // 8 MB buffer
    overlap_merger_.p_stream_.rdbuf()->pubsetbuf(nullptr, buffer_size);

    // Initialize p_file with 0
    write_initial_p_zero();
}

void LDZipConcatenator::write_initial_p_zero() {
    uint64_t zero = 0;
    overlap_merger_.p_stream_.write(reinterpret_cast<const char*>(&zero), sizeof(zero));
    overlap_merger_.p_stream_.flush();  // Ensure it's written
}

void LDZipConcatenator::init_p_for_column(uint32_t col, uint64_t value) {
    overlap_merger_.m_.p_[col] = value;
}

uint64_t LDZipConcatenator::get_p_value(uint32_t col) const {
    return overlap_merger_.m_.p_[col];
}

void LDZipConcatenator::close() {
    // Close all chunked writers (flush buffered data to disk)
    overlap_merger_.i_chunked_writer_->close();
    for (Stat s : overlap_merger_.m_.stats_available_) {
        overlap_merger_.x_chunked_writers_[s]->close();
    }
    overlap_merger_.p_stream_.close();

    // Set total nnz from tracked value
    overlap_merger_.m_.set_nnz(current_nnz_);

    // Write metadata after all data files are closed
    MetaInfo meta = overlap_merger_.m_.metaInfo();
    meta.chunk_size = overlap_merger_.chunk_size_;
    write_metadata_json(overlap_merger_.m_.metaFile(), meta);
}

void LDZipConcatenator::push_decoded_column(uint32_t cidx) {
    EnumArray<float, Stat> Statvalues(-999.0f);
    for (size_t k = 0; k < merged_rows_.size(); ++k) {
        for (Stat s : All_Stats())
            if (overlap_merger_.m_.has_stat_[s])
                Statvalues[s] = merged_x_[s][k];

        overlap_merger_.push_value_(merged_rows_[k], cidx, Statvalues);
    }
    overlap_merger_.active_column_ = cidx;
    writeMergedColumn();
}

void LDZipConcatenator::writeMergedColumn() {
    // Update p_vector before writing out column
    overlap_merger_.m_.p_[overlap_merger_.active_column_ + 1] =
        overlap_merger_.m_.p_[overlap_merger_.active_column_] +
        overlap_merger_.m_.i_[overlap_merger_.active_column_].size();

    // Write p-value immediately to stream (already in absolute coordinates)
    uint64_t p_val = overlap_merger_.m_.p_[overlap_merger_.active_column_ + 1];
    overlap_merger_.p_stream_.write(reinterpret_cast<const char*>(&p_val), sizeof(p_val));

    // Write i and x data
    for (Stat s : All_Stats()) if (overlap_merger_.m_.has_stat_[s]) overlap_merger_.write_x(s);
    overlap_merger_.write_i();
}

void LDZipConcatenator::clear_merge_buffers(const LDZipMatrix& mat) {
    above_rows_.clear(); above_idx_.clear();
    below_rows_.clear(); below_idx_.clear();
    merged_rows_.clear();

    for (Stat s : All_Stats()) {
        if (mat.has_stat(s)) {
            above_x_[s].clear();
            below_x_[s].clear();
            merged_x_[s].clear();
        }
    }
}

void LDZipConcatenator::extract_above_column(const LDZipMatrix& mat, uint32_t local_col, uint32_t threshold, uint32_t global_offset) {
    auto rows = mat.get_i(local_col);
    for (uint32_t k = 0; k < rows.size(); k++) {
        if (rows[k] <= threshold) {
            above_rows_.push_back(rows[k] + global_offset);
            above_idx_.push_back(k);
        }
    }

    for (Stat s : All_Stats())
        if (mat.has_stat(s)) {
            auto x = mat.get_x(local_col, s);
            for (uint32_t k : above_idx_) above_x_[s].push_back(x[k]);
        }
}

void LDZipConcatenator::extract_below_column(const LDZipMatrix& mat, uint32_t local_col, uint32_t threshold, uint32_t global_offset) {
    auto rows = mat.get_i(local_col);
    for (uint32_t k = 0; k < rows.size(); k++) {
        if (rows[k] > threshold) {
            below_rows_.push_back(rows[k] + global_offset);
            below_idx_.push_back(k);
        }
    }

    for (Stat s : All_Stats())
        if (mat.has_stat(s)) {
            auto x = mat.get_x(local_col, s);
            for (uint32_t k : below_idx_) below_x_[s].push_back(x[k]);
        }
}

void LDZipConcatenator::extract_full_column(const LDZipMatrix& mat, uint32_t local_col, uint32_t global_offset) {
    auto rows = mat.get_i(local_col);
    for (uint32_t k = 0; k < rows.size(); k++) {
        merged_rows_.push_back(rows[k] + global_offset);
    }

    for (Stat s : All_Stats())
        if (mat.has_stat(s)) {
            auto x = mat.get_x(local_col, s);
            for (uint32_t k = 0; k < x.size(); k++)
                merged_x_[s].push_back(x[k]);
        }
}

void LDZipConcatenator::merge_above_below_buffers(const LDZipMatrix& mat) {
    merged_rows_.insert(merged_rows_.end(), above_rows_.begin(), above_rows_.end());
    merged_rows_.insert(merged_rows_.end(), below_rows_.begin(), below_rows_.end());

    for (Stat s : All_Stats())
        if (mat.has_stat(s)) {
            merged_x_[s].insert(merged_x_[s].end(), above_x_[s].begin(), above_x_[s].end());
            merged_x_[s].insert(merged_x_[s].end(), below_x_[s].begin(), below_x_[s].end());
        }
}

void LDZipConcatenator::process_overlap_columns(const LDZipMatrix& current_mat,
                                                 const LDZipMatrix& next_mat,
                                                 const ChunkBoundary& current_chunk,
                                                 const ChunkBoundary& next_chunk) {
    uint32_t gs_current = current_chunk.global_start;
    uint32_t gs_next = next_chunk.global_start;

    // Initialize p_[] for the first overlapping column
    uint32_t first_overlap_col = gs_current + current_chunk.second_overlap_start;
    init_p_for_column(first_overlap_col, current_nnz_);

    for (uint32_t lc_current = current_chunk.second_overlap_start; lc_current < current_chunk.n_variants; lc_current++) {
        uint32_t global_col = gs_current + lc_current;
        uint32_t lc_next    = lc_current - current_chunk.second_overlap_start;

        clear_merge_buffers(current_mat);

        // Above/diagonal from current chunk (row <= lc_current)
        extract_above_column(current_mat, lc_current, lc_current, gs_current);

        // Below diagonal from next chunk (row > lc_next)
        extract_below_column(next_mat, lc_next, lc_next, gs_next);

        // Merge above + below
        merge_above_below_buffers(current_mat);

        push_decoded_column(global_col);
    }

    // Calculate overlapping columns nnz and track it
    uint32_t last_overlap_col = gs_current + current_chunk.n_variants;
    uint64_t overlap_nnz = get_p_value(last_overlap_col) - get_p_value(first_overlap_col);
    current_nnz_ += overlap_nnz;
}

void LDZipConcatenator::process_exclusive_columns(const LDZipMatrix& current_mat,
                                                   const ChunkBoundary& current_chunk) {
    size_t lc_start = current_chunk.first_overlap_end;
    size_t lc_end   = current_chunk.second_overlap_start;

    if (lc_start >= lc_end) return;  // No exclusive columns

    uint32_t gs_current = current_chunk.global_start;

    // Initialize p_[] for first exclusive column if not already set
    uint32_t first_excl_col = gs_current + lc_start;
    if (overlap_merger_.m_.p_[first_excl_col] == 0) {
        init_p_for_column(first_excl_col, current_nnz_);
    }

    // Decode and write each exclusive column
    for (size_t lc = lc_start; lc < lc_end; lc++) {
        uint32_t global_col = gs_current + lc;

        clear_merge_buffers(current_mat);

        // Extract full column with global offset
        extract_full_column(current_mat, lc, gs_current);

        push_decoded_column(global_col);
    }

    uint64_t excl_nnz = overlap_merger_.m_.p_[gs_current + lc_end] -
                       overlap_merger_.m_.p_[first_excl_col];
    current_nnz_ += excl_nnz;
}

namespace {
bool read_next_variant_id(std::ifstream& in, std::string& id) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t first_tab = line.find('\t');
        if (first_tab != std::string::npos) {
            std::istringstream ss(line);
            std::string chrom, pos;
            if (std::getline(ss, chrom, '\t') && std::getline(ss, pos, '\t') && std::getline(ss, id, '\t'))
                return true;
        }
        id = line;
        return true;
    }
    return false;
}
} // anonymous namespace

OverlapVariantInfo read_overlapping_variant_order(const std::vector<std::string>& var_files, bool /* check_overlap */) {
    if (var_files.empty()) throw std::runtime_error("No variant files provided");

    size_t n = var_files.size();
    OverlapVariantInfo info;
    info.chunks.resize(n);

    // --- Single file: trivial case ---
    if (n == 1) {
        std::ifstream f(var_files[0]);
        if (!f) throw std::runtime_error("Cannot open variant file: " + var_files[0]);
        std::string id;
        size_t count = 0;
        while (read_next_variant_id(f, id)) ++count;
        info.chunks[0] = {count, 0, count, 0, count};
        info.total_variants = count;
        return info;
    }

    std::ifstream first(var_files[0]);
    if (!first) throw std::runtime_error("Cannot open variant file: " + var_files[0]);
    std::ifstream second(var_files[1]);
    if (!second) throw std::runtime_error("Cannot open variant file: " + var_files[1]);

    std::vector<size_t> local_counts(n, 0);
    size_t global_cursor = 0;

    info.chunks[0].global_start      = 0;
    info.chunks[0].first_overlap_end = 0;

    for (size_t i = 0; i + 1 < n; i++) {
        // Step 1: read pivot = first variant of second (chunk[i+1])
        std::string pivot;
        if (!read_next_variant_id(second, pivot))
            throw std::runtime_error("Variant file is empty: " + var_files[i + 1]);
        ++local_counts[i + 1];

        // Step 2: read first (chunk[i]) until we reach the pivot
        std::string v;
        bool found_pivot = false;
        while (read_next_variant_id(first, v)) {
            ++local_counts[i];
            if (v == pivot) {
                found_pivot = true;
                break;
            }
        }

        size_t overlap_start_local_i;
        size_t new_variants_of_i;
        size_t overlap_count;

        if (found_pivot) {
            // v == pivot; local_counts[i] now includes the pivot itself
            overlap_start_local_i = local_counts[i] - 1;  // 0-based index of pivot in chunk[i]
            new_variants_of_i     = overlap_start_local_i - info.chunks[i].first_overlap_end;

            info.chunks[i + 1].global_start = global_cursor + new_variants_of_i;

            // Step 3: overlap loop — both streams are now positioned just past the pivot;
            // read one from each and verify they match until first is exhausted
            overlap_count = 1;  // pivot already matched
            std::string v0, v1;
            while (read_next_variant_id(first, v0)) {
                ++local_counts[i];
                if (!read_next_variant_id(second, v1))
                    throw std::runtime_error(
                        "File '" + var_files[i + 1] + "' exhausted during overlap with '" +
                        var_files[i] + "'");
                ++local_counts[i + 1];
                if (v0 != v1)
                    throw std::runtime_error(
                        "Overlap mismatch between '" + var_files[i] + "' and '" +
                        var_files[i + 1] + "': '" + v0 + "' vs '" + v1 + "'");
                ++overlap_count;
            }
        } else {
            // No overlap found - reset second file stream to beginning
            overlap_count = 0;
            overlap_start_local_i = local_counts[i];
            new_variants_of_i = local_counts[i] - info.chunks[i].first_overlap_end;
            info.chunks[i + 1].global_start = global_cursor + new_variants_of_i;
            // Reset second stream since we read the pivot but didn't match it
            second.close();
            second.open(var_files[i + 1]);
            if (!second) throw std::runtime_error("Cannot re-open variant file: " + var_files[i + 1]);
            local_counts[i + 1] = 0;
        }

        std::cout << " Overlap between chunk " << i << " and " << (i + 1)
                  << ": " << overlap_count << " variants\n";

        // Finalize chunk[i]
        info.chunks[i].n_variants           = local_counts[i];
        info.chunks[i].second_overlap_start = overlap_start_local_i;
        info.chunks[i].global_end           = global_cursor + new_variants_of_i + overlap_count;

        info.chunks[i + 1].first_overlap_end = local_counts[i + 1];  // == overlap_count

        global_cursor += new_variants_of_i + overlap_count;

        // Swap: first <- second (positioned after the overlap), open next file into second
        first = std::move(second);
        if (i + 2 < n) {
            second.open(var_files[i + 2]);
            if (!second) throw std::runtime_error("Cannot open variant file: " + var_files[i + 2]);
        }
    }

    // After the loop, first points to var_files[n-1] positioned right after the last overlap.
    // Drain the remaining tail.
    {
        std::string id;
        while (read_next_variant_id(first, id)) ++local_counts[n - 1];
    }

    size_t tail = local_counts[n - 1] - info.chunks[n - 1].first_overlap_end;
    info.chunks[n - 1].n_variants           = local_counts[n - 1];
    info.chunks[n - 1].second_overlap_start = local_counts[n - 1];  // no next overlap
    info.chunks[n - 1].global_end           = global_cursor + tail;
    info.total_variants                     = info.chunks[n - 1].global_end;

    return info;
}

void LDZipConcatenator::concat_naive(const std::vector<std::string>& prefixes, const std::string& out_prefix) {
    if (prefixes.empty()) throw std::runtime_error("No input chunks provided");

    std::vector<MetaInfo> metas = validate_and_load_metadata(prefixes);

    // Check for overlaps using existing function
    std::vector<std::string> var_files;
    for (const auto& pref : prefixes) var_files.push_back(pref + ".vars.txt");

    OverlapVariantInfo ov = read_overlapping_variant_order(var_files, false);
    for (size_t i = 0; i < ov.chunks.size(); i++) {
        if (ov.chunks[i].first_overlap_end > 0 || ov.chunks[i].second_overlap_start < ov.chunks[i].n_variants)
            throw std::runtime_error("Naive concat detected overlap in chunk " + std::to_string(i) + ". Use default mode (without --naive)");
    }
    size_t total_variants = ov.total_variants;

    std::cout << "Naive concat: " << prefixes.size() << " chunks, " << total_variants << " variants\n";

    NaiveConcatState state = initialize_naive_concat(out_prefix, metas[0]);

    uint64_t cumulative_nnz = 0, cumulative_i_bytes = 0, cumulative_columns = 0;
    std::map<Stat, uint64_t> cumulative_x_bytes;
    for (Stat s : All_Stats()) {
        if (metas[0].has_stat[s]) cumulative_x_bytes[s] = 0;
    }

    // Write p[0] = 0
    state.p_out.write(reinterpret_cast<const char*>(&cumulative_nnz), sizeof(uint64_t));

    // Concatenate chunks
    for (size_t chunk_idx = 0; chunk_idx < prefixes.size(); chunk_idx++) {
        const std::string& prefix = prefixes[chunk_idx];
        const MetaInfo& meta = metas[chunk_idx];

        std::cout << " Concatenating chunk " << chunk_idx << ": " << prefix << "\n";

        // Adjust and write p values
        std::ifstream p_in(prefix + ".p.bin", std::ios::binary);
        std::vector<uint64_t> p_vals(meta.cols + 1);
        p_in.read(reinterpret_cast<char*>(p_vals.data()), (meta.cols + 1) * sizeof(uint64_t));
        for (size_t c = 1; c <= meta.cols; c++) {
            uint64_t adjusted = p_vals[c] + cumulative_nnz;
            state.p_out.write(reinterpret_cast<const char*>(&adjusted), sizeof(uint64_t));
        }

        // Binary copy and adjust indices
        cumulative_i_bytes += copy_and_adjust_index(prefix + ".i.bin", prefix + ".i.bin.index", state.i_out, state.i_byte_offsets, state.i_column_counts, cumulative_i_bytes, cumulative_columns, chunk_idx > 0);
        for (Stat s : All_Stats()) {
            if (!meta.has_stat[s]) continue;
            std::string x_file = prefix + ".x." + stat_to_string(s) + ".bin";
            cumulative_x_bytes[s] += copy_and_adjust_index(x_file, x_file + ".index", state.x_outs[s], state.x_byte_offsets[s], state.x_column_counts[s], cumulative_x_bytes[s], cumulative_columns, chunk_idx > 0);
        }

        cumulative_nnz += meta.nnz;
        cumulative_columns += meta.cols;
    }

    state.p_out.close();
    state.i_out.close();
    for (Stat s : All_Stats()) {
        if (metas[0].has_stat[s]) state.x_outs[s].close();
    }

    // Write index files
    write_index(out_prefix + ".i.bin.index", state.i_byte_offsets, state.i_column_counts);
    for (Stat s : All_Stats()) {
        if (metas[0].has_stat[s])
            write_index(out_prefix + ".x." + stat_to_string(s) + ".bin.index",
                        state.x_byte_offsets[s], state.x_column_counts[s]);
    }

    // Concatenate vars.txt
    std::ofstream vars_out(out_prefix + ".vars.txt");
    for (size_t i = 0; i < var_files.size(); i++) {
        std::ifstream f(var_files[i]);
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') {
                if (i == 0) vars_out << line << '\n';
            } else {
                vars_out << line << '\n';
            }
        }
    }

    // Write metadata
    MetaInfo out_meta(total_variants, total_variants, cumulative_nnz,
                     bits_to_int(metas[0].bits), metas[0].format, "3.0", metas[0].chunk_size);
    out_meta.has_stat = metas[0].has_stat;
    write_metadata_json(out_prefix + ".meta.json", out_meta);

    std::cout << "Naive concat complete: " << out_prefix << "\n";
}

} // namespace ldzip
