#pragma once
#include <zstd.h>
#include <vector>
#include <fstream>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <memory>
#include <unordered_map>

namespace ldzip {

    // Chunked writer for zstd-compressed column data
    class ChunkedWriter {
        private:
            std::ofstream chunks_file_;
            std::ofstream index_file_;
            std::vector<uint8_t> buffer_;
            std::vector<uint8_t> compressed_buffer_;  // Reusable compression buffer
            size_t chunk_size_;          // Number of columns per chunk
            size_t current_chunk_cols_;  // Columns buffered in current chunk
            uint64_t current_offset_;    // Current byte offset in chunks file
            int compression_level_;
            bool closed_ = false;

        public:
            ChunkedWriter(const std::string& chunks_path,
                         const std::string& index_path,
                         size_t chunk_size,
                         int compression_level = 6)
                : chunk_size_(chunk_size),
                  current_chunk_cols_(0),
                  current_offset_(0),
                  compression_level_(compression_level) {

                chunks_file_.open(chunks_path, std::ios::out | std::ios::binary | std::ios::trunc);
                if (!chunks_file_) {
                    throw std::runtime_error("Cannot open chunks file for write: " + chunks_path);
                }

                index_file_.open(index_path, std::ios::out | std::ios::binary | std::ios::trunc);
                if (!index_file_) {
                    throw std::runtime_error("Cannot open index file for write: " + index_path);
                }

                // Write initial offset (0)
                index_file_.write(reinterpret_cast<const char*>(&current_offset_), sizeof(uint64_t));
            }

            // Append raw bytes for one column to the buffer
            void appendColumn(const void* data, size_t size) {
                const uint8_t* bytes = static_cast<const uint8_t*>(data);
                buffer_.insert(buffer_.end(), bytes, bytes + size);
                current_chunk_cols_++;

                // Flush when chunk is full
                if (current_chunk_cols_ >= chunk_size_) {
                    flush();
                }
            }

            // Compress and write current buffer
            void flush() {
                if (buffer_.empty()) return;

                // Compress buffer with zstd (reuse compression buffer)
                size_t compressed_bound = ZSTD_compressBound(buffer_.size());
                if (compressed_buffer_.size() < compressed_bound) {
                    compressed_buffer_.resize(compressed_bound);
                }

                size_t compressed_size = ZSTD_compress(
                    compressed_buffer_.data(), compressed_buffer_.size(),
                    buffer_.data(), buffer_.size(),
                    compression_level_
                );

                if (ZSTD_isError(compressed_size)) {
                    throw std::runtime_error("ZSTD compression failed: " +
                                           std::string(ZSTD_getErrorName(compressed_size)));
                }

                // Write compressed data
                chunks_file_.write(reinterpret_cast<const char*>(compressed_buffer_.data()), compressed_size);
                current_offset_ += compressed_size;

                // Write new offset to index
                index_file_.write(reinterpret_cast<const char*>(&current_offset_), sizeof(uint64_t));

                // Clear buffer
                buffer_.clear();
                current_chunk_cols_ = 0;
            }

            void close() {
                if (closed_) return;
                flush(); // Write any remaining data
                if (chunks_file_.is_open()) chunks_file_.close();
                if (index_file_.is_open()) index_file_.close();
                closed_ = true;
            }

            ~ChunkedWriter() {
                close();
            }
    };

    // Chunked reader with LRU cache for decompressed chunks
    class ChunkedReader {
        private:
            std::string chunks_path_;
            std::string index_path_;
            mutable std::ifstream chunks_file_;
            mutable std::ifstream index_file_;
            size_t chunk_size_;
            size_t num_chunks_;

            // Reusable buffers to avoid reallocation on every read
            mutable std::vector<uint8_t> compressed_buffer_;
            mutable std::vector<uint8_t> decompressed_buffer_;

            // Cache: chunk_id -> decompressed data
            struct CacheEntry {
                std::vector<uint8_t> data;
                size_t access_count;
            };
            mutable std::unordered_map<size_t, CacheEntry> cache_;
            mutable size_t cache_access_counter_;
            static constexpr size_t MAX_CACHE_SIZE = 10; // Cache up to 10 chunks

        public:
            ChunkedReader(const std::string& chunks_path,
                         const std::string& index_path,
                         size_t chunk_size)
                : chunks_path_(chunks_path),
                  index_path_(index_path),
                  chunk_size_(chunk_size),
                  cache_access_counter_(0) {

                // Open index file and count chunks
                index_file_.open(index_path_, std::ios::binary | std::ios::ate);
                if (!index_file_) {
                    throw std::runtime_error("Cannot open index file: " + index_path_);
                }
                size_t index_size = index_file_.tellg();
                num_chunks_ = (index_size / sizeof(uint64_t)) - 1;

                chunks_file_.open(chunks_path_, std::ios::in | std::ios::binary);
                if (!chunks_file_) {
                    throw std::runtime_error("Cannot open chunks file: " + chunks_path_);
                }
            }

            // Read and decompress a specific chunk
            const std::vector<uint8_t>& readChunk(size_t chunk_id) const {
                // Check cache first
                auto it = cache_.find(chunk_id);
                if (it != cache_.end()) {
                    it->second.access_count = ++cache_access_counter_;
                    return it->second.data;
                }

                // Read chunk offsets from index file
                index_file_.seekg(chunk_id * sizeof(uint64_t));
                uint64_t start_offset, end_offset;
                index_file_.read(reinterpret_cast<char*>(&start_offset), sizeof(uint64_t));
                index_file_.read(reinterpret_cast<char*>(&end_offset), sizeof(uint64_t));

                size_t compressed_size = end_offset - start_offset;

                // Read compressed data (reuse buffer to avoid reallocation)
                if (compressed_buffer_.size() < compressed_size) {
                    compressed_buffer_.resize(compressed_size);
                }
                chunks_file_.clear();
                chunks_file_.seekg(start_offset);
                chunks_file_.read(reinterpret_cast<char*>(compressed_buffer_.data()), compressed_size);

                // Decompress
                unsigned long long decompressed_size = ZSTD_getFrameContentSize(
                    compressed_buffer_.data(), compressed_size
                );

                if (decompressed_size == ZSTD_CONTENTSIZE_ERROR ||
                    decompressed_size == ZSTD_CONTENTSIZE_UNKNOWN) {
                    throw std::runtime_error("Cannot determine decompressed size");
                }

                // Reuse decompression buffer to avoid reallocation
                if (decompressed_buffer_.size() < decompressed_size) {
                    decompressed_buffer_.resize(decompressed_size);
                }
                size_t result = ZSTD_decompress(
                    decompressed_buffer_.data(), decompressed_size,
                    compressed_buffer_.data(), compressed_size
                );

                if (ZSTD_isError(result)) {
                    throw std::runtime_error("ZSTD decompression failed: " +
                                           std::string(ZSTD_getErrorName(result)));
                }

                // Add to cache (evict oldest if full)
                if (cache_.size() >= MAX_CACHE_SIZE) {
                    // Evict LRU entry
                    auto oldest = cache_.begin();
                    for (auto it = cache_.begin(); it != cache_.end(); ++it) {
                        if (it->second.access_count < oldest->second.access_count) {
                            oldest = it;
                        }
                    }
                    cache_.erase(oldest);
                }

                // Copy only the used portion of the buffer into cache
                std::vector<uint8_t> cached_data(decompressed_buffer_.begin(),
                                                  decompressed_buffer_.begin() + decompressed_size);
                cache_[chunk_id] = {std::move(cached_data), ++cache_access_counter_};
                return cache_[chunk_id].data;
            }

            size_t getChunkSize() const { return chunk_size_; }
            size_t getNumChunks() const { return num_chunks_; }

            ~ChunkedReader() {
                if (chunks_file_.is_open()) chunks_file_.close();
            }
    };

} // namespace ldzip
