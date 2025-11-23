#pragma once

#include "types.hpp"
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <map>

namespace agp {

/**
 * Sequence record from FASTA/FASTQ file
 */
struct SequenceRecord {
    std::string id;
    std::string description;
    std::string sequence;
    std::string quality;  // Only for FASTQ
};

/**
 * FASTA/FASTQ file reader
 *
 * Supports:
 * - Uncompressed and gzip-compressed files
 * - FASTA and FASTQ formats
 * - Memory-mapped I/O for large files
 * - Iterator-based and callback-based processing
 */
class SequenceReader {
public:
    /**
     * Open a sequence file (auto-detects format)
     */
    explicit SequenceReader(const std::string& filename);
    ~SequenceReader();

    /**
     * Read next sequence
     * Returns false when end of file is reached
     */
    bool read_next(SequenceRecord& record);

    /**
     * Process all sequences with a callback
     */
    void for_each(std::function<void(const SequenceRecord&)> callback);

    /**
     * Read all sequences into memory
     */
    std::vector<SequenceRecord> read_all();

    /**
     * Check if file is open and valid
     */
    bool is_open() const;

    /**
     * Get file format
     */
    enum class Format { FASTA, FASTQ, UNKNOWN };
    Format get_format() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/**
 * FASTA/FASTQ file writer
 */
class SequenceWriter {
public:
    /**
     * Create a sequence file writer
     */
    explicit SequenceWriter(const std::string& filename, bool compress = false);
    ~SequenceWriter();

    /**
     * Write a sequence record
     */
    void write(const SequenceRecord& record);

    /**
     * Write multiple records
     */
    void write_all(const std::vector<SequenceRecord>& records);

    /**
     * Flush and close file
     */
    void close();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/**
 * Gene output writer (GFF3 format)
 */
class GeneWriter {
public:
    /**
     * Create a GFF3 writer
     */
    explicit GeneWriter(const std::string& filename);
    ~GeneWriter();

    /**
     * Write genes for a sequence
     */
    void write_genes(const std::string& sequence_id,
                    const std::vector<Gene>& genes,
                    const std::string& source = "AncientGenePredictor");

    /**
     * Write a single gene
     */
    void write_gene(const std::string& sequence_id,
                   const Gene& gene,
                   int gene_num,
                   const std::string& source = "AncientGenePredictor");

    /**
     * Close file
     */
    void close();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/**
 * FASTA sequence writer
 * Outputs predicted genes as nucleotide or amino acid sequences
 */
class FastaWriter {
public:
    explicit FastaWriter(const std::string& filename);
    ~FastaWriter();

    /**
     * Write a single sequence to FASTA
     */
    void write_sequence(const std::string& id,
                       const std::string& description,
                       const std::string& sequence);

    /**
     * Write genes as nucleotide FASTA
     */
    void write_genes_nucleotide(const std::string& sequence_id,
                                const std::vector<Gene>& genes);

    /**
     * Write genes as amino acid FASTA
     */
    void write_genes_protein(const std::string& sequence_id,
                            const std::vector<Gene>& genes);

    void close();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/**
 * Statistics writer for damage model and codon usage
 */
class StatsWriter {
public:
    explicit StatsWriter(const std::string& filename);
    ~StatsWriter();

    /**
     * Write codon frequency statistics
     */
    void write_codon_frequencies(const std::map<std::string, float>& frequencies);

    /**
     * Write damage model parameters
     */
    void write_damage_parameters(float lambda_5p, float lambda_3p,
                                 float delta_max, float delta_bg);

    /**
     * Write damage profile (position-specific rates)
     */
    void write_damage_profile(const std::vector<float>& ct_profile,
                             const std::vector<float>& ga_profile);

    /**
     * Write organism statistics (GC content, ENC, etc.)
     */
    void write_organism_stats(float gc_content, float enc,
                             const std::map<std::string, float>& rscu);

    void close();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/**
 * Sequence utilities
 */
class SequenceUtils {
public:
    /**
     * Reverse complement a DNA sequence
     */
    static std::string reverse_complement(const std::string& seq);

    /**
     * Translate DNA to protein (single frame)
     */
    static std::string translate(const std::string& seq, int frame = 0);

    /**
     * Find all ORFs in a sequence (both strands)
     */
    static std::vector<Gene> find_orfs(const std::string& seq,
                                      size_t min_length = 30);

    /**
     * Calculate GC content
     */
    static float gc_content(const std::string& seq);

    /**
     * Validate DNA sequence (only ACGT or N)
     */
    static bool is_valid_dna(const std::string& seq);

    /**
     * Clean sequence (remove whitespace, convert to uppercase)
     */
    static std::string clean(const std::string& seq);

    /**
     * Convert to encoded nucleotide vector
     */
    static std::vector<Nucleotide> encode(const std::string& seq);

    /**
     * Convert from encoded nucleotide vector
     */
    static std::string decode(const std::vector<Nucleotide>& seq);

    /**
     * Calculate sequence composition
     */
    static std::array<size_t, 4> composition(const std::string& seq);

private:
    static char complement(char nt);
};

/**
 * Batch sequence processor
 *
 * Efficiently processes large sequence files using parallel processing
 */
class BatchProcessor {
public:
    /**
     * Create batch processor
     * @param num_threads Number of worker threads (0 = auto)
     * @param batch_size Number of sequences per batch
     */
    BatchProcessor(size_t num_threads = 0, size_t batch_size = 100);
    ~BatchProcessor();

    /**
     * Process sequences from a file
     */
    void process_file(
        const std::string& input_file,
        const std::string& output_file,
        std::function<std::vector<Gene>(const SequenceRecord&)> processor);

    /**
     * Set progress callback
     */
    void set_progress_callback(
        std::function<void(size_t processed, size_t total)> callback);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace agp
