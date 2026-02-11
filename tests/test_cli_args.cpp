// Unit tests for CLI argument parsing
// Compile: g++ -std=c++20 -I../include -I../src -o test_cli_args test_cli_args.cpp ../src/cli/args.cpp

#include "cli/args.hpp"
#include <iostream>
#include <cstring>
#include <cassert>
#include <vector>

// Helper to create argv from strings
class ArgvBuilder {
public:
    ArgvBuilder& add(const char* arg) {
        args.push_back(strdup(arg));
        return *this;
    }

    int argc() const { return static_cast<int>(args.size()); }
    char** argv() { return args.data(); }

    ~ArgvBuilder() {
        for (auto* arg : args) free(arg);
    }

private:
    std::vector<char*> args;
};

// Test basic argument parsing
void test_basic_args() {
    std::cout << "Testing basic arguments... ";

    ArgvBuilder builder;
    builder.add("agp")
           .add("-i").add("input.fa")
           .add("-o").add("output.gff");

    auto opts = agp::cli::parse_args(builder.argc(), builder.argv());

    assert(opts.input_file == "input.fa");
    assert(opts.output_file == "output.gff");
    assert(opts.use_damage == true);  // Default
    assert(opts.verbose == false);    // Default

    std::cout << "PASSED\n";
}

// Test default values
void test_defaults() {
    std::cout << "Testing default values... ";

    ArgvBuilder builder;
    builder.add("agp").add("-i").add("test.fa");

    auto opts = agp::cli::parse_args(builder.argc(), builder.argv());

    assert(opts.output_file == "predictions.gff");
    assert(opts.min_length == 30);
    assert(opts.min_coding_prob == 0.3f);
    assert(opts.num_threads == 0);
    assert(opts.dual_strand == true);
    assert(opts.best_strand_only == true);
    assert(opts.aggregate_damage == true);
    assert(opts.domain_name == "gtdb");
    assert(opts.forced_library_type == agp::cli::LibraryType::UNKNOWN);

    std::cout << "PASSED\n";
}

// Test boolean flags
void test_boolean_flags() {
    std::cout << "Testing boolean flags... ";

    ArgvBuilder builder;
    builder.add("agp")
           .add("-i").add("test.fa")
           .add("--no-damage")
           .add("--verbose")
           .add("--both-strands")
           .add("--iterative-damage");

    auto opts = agp::cli::parse_args(builder.argc(), builder.argv());

    assert(opts.use_damage == false);
    assert(opts.verbose == true);
    assert(opts.best_strand_only == false);
    assert(opts.iterative_damage == true);

    std::cout << "PASSED\n";
}

// Test numeric arguments
void test_numeric_args() {
    std::cout << "Testing numeric arguments... ";

    ArgvBuilder builder;
    builder.add("agp")
           .add("-i").add("test.fa")
           .add("--min-length").add("50")
           .add("--min-coding-prob").add("0.5")
           .add("--threads").add("8");

    auto opts = agp::cli::parse_args(builder.argc(), builder.argv());

    assert(opts.min_length == 50);
    assert(opts.min_coding_prob == 0.5f);
    assert(opts.num_threads == 8);

    std::cout << "PASSED\n";
}

// Test library type parsing
void test_library_type() {
    std::cout << "Testing library type... ";

    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--library-type").add("ds");
        auto opts = agp::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.forced_library_type == agp::cli::LibraryType::DOUBLE_STRANDED);
    }

    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--library-type").add("ss");
        auto opts = agp::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.forced_library_type == agp::cli::LibraryType::SINGLE_STRANDED);
    }

    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--library-type").add("auto");
        auto opts = agp::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.forced_library_type == agp::cli::LibraryType::UNKNOWN);
    }

    std::cout << "PASSED\n";
}

// Test domain and metagenome mode
void test_domain() {
    std::cout << "Testing domain options... ";

    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--domain").add("fungi");
        auto opts = agp::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.domain_name == "fungi");
        assert(opts.metagenome_mode == false);
    }

    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--domain").add("meta");
        auto opts = agp::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.domain_name == "meta");
        assert(opts.metagenome_mode == true);
    }

    std::cout << "PASSED\n";
}

// Test adaptive ORF mode
void test_adaptive_orf() {
    std::cout << "Testing adaptive ORF options... ";

    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--adaptive-orf");
        auto opts = agp::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.orf_confidence_threshold == 0.3f);
        assert(opts.best_strand_only == false);
    }

    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--orf-confidence").add("0.5");
        auto opts = agp::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.orf_confidence_threshold == 0.5f);
        assert(opts.best_strand_only == false);
    }

    std::cout << "PASSED\n";
}

// Test strand weight clamping
void test_strand_weight() {
    std::cout << "Testing strand weight clamping... ";

    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--strand-weight").add("1.5");
        auto opts = agp::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.strand_weight == 1.0f);  // Clamped to max
    }

    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--strand-weight").add("-0.5");
        auto opts = agp::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.strand_weight == 0.0f);  // Clamped to min
    }

    std::cout << "PASSED\n";
}

// Test all-frames mode
void test_all_frames() {
    std::cout << "Testing all-frames mode... ";

    ArgvBuilder builder;
    builder.add("agp").add("-i").add("test.fa").add("--all-frames");
    auto opts = agp::cli::parse_args(builder.argc(), builder.argv());

    assert(opts.all_frames == true);
    assert(opts.best_strand_only == false);  // All-frames implies both strands

    std::cout << "PASSED\n";
}

int main() {
    std::cout << "\n=== CLI Argument Parsing Tests ===\n\n";

    test_basic_args();
    test_defaults();
    test_boolean_flags();
    test_numeric_args();
    test_library_type();
    test_domain();
    test_adaptive_orf();
    test_strand_weight();
    test_all_frames();

    std::cout << "\nAll tests passed!\n";
    return 0;
}
