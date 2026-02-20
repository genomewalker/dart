// Unit tests for CLI argument parsing
// Compile: g++ -std=c++20 -I../include -I../src -o test_cli_args test_cli_args.cpp ../src/cli/args.cpp

#include "cli/args.hpp"
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>

class ArgvBuilder {
public:
    ArgvBuilder& add(const char* arg) {
        args_.push_back(strdup(arg));
        return *this;
    }

    int argc() const { return static_cast<int>(args_.size()); }
    char** argv() { return args_.data(); }

    ~ArgvBuilder() {
        for (char* arg : args_) {
            std::free(arg);
        }
    }

private:
    std::vector<char*> args_;
};

static void expect_parse_exit(int expected_code, ArgvBuilder& builder) {
    bool threw = false;
    try {
        (void)dart::cli::parse_args(builder.argc(), builder.argv());
    } catch (const dart::cli::ParseArgsExit& e) {
        threw = true;
        assert(e.exit_code() == expected_code);
    }
    assert(threw);
}

void test_basic_args() {
    std::cout << "Testing basic args... ";
    ArgvBuilder builder;
    builder.add("agp").add("-i").add("input.fa").add("-o").add("output.gff");
    auto opts = dart::cli::parse_args(builder.argc(), builder.argv());
    assert(opts.input_file == "input.fa");
    assert(opts.output_file == "output.gff");
    std::cout << "PASSED\n";
}

void test_defaults() {
    std::cout << "Testing defaults... ";
    ArgvBuilder builder;
    builder.add("agp").add("-i").add("test.fa");
    auto opts = dart::cli::parse_args(builder.argc(), builder.argv());
    assert(opts.output_file == "predictions.gff");
    assert(opts.min_length == 30);
    assert(opts.num_threads == 0);
    assert(opts.use_damage == true);
    assert(opts.aggregate_damage == true);
    assert(opts.damage_only == false);
    assert(opts.domain_name == "gtdb");
    assert(opts.orf_min_aa == 10);
    assert(opts.adaptive_orf == false);
    assert(opts.forced_library_type == dart::cli::LibraryType::UNKNOWN);
    std::cout << "PASSED\n";
}

void test_boolean_flags() {
    std::cout << "Testing boolean flags... ";
    ArgvBuilder builder;
    builder.add("agp")
           .add("-i").add("test.fa")
           .add("--no-damage")
           .add("--no-aggregate")
           .add("--damage-only")
           .add("--adaptive")
           .add("--verbose");
    auto opts = dart::cli::parse_args(builder.argc(), builder.argv());
    assert(opts.use_damage == false);
    assert(opts.aggregate_damage == false);
    assert(opts.damage_only == true);
    assert(opts.adaptive_orf == true);
    assert(opts.verbose == true);
    std::cout << "PASSED\n";
}

void test_numeric_args() {
    std::cout << "Testing numeric args... ";
    ArgvBuilder builder;
    builder.add("agp")
           .add("-i").add("test.fa")
           .add("--min-length").add("50")
           .add("--threads").add("8")
           .add("--orf-min-aa").add("20");
    auto opts = dart::cli::parse_args(builder.argc(), builder.argv());
    assert(opts.min_length == 50);
    assert(opts.num_threads == 8);
    assert(opts.orf_min_aa == 20);
    std::cout << "PASSED\n";
}

void test_output_options() {
    std::cout << "Testing output options... ";
    ArgvBuilder builder;
    builder.add("agp")
           .add("-i").add("test.fa")
           .add("--fasta-nt").add("out_nt.fa")
           .add("--fasta-nt-corrected").add("out_nt_corr.fa")
           .add("--fasta-aa").add("out_aa.fa")
           .add("--fasta-aa-masked").add("out_aa_masked.fa")
           .add("--summary").add("summary.json")
           .add("--damage-index").add("out.agd");
    auto opts = dart::cli::parse_args(builder.argc(), builder.argv());
    assert(opts.fasta_nt == "out_nt.fa");
    assert(opts.fasta_nt_corrected == "out_nt_corr.fa");
    assert(opts.fasta_aa == "out_aa.fa");
    assert(opts.fasta_aa_masked == "out_aa_masked.fa");
    assert(opts.summary_file == "summary.json");
    assert(opts.damage_index == "out.agd");
    std::cout << "PASSED\n";
}

void test_library_type() {
    std::cout << "Testing library type... ";
    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--library-type").add("ds");
        auto opts = dart::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.forced_library_type == dart::cli::LibraryType::DOUBLE_STRANDED);
    }
    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--library-type").add("ss");
        auto opts = dart::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.forced_library_type == dart::cli::LibraryType::SINGLE_STRANDED);
    }
    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--library-type").add("auto");
        auto opts = dart::cli::parse_args(builder.argc(), builder.argv());
        assert(opts.forced_library_type == dart::cli::LibraryType::UNKNOWN);
    }
    std::cout << "PASSED\n";
}

void test_validation_errors() {
    std::cout << "Testing validation errors... ";
    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--threads").add("0");
        expect_parse_exit(1, builder);
    }
    {
        ArgvBuilder builder;
        builder.add("agp").add("-i").add("test.fa").add("--threads").add("-2");
        expect_parse_exit(1, builder);
    }
    {
        ArgvBuilder builder;
        builder.add("agp").add("--min-length").add("x");
        expect_parse_exit(1, builder);
    }
    {
        ArgvBuilder builder;
        builder.add("agp");
        expect_parse_exit(1, builder);
    }
    std::cout << "PASSED\n";
}

void test_controlled_exits() {
    std::cout << "Testing help/version controlled exits... ";
    {
        ArgvBuilder builder;
        builder.add("agp").add("--help");
        expect_parse_exit(0, builder);
    }
    {
        ArgvBuilder builder;
        builder.add("agp").add("--version");
        expect_parse_exit(0, builder);
    }
    std::cout << "PASSED\n";
}

int main() {
    std::cout << "\n=== CLI Argument Parsing Tests ===\n\n";
    test_basic_args();
    test_defaults();
    test_boolean_flags();
    test_numeric_args();
    test_output_options();
    test_library_type();
    test_validation_errors();
    test_controlled_exits();
    std::cout << "\nAll tests passed!\n";
    return 0;
}
