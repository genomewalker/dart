// Native bounded random CLI parser target: default/override precedence,
// byte-preserving paths and rejection of truncated option values.
#include "cli/args.hpp"
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <random>
#include <string>

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

// Fixed seed makes every parser regression reproducible. Generated filenames are
// arbitrary non-NUL bytes; they must round-trip and cannot affect default selection.
static void fuzz_decoy_defaults() {
    std::mt19937 rng(0xdec010);
    for (size_t trial = 0; trial < 4096; ++trial) {
        ArgvBuilder args;
        args.add("dart").add("-i").add("reads.fa");
        std::string aa, masked;
        bool explicit_interval = false;
        int interval = 0;
        for (size_t step = 0, n = rng() % 12; step < n; ++step) {
            const unsigned op = rng() % 3;
            if (op < 2) {
                std::string path;
                for (size_t j = 0, len = 1 + rng() % 64; j < len; ++j)
                    path.push_back(static_cast<char>(1 + rng() % 255));
                args.add(op == 0 ? "--fasta-aa" : "--fasta-aa-masked").add(path.c_str());
                (op == 0 ? aa : masked) = path;
            } else {
                interval = rng() % 100;
                explicit_interval = true;
                args.add("--decoy-every").add(std::to_string(interval).c_str());
            }
        }
        const auto parsed = dart::cli::parse_args(args.argc(), args.argv());
        assert(parsed.fasta_aa == aa && parsed.fasta_aa_masked == masked);
        assert(parsed.decoy_every == (explicit_interval ? interval : aa.empty() ? 0 : 10));
        // Every value-taking option must reject a truncated argument vector.
        args.add(trial % 2 ? "--decoy-every" : "--fasta-aa");
        try {
            dart::cli::parse_args(args.argc(), args.argv());
            assert(false && "truncated option accepted");
        } catch (const dart::cli::ParseArgsExit& ex) { assert(ex.exit_code() == 1); }
    }
    std::cout << "Decoy defaults fuzz: 4096 cases, up to 11 operations, paths <= 64 bytes PASSED\n";
}

int main() { fuzz_decoy_defaults(); }
