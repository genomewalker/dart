// dart profile: reference-free ancient DNA damage profiling.
//
// Argument parsing and FASTQ input for the libtaph library profiler
// (taph::run_library_profile). `fqdup profile` runs the same profiler, so the JSON,
// .tdmg and TSV outputs are byte-identical for the same reads and options.

#include "subcommand.hpp"
#include "dart/sequence_io.hpp"
#include "dart/version.h"
#include <taph/library_profile.hpp>

#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>

namespace dart {
namespace cli {

namespace {

class SequenceSeqReader : public taph::SeqReader {
public:
    // ISA-L like fqdup: rapidgzip's parallel chunk cache cost ~0.7 GB here for no speed gain.
    explicit SequenceSeqReader(const std::string& path) : reader_(path, true) {}
    bool next(std::string& seq) override {
        if (!reader_.read_next(rec_)) return false;
        seq = std::move(rec_.sequence);
        return true;
    }
private:
    SequenceReader reader_;
    SequenceRecord rec_;
};

void print_profile_usage() {
    std::cout <<
        "Usage: dart profile (-i FILE | FILE | -1 R1.fq -2 R2.fq) [options]\n\n"
        "Reference-free ancient DNA damage profiling. Same profiler and outputs as\n"
        "`fqdup profile`; the .tdmg is what `dart predict --library-profile` reads.\n\n"
        "Input:\n"
        "  -i, --input FILE          Single-end or merged FASTQ (.gz or plain)\n"
        "  -1 FILE -2 FILE           Paired raw reads (R1 gives 5' counts, R2 3' counts)\n"
        "Output:\n"
        "  -o, --json FILE           Damage profile JSON\n"
        "  --json-level LEVEL        summary | standard | full (default: standard)\n"
        "  --damage-summary FILE     Library profile sidecar (.tdmg) for predict --library-profile\n"
        "  --tsv FILE                Per-position table\n"
        "  --count-table-json FILE   Stop-channel count tables (JSONL)\n"
        "  --html FILE               Interactive report\n"
        "Options:\n"
        "  -t, --threads N           Worker threads (default: all cores)\n"
        "  --library-type TYPE       auto | ds | ss (default: auto)\n"
        "  --length-bins SPEC        auto | off | N | e1,e2,... (default: auto)\n"
        "  --min-length N            Below 30, also scores 16-29 bp reads (floor max(16, N))\n"
        "  --max-length N            Ignore reads longer than N (default: off)\n"
        "  --adapter-scan-reads N    Reads scanned for adapter stubs (default: 1000000; 0 = all)\n"
        "  --mask-threshold FLOAT    Damage mask cutoff (default: 0.05)\n"
        "  --subst-in FILE[,FILE]    fqdup merge .bsubst files\n"
        "  --no-oxog                 Skip the 8-oxoG channels\n";
}

}  // namespace

int cmd_profile(int argc, char* argv[]) {
    taph::LibraryProfileOptions o;
    o.n_threads = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        const bool has_value = i + 1 < argc;
        if ((arg == "-i" || arg == "--input") && has_value) {
            o.in_path = argv[++i];
        } else if ((arg == "-1" || arg == "--read1") && has_value) {
            o.r1_path = argv[++i];
        } else if ((arg == "-2" || arg == "--read2") && has_value) {
            o.r2_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output" || arg == "--json") && has_value) {
            o.json_path = argv[++i];
        } else if (arg == "--json-level" && has_value) {
            const std::string lvl(argv[++i]);
            if (lvl == "summary")       o.json_level = 1;
            else if (lvl == "standard") o.json_level = 2;
            else if (lvl == "full")     o.json_level = 3;
            else { std::cerr << "Error: --json-level must be summary, standard or full\n"; return 1; }
        } else if (arg == "--damage-summary" && has_value) {
            o.damage_summary_path = argv[++i];
        } else if (arg == "--tsv" && has_value) {
            o.tsv_path = argv[++i];
        } else if (arg == "--count-table-json" && has_value) {
            o.count_table_path = argv[++i];
        } else if (arg == "--html" && has_value) {
            o.html_path = argv[++i];
        } else if ((arg == "-t" || arg == "--threads" || arg == "-p") && has_value) {
            o.n_threads = std::max(1, std::stoi(argv[++i]));
        } else if (arg == "--library-type" && has_value) {
            const std::string lt(argv[++i]);
            if (lt == "ss" || lt == "single-stranded")
                o.forced_lib = taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            else if (lt == "ds" || lt == "double-stranded")
                o.forced_lib = taph::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            else if (lt != "auto") { std::cerr << "Error: --library-type must be auto, ds or ss\n"; return 1; }
        } else if (arg == "--length-bins" && has_value) {
            bool ok = true;
            o.lb_opts = taph::parse_length_bins(argv[++i], ok);
            if (!ok) { std::cerr << "Error: invalid --length-bins: " << argv[i] << "\n"; return 1; }
        } else if (arg == "--min-length" && has_value) {
            o.min_length = std::stoi(argv[++i]);
        } else if (arg == "--max-length" && has_value) {
            o.max_length = std::stoi(argv[++i]);
        } else if (arg == "--adapter-scan-reads" && has_value) {
            const long long v = std::stoll(argv[++i]);
            if (v < 0) { std::cerr << "Error: --adapter-scan-reads must be >= 0\n"; return 1; }
            o.adapter_scan_reads = v;
        } else if (arg == "--mask-threshold" && has_value) {
            o.mask_threshold = std::stod(argv[++i]);
        } else if (arg == "--subst-in" && has_value) {
            std::stringstream ss(argv[++i]);
            for (std::string p; std::getline(ss, p, ',');) if (!p.empty()) o.subst_in_paths.push_back(p);
        } else if (arg == "--no-oxog") {
            o.run_oxog = false;
        } else if (arg == "-h" || arg == "--help") {
            print_profile_usage();
            return 0;
        } else if (!arg.empty() && arg[0] != '-' && o.in_path.empty()) {
            o.in_path = arg;
        } else {
            std::cerr << "Error: unknown argument: " << arg << "\n";
            print_profile_usage();
            return 1;
        }
    }
    if (o.r1_path.empty() != o.r2_path.empty()) {
        std::cerr << "Error: -1 and -2 must both be given for paired input\n";
        return 1;
    }
    if (o.in_path.empty() && o.r1_path.empty()) {
        std::cerr << "Error: give an input FASTQ (-i FILE) or paired reads (-1/-2)\n";
        print_profile_usage();
        return 1;
    }
    o.producer = std::string("dart/") + DART_VERSION;
    o.log_info = [](const std::string& m) { std::cerr << m << "\n"; };
    o.log_warn = [](const std::string& m) { std::cerr << "Warning: " << m << "\n"; };

    auto open_reader = [](const std::string& path, size_t) -> std::unique_ptr<taph::SeqReader> {
        return std::make_unique<SequenceSeqReader>(path);
    };
    return taph::run_library_profile(o, open_reader);
}

namespace {
struct ProfileRegistrar {
    ProfileRegistrar() {
        SubcommandRegistry::instance().register_command(
            "profile", "Reference-free ancient DNA damage profiling", cmd_profile, 10);
        SubcommandRegistry::instance().register_command(
            "damage", "Alias for 'profile'", cmd_profile, 11);
    }
} profile_registrar;
}  // namespace

}  // namespace cli
}  // namespace dart
