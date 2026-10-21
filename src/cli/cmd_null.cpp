// dart null: fit the decoy null sidecar (.null) from decoy hits, or pool existing sidecars.
//
// hits2emi --decoy-hits does this per run. This command is for a null shared across runs
// of one database and search configuration: fit once from a few samples' decoy hits (or
// pool their sidecars), then pass the file to every damage-annotate --em-null-model.

#include "subcommand.hpp"
#include "dart/version.h"
#include "dart/null_model.hpp"

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace dart {
namespace cli {

static void null_usage() {
    std::cout << "Usage: dart null [-o FILE] [--db-label STR] INPUT...\n\n"
              << "Fit the decoy null model (Gumbel per aligned-length bin) and write a .null sidecar.\n\n"
              << "INPUT is a decoy hits TSV (16-column MMseqs2 output, plain or .gz) or a .null sidecar\n"
              << "written by hits2emi --decoy-hits or this command; any mix, all pooled into one fit.\n"
              << "In a hits TSV, rows whose query id starts with DECOY_ are the decoys; a file with no\n"
              << "such row is taken as all-decoy, and the count taken is printed either way.\n\n"
              << "Options:\n"
              << "  -o, --output FILE   Sidecar to write (default: <INPUT>.null with a single input)\n"
              << "  --db-label STR      Database label recorded in the sidecar; damage-annotate warns\n"
              << "                      when its own --db-label differs\n"
              << "  -h, --help          Show this help\n";
}

int cmd_null(int argc, char* argv[]) {
    std::string out, db_label;
    std::vector<std::string> inputs;
    for (int i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            out = argv[++i];
        } else if (strcmp(argv[i], "--db-label") == 0 && i + 1 < argc) {
            db_label = argv[++i];
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            null_usage();
            return 0;
        } else if (argv[i][0] == '-') {
            std::cerr << "Unknown option: " << argv[i] << "\nRun 'dart null --help' for usage.\n";
            return 1;
        } else {
            inputs.emplace_back(argv[i]);
        }
    }
    if (inputs.empty()) {
        std::cerr << "Error: at least one decoy hits TSV or .null sidecar is required\n";
        null_usage();
        return 1;
    }
    if (out.empty()) {
        if (inputs.size() != 1) {
            std::cerr << "Error: -o FILE is required with more than one input\n";
            return 1;
        }
        out = inputs[0] + ".null";
    }

    DecoyStats stats;
    NullProvenance prov;
    prov.dart_version = DART_VERSION;
    prov.db_label = db_label;
    try {
        for (const auto& in : inputs) {
            if (is_null_sidecar(in)) {
                NullModel m; NullProvenance p; std::string err;
                if (!read_null_model(in, m, &p, &err)) { std::cerr << "Error: " << err << "\n"; return 1; }
                if (p.legacy_tsv || m.sum.empty()) {
                    std::cerr << "Error: " << in << " carries no per-bin moments (deprecated TSV sidecar); refit it from its decoy hits\n";
                    return 1;
                }
                stats.add(m);
                stats.n_hits += p.n_decoy_hits;
                prov.n_decoy_queries += p.n_decoy_queries;
                prov.n_decoy_hits += p.n_decoy_hits;
                prov.inputs.insert(prov.inputs.end(), p.inputs.begin(), p.inputs.end());
                std::cerr << in << ": pooled " << p.n_decoy_queries << " decoy queries from sidecar\n";
            } else {
                uint64_t skipped = 0;
                bool by_prefix = false;
                const uint64_t taken = accumulate_decoy_hits(in, stats, prov, &skipped, &by_prefix);
                std::cerr << in << ": " << taken << " decoy queries taken"
                          << (by_prefix ? " (DECOY_ rows" + (skipped ? "; " + std::to_string(skipped) + " queries without the prefix ignored" : std::string()) + ")"
                                        : " (no DECOY_ prefix in the file: every query taken as a decoy)")
                          << "\n";
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    if (stats.n_queries == 0) {
        std::cerr << "Error: no decoy queries in the inputs\n";
        return 1;
    }
    const NullModel m = fit_null_model(stats);
    std::string err;
    if (!write_null_model(out, m, prov, &err)) { std::cerr << "Error: " << err << "\n"; return 1; }
    NullProvenance written; NullModel check;
    if (!read_null_model(out, check, &written, &err)) { std::cerr << "Error: " << err << "\n"; return 1; }
    print_null_model(std::cout, check, written);
    std::cerr << "Wrote " << out << "\n";
    return 0;
}

namespace {
    struct NullRegistrar {
        NullRegistrar() {
            SubcommandRegistry::instance().register_command(
                "null",
                "Fit or pool the decoy null model sidecar (.null) for damage-annotate",
                cmd_null, 32);
        }
    } null_registrar;
}

}  // namespace cli
}  // namespace dart
