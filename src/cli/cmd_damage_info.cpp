/**
 * @file cmd_damage_info.cpp
 * @brief Display damage profile information from AGD index files.
 */

#include "subcommand.hpp"
#include "dart/damage_index.hpp"
#include "dart/damage_index_reader.hpp"
#include "dart/version.h"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " damage-info <file.agd> [options]\n\n"
              << "Display damage profile information from AGD index file.\n\n"
              << "Options:\n"
              << "  --json           Output as JSON (default: human-readable)\n"
              << "  -o, --output F   Write to file instead of stdout\n"
              << "  -h, --help       Show this help message\n";
}

std::string library_type_str(uint8_t lt) {
    switch (lt) {
        case 1: return "single-stranded";
        case 2: return "double-stranded";
        default: return "unknown";
    }
}

void print_human_readable(const dart::DamageIndexReader& reader, std::ostream& os) {
    os << "=== DART Damage Index (AGD) ===\n\n";
    os << "Records: " << reader.record_count() << "\n\n";

    os << "--- Core Damage Metrics ---\n";
    os << std::fixed << std::setprecision(2);
    os << "d_max:                  " << (reader.d_max() * 100.0f) << "%\n";
    os << std::setprecision(4);
    os << "lambda (decay rate):    " << reader.lambda() << "\n";
    os << "Library type:           " << library_type_str(reader.library_type()) << "\n\n";

    os << "--- Two-Channel Validation ---\n";
    os << "damage_validated:       " << (reader.damage_validated() ? "yes" : "no") << "\n";
    os << "damage_artifact:        " << (reader.damage_artifact() ? "yes" : "no") << "\n";
    os << "channel_b_valid:        " << (reader.channel_b_valid() ? "yes" : "no") << "\n";
    os << std::setprecision(2);
    os << "stop_decay_llr:         " << reader.stop_decay_llr() << "\n";
    os << std::setprecision(4);
    os << "terminal_shift:         " << reader.terminal_shift() << "\n";
    os << "damage_informative:     " << (reader.damage_informative() ? "yes" : "no") << "\n\n";

    os << "--- Channel C: Oxidation (8-oxoG, G→T) ---\n";
    os << "ox_rate_terminal:       " << reader.ox_rate_terminal() << "\n";
    os << "ox_rate_interior:       " << reader.ox_rate_interior() << "\n";
    os << std::setprecision(2);
    os << "ox_uniformity_ratio:    " << reader.ox_uniformity_ratio() << "\n";
    os << "ox_damage_detected:     " << (reader.ox_damage_detected() ? "yes" : "no") << "\n\n";

    os << "--- Channel D: Depurination ---\n";
    os << std::setprecision(4);
    os << "purine_enrichment_5':   " << reader.purine_enrichment_5prime() << "\n";
    os << "purine_enrichment_3':   " << reader.purine_enrichment_3prime() << "\n";
    os << "depurination_detected:  " << (reader.depurination_detected() ? "yes" : "no") << "\n";
}

void print_json(const dart::DamageIndexReader& reader, std::ostream& os) {
    os << "{\n";
    os << "  \"version\": \"" << DART_VERSION << "\",\n";
    os << "  \"records\": " << reader.record_count() << ",\n";
    os << "  \"damage\": {\n";

    os << std::fixed << std::setprecision(2);
    os << "    \"d_max\": " << (reader.d_max() * 100.0f) << ",\n";
    os << std::setprecision(4);
    os << "    \"lambda\": " << reader.lambda() << ",\n";
    os << "    \"library_type\": " << static_cast<int>(reader.library_type()) << ",\n";
    os << "    \"library_type_str\": \"" << library_type_str(reader.library_type()) << "\",\n";

    os << "    \"damage_validated\": " << (reader.damage_validated() ? "true" : "false") << ",\n";
    os << "    \"damage_artifact\": " << (reader.damage_artifact() ? "true" : "false") << ",\n";
    os << "    \"channel_b_valid\": " << (reader.channel_b_valid() ? "true" : "false") << ",\n";
    os << std::setprecision(2);
    os << "    \"stop_decay_llr\": " << reader.stop_decay_llr() << ",\n";
    os << std::setprecision(4);
    os << "    \"terminal_shift\": " << reader.terminal_shift() << ",\n";
    os << "    \"damage_informative\": " << (reader.damage_informative() ? "true" : "false") << ",\n";

    os << "    \"ox_rate_terminal\": " << reader.ox_rate_terminal() << ",\n";
    os << "    \"ox_rate_interior\": " << reader.ox_rate_interior() << ",\n";
    os << std::setprecision(2);
    os << "    \"ox_uniformity_ratio\": " << reader.ox_uniformity_ratio() << ",\n";
    os << "    \"ox_damage_detected\": " << (reader.ox_damage_detected() ? "true" : "false") << ",\n";

    os << std::setprecision(4);
    os << "    \"purine_enrichment_5prime\": " << reader.purine_enrichment_5prime() << ",\n";
    os << "    \"purine_enrichment_3prime\": " << reader.purine_enrichment_3prime() << ",\n";
    os << "    \"depurination_detected\": " << (reader.depurination_detected() ? "true" : "false") << "\n";

    os << "  }\n";
    os << "}\n";
}

}  // namespace

int cmd_damage_info(int argc, char* argv[]) {
    std::string agd_file;
    std::string output_file;
    bool json_output = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--json") {
            json_output = true;
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg[0] != '-' && agd_file.empty()) {
            agd_file = arg;
        }
    }

    if (agd_file.empty()) {
        std::cerr << "Error: AGD file required\n\n";
        print_usage(argv[0]);
        return 1;
    }

    dart::DamageIndexReader reader(agd_file, false);
    if (!reader.is_valid()) {
        std::cerr << "Error: Cannot open AGD file: " << agd_file << "\n";
        return 1;
    }

    if (!output_file.empty()) {
        std::ofstream ofs(output_file);
        if (!ofs) {
            std::cerr << "Error: Cannot open output file: " << output_file << "\n";
            return 1;
        }
        json_output ? print_json(reader, ofs) : print_human_readable(reader, ofs);
    } else {
        json_output ? print_json(reader, std::cout) : print_human_readable(reader, std::cout);
    }

    return 0;
}

namespace {
    struct DamageInfoRegistrar {
        DamageInfoRegistrar() {
            dart::cli::SubcommandRegistry::instance().register_command(
                "damage-info",
                "Display damage profile from AGD index file",
                cmd_damage_info, 45);
        }
    };
    static DamageInfoRegistrar registrar;
}
