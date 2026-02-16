// Main entry point for agp CLI with subcommand dispatch
//
// Usage:
//   agp predict <input.fq> [options]     Gene prediction (main function)
//   agp sample-damage <input.fq>         Quick sample-level damage check
//   agp damage-annotate --emi <hits.emi> Post-mapping damage annotation
//   agp validate <pred.gff> <ref.gff>    Validate predictions

#include "subcommand.hpp"
#include "agp/version.h"
#include <iostream>
#include <cstring>

int main(int argc, char* argv[]) {
    auto& registry = agp::cli::SubcommandRegistry::instance();

    // Handle no arguments
    if (argc < 2) {
        registry.print_help(argv[0]);
        return 1;
    }

    const char* first_arg = argv[1];

    // Handle --help and --version at top level
    if (strcmp(first_arg, "--help") == 0 || strcmp(first_arg, "-h") == 0) {
        registry.print_help(argv[0]);
        return 0;
    }

    if (strcmp(first_arg, "--version") == 0 || strcmp(first_arg, "-V") == 0) {
        std::cout << "agp " << AGP_VERSION << "\n";
        return 0;
    }

    // Dispatch to subcommand
    if (registry.has_command(first_arg)) {
        return registry.run_command(first_arg, argc - 1, argv + 1);
    }

    // Unknown command
    std::cerr << "Unknown command: " << first_arg << "\n";
    std::cerr << "Run 'agp --help' for usage information.\n";
    return 1;
}
