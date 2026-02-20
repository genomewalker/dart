// Main entry point for dart CLI with subcommand dispatch
//
// Usage:
//   dart predict <input.fq> [options]     Gene prediction (main function)
//   dart sample-damage <input.fq>         Quick sample-level damage check
//   dart damage-annotate --emi <hits.emi> Post-mapping damage annotation
//   dart validate <pred.gff> <ref.gff>    Validate predictions

#include "subcommand.hpp"
#include "dart/version.h"
#include <iostream>
#include <cstring>

int main(int argc, char* argv[]) {
    auto& registry = dart::cli::SubcommandRegistry::instance();

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
        std::cout << "dart " << DART_VERSION << "\n";
        return 0;
    }

    // Dispatch to subcommand
    if (registry.has_command(first_arg)) {
        return registry.run_command(first_arg, argc - 1, argv + 1);
    }

    // Unknown command
    std::cerr << "Unknown command: " << first_arg << "\n";
    std::cerr << "Run 'dart --help' for usage information.\n";
    return 1;
}
