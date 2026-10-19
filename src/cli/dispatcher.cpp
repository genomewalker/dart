// Main entry point for dart CLI with subcommand dispatch
//
// Usage:
//   dart predict <input.fq> [options]     Gene prediction (main function)
//   dart profile <input.fq>               Reference-free damage profiling
//   dart damage <input.fq>               Alias for 'profile'
//   dart damage-annotate --emi <hits.emi> Post-mapping damage annotation
//   dart validate <pred.gff> <ref.gff>    Validate predictions

#include "subcommand.hpp"
#include "dart/version.h"
#include <iostream>
#include <cstring>
#include <stdexcept>

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
        std::cout << "libtaph " << DART_LIBTAPH_SHA << "\n";
        return 0;
    }

    // Dispatch to subcommand
    // Uncaught exceptions (a file that cannot be opened, std::stoi on a non-number)
    // would otherwise abort with a core dump; report them and exit 1.
    if (registry.has_command(first_arg)) {
        try {
            return registry.run_command(first_arg, argc - 1, argv + 1);
        } catch (const std::invalid_argument& e) {
            std::cerr << "dart " << first_arg << ": invalid numeric value (" << e.what()
                      << "); check the option values.\n";
        } catch (const std::out_of_range& e) {
            std::cerr << "dart " << first_arg << ": value out of range (" << e.what() << ").\n";
        } catch (const std::exception& e) {
            std::cerr << "dart " << first_arg << ": " << e.what() << "\n";
        }
        std::cerr << "Run 'dart " << first_arg << " --help' for usage.\n";
        return 1;
    }

    // Unknown command
    std::cerr << "Unknown command: " << first_arg << "\n";
    std::cerr << "Run 'dart --help' for usage information.\n";
    return 1;
}
