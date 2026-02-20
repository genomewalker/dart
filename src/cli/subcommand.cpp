#include "subcommand.hpp"
#include "dart/version.h"
#include <iostream>
#include <algorithm>

namespace dart {
namespace cli {

SubcommandRegistry& SubcommandRegistry::instance() {
    static SubcommandRegistry registry;
    return registry;
}

void SubcommandRegistry::register_command(const std::string& name,
                                          const std::string& description,
                                          SubcommandFn fn,
                                          int order) {
    handlers_[name] = fn;
    command_list_.push_back({name, description, order});
}

bool SubcommandRegistry::has_command(const std::string& name) const {
    return handlers_.find(name) != handlers_.end();
}

int SubcommandRegistry::run_command(const std::string& name, int argc, char* argv[]) const {
    auto it = handlers_.find(name);
    if (it == handlers_.end()) {
        std::cerr << "Unknown command: " << name << "\n";
        std::cerr << "Run 'dart --help' for usage information.\n";
        return 1;
    }
    return it->second(argc, argv);
}

void SubcommandRegistry::print_help(const char* program_name) const {
    std::cout << "DART v" << DART_VERSION << "\n\n";
    std::cout << "Usage: " << program_name << " <command> [options]\n\n";
    std::cout << "Commands:\n";

    // Sort by workflow order
    auto sorted = command_list_;
    std::sort(sorted.begin(), sorted.end(),
              [](const CommandEntry& a, const CommandEntry& b) {
                  return a.order < b.order;
              });

    size_t max_len = 0;
    for (const auto& cmd : sorted) {
        max_len = std::max(max_len, cmd.name.length());
    }

    for (const auto& cmd : sorted) {
        std::cout << "  " << cmd.name;
        for (size_t i = cmd.name.length(); i < max_len + 2; ++i) {
            std::cout << ' ';
        }
        std::cout << cmd.description << "\n";
    }

    std::cout << "\nFor help on a specific command: " << program_name << " <command> --help\n";
}

}  // namespace cli
}  // namespace dart
