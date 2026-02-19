#ifndef AGP_CLI_SUBCOMMAND_HPP
#define AGP_CLI_SUBCOMMAND_HPP

#include <string>
#include <vector>
#include <functional>
#include <unordered_map>

namespace agp {
namespace cli {

// Subcommand handler function type
using SubcommandFn = std::function<int(int argc, char* argv[])>;

// Registry of available subcommands
class SubcommandRegistry {
public:
    struct CommandEntry {
        std::string name;
        std::string description;
        int order;
    };

    static SubcommandRegistry& instance();

    void register_command(const std::string& name,
                         const std::string& description,
                         SubcommandFn fn,
                         int order = 99);

    bool has_command(const std::string& name) const;
    int run_command(const std::string& name, int argc, char* argv[]) const;

    void print_help(const char* program_name) const;

    const std::vector<CommandEntry>& commands() const {
        return command_list_;
    }

private:
    SubcommandRegistry() = default;
    std::unordered_map<std::string, SubcommandFn> handlers_;
    std::vector<CommandEntry> command_list_;
};

// Forward declarations of subcommand entry points
int cmd_predict(int argc, char* argv[]);
int cmd_validate(int argc, char* argv[]);
int cmd_sample_damage(int argc, char* argv[]);
int cmd_damage_annotate(int argc, char* argv[]);
int cmd_damage_profile(int argc, char* argv[]);
int cmd_hits2emi(int argc, char* argv[]);

}  // namespace cli
}  // namespace agp

#endif  // AGP_CLI_SUBCOMMAND_HPP
