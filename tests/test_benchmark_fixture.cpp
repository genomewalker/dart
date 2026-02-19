#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream iss(line);
    std::string field;
    while (std::getline(iss, field, '\t')) {
        fields.push_back(field);
    }
    return fields;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <benchmark_summary.tsv>\n";
        return 2;
    }

    std::ifstream in(argv[1]);
    if (!in.is_open()) {
        std::cerr << "Cannot open benchmark fixture: " << argv[1] << "\n";
        return 1;
    }

    std::string header_line;
    if (!std::getline(in, header_line) || header_line.empty()) {
        std::cerr << "Benchmark fixture is empty\n";
        return 1;
    }

    std::vector<std::string> header = split_tsv(header_line);
    std::unordered_map<std::string, size_t> col;
    for (size_t i = 0; i < header.size(); ++i) {
        col[header[i]] = i;
    }

    const char* required[] = {
        "sample_id", "sample_name", "group", "n_total", "auc", "pct_damage_informative"
    };
    for (const char* name : required) {
        if (col.find(name) == col.end()) {
            std::cerr << "Missing required column: " << name << "\n";
            return 1;
        }
    }

    size_t rows = 0;
    bool saw_at = false;
    bool saw_gc = false;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> fields = split_tsv(line);
        if (fields.size() < header.size()) {
            std::cerr << "Malformed row: " << line << "\n";
            return 1;
        }
        ++rows;
        const std::string& group = fields[col["group"]];
        if (group == "AT-rich") saw_at = true;
        if (group == "GC-rich") saw_gc = true;
    }

    if (rows == 0) {
        std::cerr << "No data rows found in benchmark fixture\n";
        return 1;
    }
    if (!saw_at || !saw_gc) {
        std::cerr << "Expected both AT-rich and GC-rich rows in benchmark fixture\n";
        return 1;
    }

    std::cout << "Benchmark fixture OK (" << rows << " rows)\n";
    return 0;
}
