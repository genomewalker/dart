#pragma once

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>

namespace dart {
namespace log_utils {

inline std::string format_duration_ms(int64_t ms) {
    if (ms < 0) ms = 0;
    if (ms < 1000) {
        return std::to_string(ms) + " ms";
    }

    const int64_t total_seconds = ms / 1000;
    if (total_seconds < 60) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1)
            << (static_cast<double>(ms) / 1000.0) << " s";
        return oss.str();
    }

    const int64_t seconds = total_seconds % 60;
    const int64_t total_minutes = total_seconds / 60;
    if (total_minutes < 60) {
        const int64_t minutes = total_minutes;
        return std::to_string(minutes) + "m " + std::to_string(seconds) + "s";
    }

    const int64_t minutes = total_minutes % 60;
    const int64_t total_hours = total_minutes / 60;
    if (total_hours < 24) {
        const int64_t hours = total_hours;
        return std::to_string(hours) + "h " + std::to_string(minutes) + "m " +
               std::to_string(seconds) + "s";
    }

    const int64_t hours = total_hours % 24;
    const int64_t days = total_hours / 24;
    return std::to_string(days) + "d " + std::to_string(hours) + "h " +
           std::to_string(minutes) + "m";
}

template <typename Clock, typename DurA, typename DurB>
inline std::string format_elapsed(
    const std::chrono::time_point<Clock, DurA>& start,
    const std::chrono::time_point<Clock, DurB>& end) {
    const auto ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    return format_duration_ms(ms);
}

}  // namespace log_utils
}  // namespace dart

