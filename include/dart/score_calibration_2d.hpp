#pragma once

// Auto-generated 2D score-length calibration for AGP
// Converts (raw_score, length) to P(correct)
// Phase 2: Length-aware calibration

#include <cstddef>

namespace dart {

// 2D Score-Length calibration lookup table
// P(correct|score, length) with length-aware corrections
// Generated from empirical validation data

namespace detail {

constexpr size_t N_SCORE_BINS = 20;
constexpr size_t N_LENGTH_BINS = 8;

// Score bin edges
constexpr float SCORE_EDGES[] = {
    0.3000f, 0.3350f, 0.3700f, 0.4050f, 0.4400f, 0.4750f, 0.5100f, 0.5450f, 0.5800f, 0.6150f, 0.6500f, 0.6850f, 0.7200f, 0.7550f, 0.7900f, 0.8250f, 0.8600f, 0.8950f, 0.9300f, 0.9650f, 1.0000f
};

// Length bin edges
constexpr float LENGTH_EDGES[] = {
    30.0f, 45.0f, 60.0f, 80.0f, 100.0f, 130.0f, 160.0f, 200.0f, 300.0f
};

// P(correct) matrix [length_bin][score_bin]
constexpr float CALIBRATION_2D[8][20] = {
    {0.0000f, 0.1579f, 0.3167f, 0.3167f, 0.3167f, 0.3167f, 0.3340f, 0.3645f, 0.3983f, 0.4274f, 0.4663f, 0.5074f, 0.5533f, 0.6031f, 0.6633f, 0.7354f, 0.7910f, 0.8190f, 0.8358f, 0.8718f},  // length 30-45
    {0.0000f, 0.1579f, 0.2797f, 0.3240f, 0.3240f, 0.3431f, 0.3523f, 0.4001f, 0.4228f, 0.4579f, 0.4899f, 0.5365f, 0.5949f, 0.6481f, 0.7153f, 0.7857f, 0.8388f, 0.8626f, 0.8700f, 0.8954f},  // length 45-60
    {0.0000f, 0.1579f, 0.3048f, 0.3152f, 0.3152f, 0.3198f, 0.3577f, 0.4129f, 0.4362f, 0.4555f, 0.5057f, 0.5446f, 0.6254f, 0.6923f, 0.7574f, 0.8337f, 0.8743f, 0.8956f, 0.9015f, 0.9175f},  // length 60-80
    {0.0000f, 0.1579f, 0.2788f, 0.2951f, 0.3150f, 0.3150f, 0.3511f, 0.4415f, 0.4613f, 0.4613f, 0.4933f, 0.5456f, 0.6403f, 0.7308f, 0.7965f, 0.8632f, 0.9044f, 0.9214f, 0.9241f, 0.9295f},  // length 80-100
    {0.0000f, 0.1579f, 0.1753f, 0.2670f, 0.2905f, 0.2905f, 0.3518f, 0.4542f, 0.4687f, 0.4687f, 0.4687f, 0.5274f, 0.6353f, 0.7375f, 0.8015f, 0.8630f, 0.8909f, 0.9088f, 0.9132f, 0.9216f},  // length 100-130
    {0.0000f, 0.1579f, 0.1746f, 0.1875f, 0.2099f, 0.2099f, 0.2570f, 0.3324f, 0.3532f, 0.3532f, 0.3532f, 0.3630f, 0.4444f, 0.5141f, 0.5629f, 0.6072f, 0.6233f, 0.6285f, 0.6285f, 0.6427f},  // length 130-160
    {0.0000f, 0.1579f, 0.1935f, 0.1935f, 0.2403f, 0.2403f, 0.2510f, 0.3049f, 0.3356f, 0.3356f, 0.3523f, 0.3711f, 0.4571f, 0.5108f, 0.5677f, 0.6024f, 0.6146f, 0.6295f, 0.6295f, 0.6295f},  // length 160-200
    {0.0000f, 0.1579f, 0.2626f, 0.2900f, 0.2939f, 0.3091f, 0.3387f, 0.3879f, 0.4140f, 0.4380f, 0.4766f, 0.5207f, 0.5822f, 0.6447f, 0.7123f, 0.7869f, 0.8330f, 0.8546f, 0.8630f, 0.8872f},  // length 200-300
};

} // namespace detail

inline float calibrate_score_2d(float raw_score, size_t length) {
    // Find length bin
    size_t l_bin = 0;
    for (size_t i = 1; i < detail::N_LENGTH_BINS; ++i) {
        if (static_cast<float>(length) >= detail::LENGTH_EDGES[i]) l_bin = i;
        else break;
    }
    if (l_bin >= detail::N_LENGTH_BINS) l_bin = detail::N_LENGTH_BINS - 1;

    // Find score bin
    size_t s_bin = 0;
    for (size_t i = 1; i < detail::N_SCORE_BINS; ++i) {
        if (raw_score >= detail::SCORE_EDGES[i]) s_bin = i;
        else break;
    }
    if (s_bin >= detail::N_SCORE_BINS) s_bin = detail::N_SCORE_BINS - 1;

    return detail::CALIBRATION_2D[l_bin][s_bin];
}

} // namespace dart
