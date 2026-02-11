#pragma once

// Auto-generated score calibration for AGP
// Converts raw prediction scores to P(correct)

namespace agp {

// Score calibration lookup table (isotonic regression)
// Generated from 49 calibration points
// Usage: p_correct = calibrate_score(raw_score)

struct CalibrationPoint {
    float threshold;
    float p_correct;
};

constexpr CalibrationPoint SCORE_CALIBRATION[] = {
    {0.331014f, 0.314378f},
    {0.537581f, 0.380378f},
    {0.573199f, 0.407038f},
    {0.598063f, 0.419830f},
    {0.618333f, 0.432483f},
    {0.635247f, 0.444746f},
    {0.649757f, 0.464004f},
    {0.662599f, 0.479815f},
    {0.674160f, 0.485636f},
    {0.684792f, 0.501468f},
    {0.694836f, 0.518421f},
    {0.704348f, 0.527403f},
    {0.713371f, 0.542602f},
    {0.722052f, 0.559332f},
    {0.730332f, 0.576325f},
    {0.738295f, 0.588264f},
    {0.745896f, 0.602088f},
    {0.753150f, 0.621343f},
    {0.760135f, 0.630918f},
    {0.767089f, 0.640488f},
    {0.773844f, 0.645680f},
    {0.780402f, 0.664852f},
    {0.786782f, 0.676182f},
    {0.792959f, 0.686953f},
    {0.799026f, 0.699812f},
    {0.804961f, 0.706547f},
    {0.810695f, 0.720201f},
    {0.816186f, 0.737823f},
    {0.821373f, 0.750867f},
    {0.826326f, 0.760597f},
    {0.831179f, 0.770481f},
    {0.835956f, 0.780109f},
    {0.840646f, 0.788566f},
    {0.845305f, 0.796049f},
    {0.849968f, 0.804453f},
    {0.854700f, 0.812462f},
    {0.859496f, 0.817053f},
    {0.864495f, 0.823717f},
    {0.869679f, 0.829409f},
    {0.875055f, 0.834829f},
    {0.880803f, 0.842641f},
    {0.887042f, 0.845614f},
    {0.893997f, 0.853936f},
    {0.901933f, 0.854524f},
    {0.911427f, 0.854524f},
    {0.922881f, 0.857012f},
    {0.936939f, 0.861011f},
    {0.954653f, 0.876463f},
    {0.979419f, 0.910779f},
};
constexpr size_t SCORE_CALIBRATION_SIZE = 49;

inline float calibrate_score(float raw_score) {
    if (raw_score <= SCORE_CALIBRATION[0].threshold) return SCORE_CALIBRATION[0].p_correct;
    if (raw_score >= SCORE_CALIBRATION[SCORE_CALIBRATION_SIZE-1].threshold) return SCORE_CALIBRATION[SCORE_CALIBRATION_SIZE-1].p_correct;
    
    // Binary search for bracket
    size_t lo = 0, hi = SCORE_CALIBRATION_SIZE - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (raw_score < SCORE_CALIBRATION[mid].threshold) hi = mid;
        else lo = mid;
    }
    
    // Linear interpolation
    float t = (raw_score - SCORE_CALIBRATION[lo].threshold) / 
              (SCORE_CALIBRATION[hi].threshold - SCORE_CALIBRATION[lo].threshold);
    return SCORE_CALIBRATION[lo].p_correct + t * (SCORE_CALIBRATION[hi].p_correct - SCORE_CALIBRATION[lo].p_correct);
}

} // namespace agp
