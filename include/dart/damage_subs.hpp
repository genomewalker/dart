#pragma once
// Damage-consistent amino-acid substitution table (single source).
// target = reference (undamaged), query = observed (potentially damaged).
// damage_class: 'C' for C->T (5' deamination), 'G' for G->A (3' deamination).
// Consumed by cmd_damage_annotate (opportunity counting) and the damage-aware
// substitution-matrix emitter (ss/ds edge set).

#include <cstddef>

namespace dart {
namespace cli {

struct DamageSubstitution {
    char target_aa;
    char query_aa;
    char damage_class;  // 'C' for C->T, 'G' for G->A
    bool high_confidence;
};

static constexpr DamageSubstitution DAMAGE_SUBS[] = {
    // C->T damage (5' end)
    // Codon position 1: C->T
    {'R', 'W', 'C', true},    // CGG->TGG
    {'R', 'C', 'C', true},    // CGC->TGC, CGT->TGT
    {'R', '*', 'C', true},    // CGA->TGA
    {'Q', '*', 'C', true},    // CAA->TAA, CAG->TAG
    {'H', 'Y', 'C', true},    // CAC->TAC, CAT->TAT
    {'P', 'S', 'C', false},   // CCN->TCN
    {'P', 'L', 'C', false},   // CCG->CTG
    // Codon position 2: C->T
    {'T', 'I', 'C', false},   // ACN->ATN (except ACG->ATG)
    {'T', 'M', 'C', false},   // ACG->ATG
    {'A', 'V', 'C', false},   // GCN->GTN
    {'S', 'F', 'C', false},   // TCC->TTC, TCT->TTT
    {'S', 'L', 'C', false},   // TCA->TTA, TCG->TTG

    // G->A damage (3' end)
    // Codon position 1: G->A
    {'E', 'K', 'G', true},    // GAA/GAG->AAA/AAG
    {'D', 'N', 'G', true},    // GAC/GAT->AAC/AAT
    {'A', 'T', 'G', false},   // GCN->ACN
    {'G', 'R', 'G', false},   // GGA/GGG->AGA/AGG
    {'G', 'S', 'G', false},   // GGC/GGT->AGC/AGT
    {'V', 'I', 'G', false},   // GTA/GTC/GTT->ATA/ATC/ATT
    {'V', 'M', 'G', false},   // GTG->ATG
    // Codon position 2: G->A
    {'R', 'K', 'G', true},    // AGA/AGG->AAA/AAG
    {'R', 'Q', 'G', true},    // CGA/CGG->CAA/CAG
    {'R', 'H', 'G', true},    // CGC/CGT->CAC/CAT
    {'G', 'E', 'G', true},    // GGA/GGG->GAA/GAG
    {'G', 'D', 'G', true},    // GGC/GGT->GAC/GAT
    {'S', 'N', 'G', false},   // AGC/AGT->AAC/AAT
    {'C', 'Y', 'G', false},   // TGC/TGT->TAC/TAT
    {'W', '*', 'G', false},   // TGG->TAG/TGA

    // X-masked positions (from search_protein): X was originally damaged AA
    // DART X-masks: (1) damage-convertible stops, (2) W/C likely from R at 5' terminus
    // High-confidence because DART only X-masks when damage probability is high
    {'Q', 'X', 'C', true},    // CAA/CAG->TAA/TAG (stop X-masked)
    {'R', 'X', 'C', true},    // CGA->TGA (stop) OR CGG->TGG (W) OR CGC/CGT->TGC/TGT (C)
    {'W', 'X', 'G', true},    // TGG->TGA/TAG (stop X-masked)
};

static constexpr size_t NUM_DAMAGE_SUBS = sizeof(DAMAGE_SUBS) / sizeof(DAMAGE_SUBS[0]);

}  // namespace cli
}  // namespace dart
