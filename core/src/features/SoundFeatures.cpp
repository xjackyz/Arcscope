#include "arcscope/features/SoundFeatures.h"
#include "arcscope/Normalization.h"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace arcscope {
namespace features {

std::vector<AudioFeatures> SoundAnalyzer::analyze(const std::vector<FeatureSample>& samples) const {
    std::vector<AudioFeatures> result;
    result.reserve(samples.size());

    for (size_t i = 0; i < samples.size(); ++i) {
        AudioFeatures af;
        af.timeSeconds = samples[i].timeSeconds;

        // Observation-level audio descriptors are computed in Bridge (AudioAnalysisEngine) and passed through.
        // Core does NOT approximate FFT-derived features from 1Hz aggregates.
        af.rms = samples[i].audioRms;
        af.loudness = samples[i].audioLoudness;
        af.spectralFlux = samples[i].audioSpectralFlux;
        af.spectralCentroid = samples[i].audioSpectralCentroid;
        af.spectralBalance = samples[i].spectralBalance;
        af.zeroCrossingRate = samples[i].audioZeroCrossingRate;

        af.bpm = samples[i].audioBpm;
        af.tempoConfidence = samples[i].audioTempoConfidence;
        af.rhythmEnergy = samples[i].audioRhythmEnergy;
        af.chroma = samples[i].audioChroma.empty() ? std::vector<double>(12, 0.0) : samples[i].audioChroma;
        af.estimatedKey = samples[i].audioEstimatedKey;
        af.harmonicTension = samples[i].audioHarmonicTension;
        af.dialogueClarity = samples[i].audioDialogueClarity;
        af.speechProbability = samples[i].audioSpeechProbability;
        af.stereoWidth = samples[i].audioStereoWidth;
        af.reverbAmount = samples[i].audioReverbAmount;
        af.dialogueEnergy = samples[i].audioDialogueEnergy;
        af.musicEnergy = samples[i].audioMusicEnergy;
        af.effectsEnergy = samples[i].audioEffectsEnergy;
        af.dialogueDominance = samples[i].audioDialogueDominance;
        af.musicDominance = samples[i].audioMusicDominance;
        af.effectsDominance = samples[i].audioEffectsDominance;

        result.push_back(af);
    }

    return result;
}

// ========== DELETED FUNCTIONS ==========
// The following functions have been removed because they attempted to estimate
// features from insufficient input data (信息论层面的缺陷):
//
// - estimateBPM / estimateBPMWithConfidence
//   Reason: BPM needs high-rate onset detection (≈86Hz), not 1Hz transients
//
// - estimateChromaFromSpectrum
//   Reason: Chroma needs FFT with pitch class resolution, not spectralBalance
//
// - estimateKey
//   Reason: Key estimation needs real chroma vectors from FFT
//
// - calculateHarmonicTension
//   Reason: Harmonic tension needs real chroma and key
//
// - computeRhythmEnergy
//   Reason: Rhythm energy needs per-beat phase alignment, not global average
//
// These features should be computed in Bridge (AudioAnalysisEngine) and passed
// to Core via an expanded FeatureSample structure.
// ========================================

}  // namespace features
}  // namespace arcscope
