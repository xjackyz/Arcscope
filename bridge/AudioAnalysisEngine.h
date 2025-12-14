//
//  AudioAnalysisEngine.h
//  ArcScope Bridge - Complete Audio Analysis with FFT
//
//  Implements arcscope.md sections 3.4.2 and 4.2
//  - Psychoacoustic Loudness (LUFS)
//  - Chroma extraction (FFT + harmonic summation)
//  - Key estimation
//  - Harmonic tension calculation
//  - Dialogue clarity (SNR in speech band)
//  - BPM/tempo detection
//

#import <Foundation/Foundation.h>
#import <Accelerate/Accelerate.h>
#import <AVFoundation/AVFoundation.h>

NS_ASSUME_NONNULL_BEGIN

/**
 * Per-second audio analysis result
 */
@interface AudioAnalysisResult : NSObject
@property (nonatomic, assign) double timeSeconds;
@property (nonatomic, assign) double rms;
@property (nonatomic, assign) double loudness;                // LUFS-like (dB), A-weighted spectrum integration
@property (nonatomic, assign) double spectralFlux;            // Frame-to-frame change
@property (nonatomic, assign) double spectralCentroid;        // Brightness (Hz)
@property (nonatomic, assign) double spectralBalance;         // HF / LF ratio
@property (nonatomic, assign) double zeroCrossingRate;        // Transient indicator
@property (nonatomic, assign) double bpm;
@property (nonatomic, assign) double tempoConfidence;
@property (nonatomic, assign) double rhythmEnergy;
@property (nonatomic, strong) NSArray<NSNumber *> *chroma;    // 12-bin chroma vector
@property (nonatomic, assign) int estimatedKey;               // 0=C, 11=B, -1=unknown
@property (nonatomic, assign) double harmonicTension;
@property (nonatomic, assign) double dialogueClarity;         // SNR in speech band
@property (nonatomic, assign) double speechProbability;
@property (nonatomic, assign) double stereoWidth;             // Mid/side energy ratio [0,1]
@property (nonatomic, assign) double reverbAmount;            // Late/early energy ratio mapped to [0,1]
@property (nonatomic, assign) double dialogueEnergy;          // Speech-band power (raw units)
@property (nonatomic, assign) double musicEnergy;             // Non-speech tonal power (raw units)
@property (nonatomic, assign) double effectsEnergy;           // Residual non-speech power (raw units)
@property (nonatomic, assign) double dialogueDominance;       // DialogueEnergy / total [0,1]
@property (nonatomic, assign) double musicDominance;          // MusicEnergy / total [0,1]
@property (nonatomic, assign) double effectsDominance;        // EffectsEnergy / total [0,1]
@end

/**
 * Complete audio analysis engine
 *
 * Uses vDSP for efficient FFT computation
 */
@interface AudioAnalysisEngine : NSObject

/**
 * Analyze audio asset and extract per-second features
 *
 * @param asset Audio/video asset
 * @param durationSec Total duration
 * @param progressBlock Progress callback (0.0 to 1.0)
 * @return Array of per-second audio features
 */
+ (nullable NSArray<AudioAnalysisResult *> *)analyzeAudio:(AVAsset *)asset
                                               durationSec:(double)durationSec
                                             progressBlock:(nullable void (^)(double progress))progressBlock
                                                     error:(NSError **)error;

/**
 * Extract chroma vector from magnitude spectrum
 *
 * @param magnitudes FFT magnitude spectrum
 * @param sampleRate Audio sample rate
 * @return 12-element chroma vector [C, C#, D, ..., B]
 */
+ (NSArray<NSNumber *> *)extractChroma:(const float *)magnitudes
                                length:(int)length
                            sampleRate:(double)sampleRate;

/**
 * Estimate musical key from chroma vectors
 *
 * @param chromaVectors Array of chroma vectors over time
 * @return Key index (0=C, 11=B) or -1 if unknown
 */
+ (int)estimateKey:(NSArray<NSArray<NSNumber *> *> *)chromaVectors;

/**
 * Calculate harmonic tension relative to key
 *
 * @param chroma 12-bin chroma vector
 * @param key Key index (0=C, 11=B)
 * @return Tension [0,1]: 0=consonant, 1=dissonant
 */
+ (double)calculateHarmonicTension:(NSArray<NSNumber *> *)chroma key:(int)key;

/**
 * Estimate BPM from spectral flux or onset envelope
 *
 * @param onsetStrength Per-frame onset strength
 * @param sampleRate Frame rate (typically 1 Hz for per-second analysis)
 * @return Estimated BPM and confidence
 */
+ (void)estimateBPM:(const double *)onsetStrength
             length:(int)length
         sampleRate:(double)sampleRate
             outBPM:(double *)outBPM
      outConfidence:(double *)outConfidence;

/**
 * Calculate dialogue clarity (SNR in speech band 300-3400 Hz)
 *
 * @param magnitudes FFT magnitude spectrum
 * @param length Spectrum length
 * @param sampleRate Audio sample rate
 * @return Clarity [0,1]: 0=noisy, 1=clear
 */
+ (double)calculateDialogueClarity:(const float *)magnitudes
                            length:(int)length
                        sampleRate:(double)sampleRate;

@end

NS_ASSUME_NONNULL_END
