//
//  AudioAnalysisEngine.mm
//  ArcScope Bridge - Complete Audio Analysis Implementation
//
//   修正版（符合 ArcScope 1.0 标准）
//  核心修正：
//  1. 强制 mono 输出（避免多通道采样率计算错误）
//  2. FFT 实信号输入修正（fftSize=1024, 与 ctoz 填法对齐，方案 A）
//  3. 使用 vDSP_zvabs 得到 magnitude（而非 zvmags 的 power）
//  4. hop 级 onset 序列用于 BPM 估计（~86Hz 采样率，支持 30-180 BPM）
//  5. key 编码支持大小调（0-11=major, 12-23=minor）
//  6. harmonic tension 根据大小调选择正确三和弦
//  7. dialogue clarity 使用功率（magnitude²）+ SNR(dB) 映射
//

#import "AudioAnalysisEngine.h"
#import <cmath>
#import <vector>
#import <algorithm>
#import <numeric>

@implementation AudioAnalysisResult
@end

@implementation AudioAnalysisEngine

#pragma mark - Main Analysis Pipeline

static inline double clamp01d(double x) {
    if (!std::isfinite(x)) return 0.0;
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

static inline double safe_div(double a, double b, double fallback = 0.0) {
    return (std::abs(b) > 1e-12) ? (a / b) : fallback;
}

namespace {

struct Biquad {
    double b0{0.0}, b1{0.0}, b2{0.0};
    double a1{0.0}, a2{0.0};
    double z1{0.0}, z2{0.0};

    inline double process(double x) {
        // Direct Form II transposed
        const double y = b0 * x + z1;
        z1 = b1 * x - a1 * y + z2;
        z2 = b2 * x - a2 * y;
        return y;
    }
};

Biquad design_highpass(double fs, double f0, double Q) {
    const double w0 = 2.0 * M_PI * f0 / fs;
    const double cosw = std::cos(w0);
    const double sinw = std::sin(w0);
    const double alpha = sinw / (2.0 * Q);

    const double b0 =  (1.0 + cosw) * 0.5;
    const double b1 = -(1.0 + cosw);
    const double b2 =  (1.0 + cosw) * 0.5;
    const double a0 =  1.0 + alpha;
    const double a1 = -2.0 * cosw;
    const double a2 =  1.0 - alpha;

    Biquad q;
    q.b0 = b0 / a0;
    q.b1 = b1 / a0;
    q.b2 = b2 / a0;
    q.a1 = a1 / a0;
    q.a2 = a2 / a0;
    return q;
}

Biquad design_highshelf(double fs, double f0, double gain_db, double slope) {
    const double A = std::pow(10.0, gain_db / 40.0);
    const double w0 = 2.0 * M_PI * f0 / fs;
    const double cosw = std::cos(w0);
    const double sinw = std::sin(w0);
    const double alpha = (sinw / 2.0) * std::sqrt((A + 1.0 / A) * (1.0 / slope - 1.0) + 2.0);

    const double b0 =      A * ((A + 1.0) + (A - 1.0) * cosw + 2.0 * std::sqrt(A) * alpha);
    const double b1 = -2.0 * A * ((A - 1.0) + (A + 1.0) * cosw);
    const double b2 =      A * ((A + 1.0) + (A - 1.0) * cosw - 2.0 * std::sqrt(A) * alpha);
    const double a0 =           (A + 1.0) - (A - 1.0) * cosw + 2.0 * std::sqrt(A) * alpha;
    const double a1 =  2.0 * ((A - 1.0) - (A + 1.0) * cosw);
    const double a2 =           (A + 1.0) - (A - 1.0) * cosw - 2.0 * std::sqrt(A) * alpha;

    Biquad q;
    q.b0 = b0 / a0;
    q.b1 = b1 / a0;
    q.b2 = b2 / a0;
    q.a1 = a1 / a0;
    q.a2 = a2 / a0;
    return q;
}

struct KWeighting {
    Biquad hpL, hsL;
    Biquad hpR, hsR;

    explicit KWeighting(double fs)
        : hpL(design_highpass(fs, /*f0=*/40.0, /*Q=*/0.5)),
          hsL(design_highshelf(fs, /*f0=*/1500.0, /*gain_db=*/4.0, /*slope=*/1.0)),
          hpR(design_highpass(fs, /*f0=*/40.0, /*Q=*/0.5)),
          hsR(design_highshelf(fs, /*f0=*/1500.0, /*gain_db=*/4.0, /*slope=*/1.0)) {}

    inline void processStereo(double inL, double inR, double& outL, double& outR) {
        outL = hsL.process(hpL.process(inL));
        outR = hsR.process(hpR.process(inR));
    }
};

double lufs_from_mean_square(double ms) {
    // ITU-R BS.1770: LUFS = -0.691 + 10*log10(mean_square)
    if (!(ms > 0.0) || !std::isfinite(ms)) return -120.0;
    return -0.691 + 10.0 * std::log10(ms);
}

double weighted_median(std::vector<std::pair<double,double>> vw) {
    if (vw.empty()) return 0.0;
    std::sort(vw.begin(), vw.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    double total = 0.0;
    for (const auto& p : vw) total += std::max(0.0, p.second);
    if (!(total > 0.0)) return vw[vw.size() / 2].first;
    double acc = 0.0;
    const double half = 0.5 * total;
    for (const auto& p : vw) {
        acc += std::max(0.0, p.second);
        if (acc >= half) return p.first;
    }
    return vw.back().first;
}

double median_copy(std::vector<double> v) {
    if (v.empty()) return 0.0;
    const size_t n = v.size();
    std::nth_element(v.begin(), v.begin() + n / 2, v.end());
    double med = v[n / 2];
    if (n % 2 == 0) {
        const auto it = std::max_element(v.begin(), v.begin() + n / 2);
        med = 0.5 * (med + *it);
    }
    return med;
}

double mad_copy(const std::vector<double>& v, double med) {
    if (v.empty()) return 0.0;
    std::vector<double> dev(v.size());
    for (size_t i = 0; i < v.size(); ++i) dev[i] = std::abs(v[i] - med);
    return median_copy(std::move(dev));
}

double quantile_copy(std::vector<double> v, double q01) {
    if (v.empty()) return 0.0;
    if (!std::isfinite(q01)) q01 = 0.5;
    q01 = std::clamp(q01, 0.0, 1.0);
    const size_t n = v.size();
    if (n == 1) return v[0];
    const double pos = q01 * (double)(n - 1);
    const size_t k0 = (size_t)std::floor(pos);
    const size_t k1 = (size_t)std::ceil(pos);
    std::nth_element(v.begin(), v.begin() + k0, v.end());
    const double a0 = v[k0];
    if (k1 == k0) return a0;
    std::nth_element(v.begin(), v.begin() + k1, v.end());
    const double a1 = v[k1];
    return a0 + (pos - (double)k0) * (a1 - a0);
}

std::vector<double> robust_z_copy(const std::vector<double>& x) {
    if (x.empty()) return {};
    const double med = median_copy(std::vector<double>(x.begin(), x.end()));
    const double mad = mad_copy(x, med);
    const double den = (std::isfinite(mad) ? mad : 0.0) + 1e-6;
    std::vector<double> z(x.size());
    for (size_t i = 0; i < x.size(); ++i) z[i] = (x[i] - med) / den;
    return z;
}

double sigmoid01(double z) {
    if (z >= 40.0) return 1.0;
    if (z <= -40.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-z));
}

std::vector<double> robust_sigmoid01_copy(const std::vector<double>& x) {
    auto z = robust_z_copy(x);
    for (auto& v : z) v = sigmoid01(v);
    return z;
}

double mean_abs_slope_1hz(const std::vector<double>& z) {
    if (z.size() < 2) return 0.0;
    double acc = 0.0;
    for (size_t i = 1; i < z.size(); ++i) acc += std::abs(z[i] - z[i - 1]);
    return acc / (double)(z.size() - 1);
}

double adaptive_window_size(double lambda) {
    // Keep in sync with core/include/arcscope/CoreContract.h (contract::adaptive_window_size).
    const double W_BASE = 12.0;
    const double W_LAMBDA_SCALE = 2.0;
    const double W_MIN = 4.0;
    const double W_MAX = 30.0;
    const double w = W_BASE * std::exp(-W_LAMBDA_SCALE * lambda);
    if (w < W_MIN) return W_MIN;
    if (w > W_MAX) return W_MAX;
    return w;
}

struct TempoEstimate {
    double bpm{0.0};
    double confidence{0.0}; // [0,1]
    double beatRatio{0.0};  // [0,1] energy share around the peak
};

TempoEstimate estimate_tempo_from_onset_window(const std::vector<double>& onset,
                                               double hopRate,
                                               double bpmMin,
                                               double bpmMax) {
    TempoEstimate out;
    if (onset.size() < 8) return out;
    if (!(hopRate > 1e-6)) return out;
    if (!(bpmMax > bpmMin && bpmMin > 0.0)) return out;

    const int minLag = std::max(1, (int)std::llround(hopRate * 60.0 / bpmMax));
    const int maxLag = std::max(minLag + 1, (int)std::llround(hopRate * 60.0 / bpmMin));
    if ((int)onset.size() <= maxLag + 2) return out;

    // Precompute prefix sums of squares for quick normalization.
    std::vector<double> x(onset.size(), 0.0);
    for (size_t i = 0; i < onset.size(); ++i) x[i] = std::max(0.0, onset[i]);
    std::vector<double> ps2(x.size() + 1, 0.0);
    for (size_t i = 0; i < x.size(); ++i) ps2[i + 1] = ps2[i] + x[i] * x[i];

    auto sumsq = [&](size_t a, size_t b) -> double {
        if (b <= a) return 0.0;
        return ps2[b] - ps2[a];
    };

    std::vector<double> ac;
    ac.reserve((size_t)(maxLag - minLag + 1));
    for (int lag = minLag; lag <= maxLag; ++lag) {
        const size_t n = x.size() - (size_t)lag;
        double dot = 0.0;
        for (size_t i = 0; i < n; ++i) dot += x[i] * x[i + (size_t)lag];
        const double n0 = sumsq(0, n);
        const double n1 = sumsq((size_t)lag, (size_t)lag + n);
        const double den = std::sqrt(std::max(1e-24, n0 * n1));
        ac.push_back(dot / den);
    }

    // Best lag = max autocorrelation.
    int bestIdx = 0;
    for (int i = 1; i < (int)ac.size(); ++i) if (ac[(size_t)i] > ac[(size_t)bestIdx]) bestIdx = i;
    const int bestLag = minLag + bestIdx;
    out.bpm = 60.0 * hopRate / (double)bestLag;

    // Confidence: peak prominence relative to robust center.
    const double med = median_copy(ac);
    const double mad = mad_copy(ac, med) + 1e-6;
    const double peak = ac[(size_t)bestIdx];
    out.confidence = clamp01d((peak - med) / (3.0 * mad));

    // Beat-neighborhood energy ratio: contiguous region around peak where robust z > 0 (adaptive width).
    const auto z = robust_z_copy(ac);
    int lo = bestIdx;
    int hi = bestIdx;
    while (lo - 1 >= 0 && z[(size_t)(lo - 1)] > 0.0) --lo;
    while (hi + 1 < (int)z.size() && z[(size_t)(hi + 1)] > 0.0) ++hi;

    double sumAll = 0.0;
    double sumPeak = 0.0;
    for (int i = 0; i < (int)ac.size(); ++i) {
        const double v = std::max(0.0, ac[(size_t)i]);
        sumAll += v;
        if (i >= lo && i <= hi) sumPeak += v;
    }
    out.beatRatio = (sumAll > 1e-12) ? clamp01d(sumPeak / sumAll) : 0.0;
    return out;
}

} // namespace

+ (nullable NSArray<AudioAnalysisResult *> *)analyzeAudio:(AVAsset *)asset
                                               durationSec:(double)durationSec
                                             progressBlock:(nullable void (^)(double progress))progressBlock
                                                     error:(NSError **)error {

    NSMutableArray<AudioAnalysisResult *> *results = [NSMutableArray array];

    // Get first audio track
    NSArray<AVAssetTrack *> *audioTracks = [asset tracksWithMediaType:AVMediaTypeAudio];
    if (audioTracks.count == 0) {
        if (error) {
            *error = [NSError errorWithDomain:@"AudioAnalysisEngine"
                                         code:-1
                                     userInfo:@{NSLocalizedDescriptionKey: @"No audio track found"}];
        }
        return nil;
    }

    AVAssetTrack *audioTrack = audioTracks[0];

    // Setup reader
    AVAssetReader *reader = [[AVAssetReader alloc] initWithAsset:asset error:error];
    if (!reader) return nil;

    // Output stereo PCM; analysis uses mono downmix where needed, but StereoWidth needs L/R.
    NSDictionary *outputSettings = @{
        AVFormatIDKey: @(kAudioFormatLinearPCM),
        AVLinearPCMBitDepthKey: @(32),
        AVLinearPCMIsFloatKey: @(YES),
        AVLinearPCMIsBigEndianKey: @(NO),
        AVLinearPCMIsNonInterleaved: @(NO),
        AVSampleRateKey: @(44100.0),
        AVNumberOfChannelsKey: @(2)
    };

    AVAssetReaderTrackOutput *output = [[AVAssetReaderTrackOutput alloc] initWithTrack:audioTrack
                                                                          outputSettings:outputSettings];
    [reader addOutput:output];

    if (![reader startReading]) {
        if (error) *error = reader.error;
        return nil;
    }

    const double sampleRate = 44100.0;
    const int channels = 2;
    const int samplesPerSecondPerChannel = (int)sampleRate;
    const int fftSize = 1024;  // ✅ 修正：与当前 ctoz 填法对齐（方案 A）
    const int hopSize = 512;

    // FFT Setup - 使用统一的cleanup路径避免资源泄漏
    FFTSetup fftSetup = NULL;
    DSPSplitComplex splitComplex;
    splitComplex.realp = NULL;
    splitComplex.imagp = NULL;
    float *window = NULL;
    float *magnitudes = NULL;
    float *previousMagnitudes = NULL;

    fftSetup = vDSP_create_fftsetup(vDSP_Length(log2(fftSize)), kFFTRadix2);
    if (!fftSetup) {
        if (error) {
            *error = [NSError errorWithDomain:@"AudioAnalysisEngine"
                                         code:-2
                                     userInfo:@{NSLocalizedDescriptionKey: @"FFT setup failed"}];
        }
        return nil;
    }

    splitComplex.realp = (float *)malloc(fftSize/2 * sizeof(float));
    splitComplex.imagp = (float *)malloc(fftSize/2 * sizeof(float));
    window = (float *)malloc(fftSize * sizeof(float));
    magnitudes = (float *)malloc(fftSize/2 * sizeof(float));
    previousMagnitudes = (float *)calloc(fftSize/2, sizeof(float));

    if (!splitComplex.realp || !splitComplex.imagp || !window || !magnitudes || !previousMagnitudes) {
        if (error) {
            *error = [NSError errorWithDomain:@"AudioAnalysisEngine"
                                         code:-3
                                     userInfo:@{NSLocalizedDescriptionKey: @"Memory allocation failed"}];
        }
        goto cleanup;
    }

    vDSP_hann_window(window, fftSize, vDSP_HANN_NORM);

    // ✅ 使用作用域块避免C++ goto跨越变量初始化的问题
	    {
	        std::vector<double> rmsPerSecond;
            std::vector<double> loudnessPerSecond; // LUFS (dB), BS.1770 K-weighting + gating (per second)
	        std::vector<double> spectralFluxPerSecond;
	        std::vector<double> spectralCentroidPerSecond;
	        std::vector<double> spectralBalancePerSecond;
	        std::vector<double> zcRatePerSecond;
	        std::vector<std::vector<double>> chromaPerSecond;
	        std::vector<double> dialogueClarityPerSecond;
	        std::vector<double> stereoWidthPerSecond;
            std::vector<double> reverbAmountPerSecond;
            std::vector<double> dialogueEnergyPerSecond;
            std::vector<double> musicEnergyPerSecond;
            std::vector<double> effectsEnergyPerSecond;
            std::vector<double> dialogueDominancePerSecond;
            std::vector<double> musicDominancePerSecond;
            std::vector<double> effectsDominancePerSecond;
            std::vector<double> rhythmEnergyPerSecond;

        // ✅ FFT buffer 移到循环外复用，避免每秒重新分配
        std::vector<float> interleaved(fftSize, 0.0f);
        std::vector<float> windowed(fftSize, 0.0f);

	        // ✅ hop级onset序列，用于BPM估计（~86Hz采样率）
	        std::vector<double> onsetFrames;
            std::vector<std::vector<double>> onsetPerSecond;
            std::vector<size_t> onsetStartIndexPerSecond;

	        int currentSecond = 0;
	        std::vector<float> secondBufferL;
	        std::vector<float> secondBufferR;
	        std::vector<float> secondBufferMono;

            // BS.1770 loudness: K-weighting + 400ms blocks (75% overlap) + gating.
            KWeighting kweight(sampleRate);
            const int blockSizeSamples = (int)std::llround(0.400 * sampleRate);
            const int blockHopSamples  = (int)std::llround(0.100 * sampleRate);
            std::vector<double> energyRing((size_t)std::max(1, blockSizeSamples), 0.0);
            size_t ringPos = 0;
            size_t ringFill = 0;
            double ringSum = 0.0;
            int hopCountdown = blockHopSamples;
            int64_t sampleIndex = 0; // per-channel sample index (stereo frame index)
            std::vector<double> blockMs;
            std::vector<double> blockTimeCenter;
            std::vector<int> blockSecondIndex;

            // Audio role separation (doc 4.2.6): hop-level descriptors for per-second aggregation.
            std::vector<double> hopTotalPower;
            std::vector<double> hopSpeechPower;
            std::vector<double> hopSpeechRatio;
            std::vector<double> hopTonalness;

    // Process audio samples
    while (reader.status == AVAssetReaderStatusReading) {
        CMSampleBufferRef sampleBuffer = [output copyNextSampleBuffer];
        if (!sampleBuffer) break;

        CMBlockBufferRef blockBuffer = CMSampleBufferGetDataBuffer(sampleBuffer);
        if (!blockBuffer) {
            CFRelease(sampleBuffer);
            continue;
        }

        size_t length = CMBlockBufferGetDataLength(blockBuffer);
        char *data = NULL;
        CMBlockBufferGetDataPointer(blockBuffer, 0, NULL, NULL, &data);

	        float *samples = (float *)data;
	        size_t numSamples = length / sizeof(float);

	        // Interleaved stereo float32: L,R,L,R,...
	        for (size_t i = 0; i + 1 < numSamples; i += (size_t)channels) {
                const double inL = (double)samples[i + 0];
                const double inR = (double)samples[i + 1];
	            secondBufferL.push_back((float)inL);
	            secondBufferR.push_back((float)inR);

                // K-weight filtering for loudness
                double yL = 0.0, yR = 0.0;
                kweight.processStereo(inL, inR, yL, yR);
                const double e = yL * yL + yR * yR;

                if (ringFill < energyRing.size()) {
                    energyRing[ringFill] = e;
                    ringSum += e;
                    ringFill++;
                } else {
                    ringSum -= energyRing[ringPos];
                    energyRing[ringPos] = e;
                    ringSum += e;
                    ringPos = (ringPos + 1) % energyRing.size();
                }

                if (ringFill == energyRing.size()) {
                    hopCountdown--;
                    if (hopCountdown <= 0) {
                        const double ms = ringSum / (double)energyRing.size();
                        const double center = ((double)sampleIndex - 0.5 * (double)(energyRing.size() - 1)) / sampleRate;
                        blockMs.push_back(ms);
                        blockTimeCenter.push_back(center);
                        blockSecondIndex.push_back((int)std::floor(center));
                        hopCountdown = blockHopSamples;
                    }
                }
                sampleIndex++;

	            // Process by second
	            if (secondBufferL.size() >= (size_t)samplesPerSecondPerChannel) {
                    const size_t n = secondBufferL.size();
                    secondBufferMono.resize(n);
                    for (size_t k = 0; k < n; ++k) {
                        secondBufferMono[k] = 0.5f * (secondBufferL[k] + secondBufferR[k]);
                    }

	                // RMS
	                double rms = [self calculateRMS:secondBufferMono.data() length:(int)secondBufferMono.size()];
	                rmsPerSecond.push_back(rms);

	                // ZCR
	                double zcr = [self calculateZeroCrossingRate:secondBufferMono.data() length:(int)secondBufferMono.size()];
	                zcRatePerSecond.push_back(zcr);

	                // FFT-based features
	                std::vector<double> chromaSum(12, 0.0);
	                double fluxSum = 0.0;
	                double centroidSum = 0.0;
	                double balanceSum = 0.0;
	                double claritySum = 0.0;
                    std::vector<double> frameEnergy;
	                int numFrames = 0;
                    std::vector<double> secondOnsets;

                    frameEnergy.clear();
                    const size_t onsetStart = onsetFrames.size();
	                for (int offset = 0; offset + fftSize <= (int)secondBufferMono.size(); offset += hopSize) {
	                    // Apply window
	                    float *frame = secondBufferMono.data() + offset;
	                    vDSP_vmul(frame, 1, window, 1, windowed.data(), 1, fftSize);

                    // ✅ 把实信号转成复数交错格式：even samples→real, odd samples→imag
                    // vDSP_fft_zrip 经典 packing: z[n] = x[2n] + i*x[2n+1]
                    for (int i = 0; i < fftSize/2; i++) {
                        interleaved[2*i]   = windowed[2*i];     // even samples
                        interleaved[2*i+1] = windowed[2*i + 1]; // odd samples
                    }

                    // Convert to split complex
                    vDSP_ctoz((DSPComplex *)interleaved.data(), 2, &splitComplex, 1, fftSize/2);

                    // Forward FFT
                    vDSP_fft_zrip(fftSetup, &splitComplex, 1, vDSP_Length(log2(fftSize)), FFT_FORWARD);

	                    // ✅ 使用 vDSP_zvabs 得到幅度（magnitude），而非 zvmags（power=magnitude²）
	                    vDSP_zvabs(&splitComplex, 1, magnitudes, 1, fftSize/2);

                        // Energy partition + tonalness (observation-level)
                        const double freqPerBin = sampleRate / (double)fftSize;
                        int speechStartBin = (int)(300.0 / freqPerBin);
                        int speechEndBin = (int)(3400.0 / freqPerBin);
                        speechStartBin = std::clamp(speechStartBin, 0, fftSize/2 - 1);
                        speechEndBin = std::clamp(speechEndBin, 0, fftSize/2 - 1);

                        double totalPower = 0.0;
                        double speechPower = 0.0;
                        double logSum = 0.0;
                        for (int k = 0; k < fftSize/2; k++) {
                            const double mag = magnitudes[k];
                            const double pwr = mag * mag;
                            totalPower += pwr;
                            logSum += std::log(std::max(pwr, 1e-20));
                            if (k >= speechStartBin && k <= speechEndBin) speechPower += pwr;
                        }
                        const double bins = std::max(1, fftSize/2);
                        const double arith = totalPower / bins;
                        const double geo = std::exp(logSum / bins);
                        const double flatness = geo / std::max(arith, 1e-20);
                        const double tonalness = clamp01d(1.0 - flatness);

                        frameEnergy.push_back(totalPower);

                    // Spectral flux（用magnitude计算）
                    double flux = 0.0;
                    for (int k = 0; k < fftSize/2; k++) {
                        double diff = magnitudes[k] - previousMagnitudes[k];
                        if (diff > 0) flux += diff;
                        previousMagnitudes[k] = magnitudes[k];
                    }
                    fluxSum += flux;

	                    // ✅ 保存hop级onset强度，用于BPM估计
                    onsetFrames.push_back(flux);
                    secondOnsets.push_back(flux);

                        // Audio role separation hop descriptors (aligned with onsetFrames order)
                        const double tp = std::max(0.0, totalPower);
                        const double sp = std::clamp(speechPower, 0.0, tp);
                        hopTotalPower.push_back(tp);
                        hopSpeechPower.push_back(sp);
                        hopSpeechRatio.push_back(safe_div(sp, tp, 0.0));
                        hopTonalness.push_back(tonalness);

                    // Spectral centroid
                    double centroid = [self calculateSpectralCentroid:magnitudes length:fftSize/2 sampleRate:sampleRate];
                    centroidSum += centroid;

                    // Spectral balance
                    double balance = [self calculateSpectralBalance:magnitudes length:fftSize/2];
                    balanceSum += balance;

                    // Chroma
                    NSArray<NSNumber *> *chroma = [self extractChroma:magnitudes length:fftSize/2 sampleRate:sampleRate];
                    for (int c = 0; c < 12; c++) {
                        chromaSum[c] += [chroma[c] doubleValue];
                    }

	                    // Dialogue clarity
	                    double clarity = [self calculateDialogueClarity:magnitudes length:fftSize/2 sampleRate:sampleRate];
	                    claritySum += clarity;

                    numFrames++;
                }
                    onsetStartIndexPerSecond.push_back(onsetStart);

	                if (numFrames > 0) {
                    spectralFluxPerSecond.push_back(fluxSum / numFrames);
                    spectralCentroidPerSecond.push_back(centroidSum / numFrames);
                    spectralBalancePerSecond.push_back(balanceSum / numFrames);
                    dialogueClarityPerSecond.push_back(claritySum / numFrames);

                    std::vector<double> avgChroma(12);
                    for (int c = 0; c < 12; c++) {
                        avgChroma[c] = chromaSum[c] / numFrames;
                    }
	                    chromaPerSecond.push_back(avgChroma);
	                } else {
                    spectralFluxPerSecond.push_back(0.0);
                    spectralCentroidPerSecond.push_back(0.0);
                    spectralBalancePerSecond.push_back(0.0);
                    dialogueClarityPerSecond.push_back(0.0);
	                    chromaPerSecond.push_back(std::vector<double>(12, 0.0));
	                }
                    onsetPerSecond.push_back(std::move(secondOnsets));

                    // Stereo width via mid/side energy ratio
                    double midEnergy = 0.0;
                    double sideEnergy = 0.0;
                    for (size_t k = 0; k < n; ++k) {
                        const double mid = 0.5 * ((double)secondBufferL[k] + (double)secondBufferR[k]);
                        const double side = 0.5 * ((double)secondBufferL[k] - (double)secondBufferR[k]);
                        midEnergy += mid * mid;
                        sideEnergy += side * side;
                    }
                    const double stereoWidth = clamp01d(safe_div(sideEnergy, midEnergy + sideEnergy, 0.0));
                    stereoWidthPerSecond.push_back(stereoWidth);

                    // Reverb amount: tail energy ratio within the second (smoothly mapped to [0,1])
                    double reverbAmount = 0.0;
                    if (frameEnergy.size() >= 4) {
                        const size_t F = frameEnergy.size();
                        const size_t headCount = std::max<size_t>(1, F / 3);
                        const size_t tailCount = std::max<size_t>(1, F / 3);
                        double head = 0.0, tail = 0.0;
                        for (size_t fi = 0; fi < headCount; ++fi) head += frameEnergy[fi];
                        for (size_t fi = F - tailCount; fi < F; ++fi) tail += frameEnergy[fi];
                        head = head / (double)headCount;
                        tail = tail / (double)tailCount;
                        const double ratio = safe_div(tail, head + 1e-12, 0.0);
                        reverbAmount = clamp01d(ratio / (1.0 + ratio));
                    }
                    reverbAmountPerSecond.push_back(reverbAmount);

	                // ✅ 只删除已处理的样本，保留剩余部分避免时间轴漂移
	                secondBufferL.erase(secondBufferL.begin(), secondBufferL.begin() + samplesPerSecondPerChannel);
	                secondBufferR.erase(secondBufferR.begin(), secondBufferR.begin() + samplesPerSecondPerChannel);
	                secondBufferMono.clear();
	                currentSecond++;

                if (progressBlock) {
                    progressBlock((double)currentSecond / durationSec);
                }
            }
        }

        CFRelease(sampleBuffer);
    }

    // Local tempo per second (doc 7.2 / 4.2.2): adaptive window from onset dynamics.
    std::vector<double> bpmPerSecond(rmsPerSecond.size(), 0.0);
    std::vector<double> tempoConfPerSecond(rmsPerSecond.size(), 0.0);
    std::vector<double> beatRatioPerSecond(rmsPerSecond.size(), 0.0);
    {
        const double hopRate = sampleRate / hopSize;  // ≈ 86.13 Hz
        const auto z = robust_z_copy(spectralFluxPerSecond);
        const double lambda = mean_abs_slope_1hz(z);
        const double W = adaptive_window_size(lambda);
        const int halfW = std::max(1, (int)std::llround(0.5 * W));

        auto frame_range_for_second = [&](int sec) -> std::pair<size_t,size_t> {
            if (sec < 0) return {0, 0};
            if (sec >= (int)onsetStartIndexPerSecond.size()) return {0, 0};
            const size_t start = onsetStartIndexPerSecond[(size_t)sec];
            const size_t end = (sec + 1 < (int)onsetStartIndexPerSecond.size())
                ? onsetStartIndexPerSecond[(size_t)sec + 1]
                : onsetFrames.size();
            return {start, std::min(end, onsetFrames.size())};
        };

        for (int sec = 0; sec < (int)rmsPerSecond.size(); ++sec) {
            const int s0 = std::max(0, sec - halfW);
            const int s1 = std::min((int)rmsPerSecond.size() - 1, sec + halfW);
            const auto r0 = frame_range_for_second(s0);
            const auto r1 = frame_range_for_second(s1);
            if (r0.second <= r0.first) continue;
            if (r1.second <= r0.first) continue;
            const size_t a = r0.first;
            const size_t b = std::max(r1.second, r0.second);
            if (b <= a || b > onsetFrames.size()) continue;

            std::vector<double> window;
            window.reserve(b - a);
            for (size_t i = a; i < b; ++i) window.push_back(onsetFrames[i]);

            const auto est = estimate_tempo_from_onset_window(window, hopRate, /*bpmMin=*/30.0, /*bpmMax=*/180.0);
            bpmPerSecond[(size_t)sec] = est.bpm;
            tempoConfPerSecond[(size_t)sec] = clamp01d(est.confidence);
            beatRatioPerSecond[(size_t)sec] = clamp01d(est.beatRatio);
        }
    }

    // RhythmDrive: SF energy * beat-neighborhood energy share (doc 4.2.2).
    {
        rhythmEnergyPerSecond.assign(rmsPerSecond.size(), 0.0);
        for (size_t sec = 0; sec < rhythmEnergyPerSecond.size(); ++sec) {
            const double sf = (sec < spectralFluxPerSecond.size()) ? std::max(0.0, spectralFluxPerSecond[sec]) : 0.0;
            const double br = (sec < beatRatioPerSecond.size()) ? beatRatioPerSecond[sec] : 0.0;
            const double tc = (sec < tempoConfPerSecond.size()) ? tempoConfPerSecond[sec] : 0.0;
            rhythmEnergyPerSecond[sec] = clamp01d(sf * br * tc);
        }
    }

    // Loudness per second from BS.1770 gating blocks
    {
        const int seconds = (int)rmsPerSecond.size();
        loudnessPerSecond.assign((size_t)seconds, -120.0);
        if (!blockMs.empty()) {
            std::vector<double> blockL(blockMs.size(), -120.0);
            for (size_t i = 0; i < blockMs.size(); ++i) blockL[i] = lufs_from_mean_square(blockMs[i]);

            // Absolute gate: -70 LUFS
            std::vector<uint8_t> keep(blockMs.size(), 0);
            for (size_t i = 0; i < keep.size(); ++i) keep[i] = (blockL[i] > -70.0) ? 1 : 0;

            auto integrated = [&](const std::vector<uint8_t>& mask) -> double {
                double sum = 0.0;
                double cnt = 0.0;
                for (size_t i = 0; i < blockMs.size(); ++i) {
                    if (!mask[i]) continue;
                    const double ms = blockMs[i];
                    if (!(ms > 0.0)) continue;
                    sum += ms;
                    cnt += 1.0;
                }
                const double msMean = (cnt > 0.0) ? (sum / cnt) : 0.0;
                return lufs_from_mean_square(msMean);
            };

            const double int1 = integrated(keep);
            const double relThr = int1 - 10.0;
            for (size_t i = 0; i < keep.size(); ++i) {
                keep[i] = (keep[i] && blockL[i] > relThr) ? 1 : 0;
            }

            // Aggregate gated blocks into per-second loudness
            std::vector<double> sumMs((size_t)seconds, 0.0);
            std::vector<int> cntMs((size_t)seconds, 0);
            for (size_t i = 0; i < blockMs.size(); ++i) {
                if (!keep[i]) continue;
                const int sec = (i < blockSecondIndex.size()) ? blockSecondIndex[i] : -1;
                if (sec < 0 || sec >= seconds) continue;
                const double ms = blockMs[i];
                if (!(ms > 0.0)) continue;
                sumMs[(size_t)sec] += ms;
                cntMs[(size_t)sec] += 1;
            }
            for (int sec = 0; sec < seconds; ++sec) {
                if (cntMs[(size_t)sec] > 0) {
                    const double msMean = sumMs[(size_t)sec] / (double)cntMs[(size_t)sec];
                    loudnessPerSecond[(size_t)sec] = lufs_from_mean_square(msMean);
                }
            }
        }
    }

    // Audio role separation per second (doc 4.2.6):
    // - DialogueEnergy: speech-band power (300–3400Hz) summed over hops
    // - MusicEnergy: non-speech power gated by tonalness (robust-adaptive, no fixed weights)
    // - EffectsEnergy: residual non-speech power
    std::vector<double> speechProbabilityPerSecond(rmsPerSecond.size(), 0.0);
    {
        const size_t seconds = rmsPerSecond.size();
        dialogueEnergyPerSecond.assign(seconds, 0.0);
        musicEnergyPerSecond.assign(seconds, 0.0);
        effectsEnergyPerSecond.assign(seconds, 0.0);
        dialogueDominancePerSecond.assign(seconds, 0.0);
        musicDominancePerSecond.assign(seconds, 0.0);
        effectsDominancePerSecond.assign(seconds, 0.0);

        if (hopTotalPower.size() == onsetFrames.size() &&
            hopSpeechPower.size() == onsetFrames.size() &&
            hopSpeechRatio.size() == onsetFrames.size() &&
            hopTonalness.size() == onsetFrames.size() &&
            onsetStartIndexPerSecond.size() == seconds) {

            const auto speechScore = robust_sigmoid01_copy(hopSpeechRatio);
            const auto tonalScore = robust_sigmoid01_copy(hopTonalness);

            std::vector<double> totalEnergyPerSecond(seconds, 0.0);

            auto frame_range_for_second = [&](size_t sec) -> std::pair<size_t,size_t> {
                const size_t start = onsetStartIndexPerSecond[sec];
                const size_t end = (sec + 1 < seconds) ? onsetStartIndexPerSecond[sec + 1] : onsetFrames.size();
                return {start, std::min(end, onsetFrames.size())};
            };

            for (size_t sec = 0; sec < seconds; ++sec) {
                const auto r = frame_range_for_second(sec);
                double sumTotal = 0.0;
                double sumDialogue = 0.0;
                double sumMusic = 0.0;
                double sumEffects = 0.0;

                double spW = 0.0;
                double spWSum = 0.0;

                for (size_t i = r.first; i < r.second; ++i) {
                    const double tp = std::max(0.0, hopTotalPower[i]);
                    const double dp = std::clamp(hopSpeechPower[i], 0.0, tp);
                    const double rem = std::max(0.0, tp - dp);
                    const double ts = clamp01d(tonalScore[i]);

                    sumTotal += tp;
                    sumDialogue += dp;
                    sumMusic += rem * ts;
                    sumEffects += rem * (1.0 - ts);

                    spWSum += tp * clamp01d(speechScore[i]);
                    spW += tp;
                }

                totalEnergyPerSecond[sec] = sumTotal;
                dialogueEnergyPerSecond[sec] = sumDialogue;
                musicEnergyPerSecond[sec] = sumMusic;
                effectsEnergyPerSecond[sec] = sumEffects;
                const double den = sumTotal + 1e-12;
                dialogueDominancePerSecond[sec] = clamp01d(sumDialogue / den);
                musicDominancePerSecond[sec] = clamp01d(sumMusic / den);
                effectsDominancePerSecond[sec] = clamp01d(sumEffects / den);
                speechProbabilityPerSecond[sec] = clamp01d(safe_div(spWSum, spW, 0.0));
            }

            // Silence-adaptive cleanup: for very low-energy seconds, mark as "effects" and set speechProbability=0.
            std::vector<double> nz;
            nz.reserve(seconds);
            for (double e : totalEnergyPerSecond) if (e > 0.0 && std::isfinite(e)) nz.push_back(e);
            const double thrSilence = nz.empty() ? 0.0 : quantile_copy(nz, 0.10);
            for (size_t sec = 0; sec < seconds; ++sec) {
                if (!(totalEnergyPerSecond[sec] > thrSilence)) {
                    dialogueEnergyPerSecond[sec] = 0.0;
                    musicEnergyPerSecond[sec] = 0.0;
                    effectsEnergyPerSecond[sec] = 0.0;
                    dialogueDominancePerSecond[sec] = 0.0;
                    musicDominancePerSecond[sec] = 0.0;
                    effectsDominancePerSecond[sec] = 1.0;
                    speechProbabilityPerSecond[sec] = 0.0;
                }
            }
        }
    }

    // Key estimation
    NSMutableArray<NSArray<NSNumber *> *> *chromaArrays = [NSMutableArray array];
    for (const auto& chroma : chromaPerSecond) {
        NSMutableArray<NSNumber *> *chromaArray = [NSMutableArray array];
        for (double c : chroma) {
            [chromaArray addObject:@(c)];
        }
        [chromaArrays addObject:chromaArray];
    }
    int estimatedKey = [self estimateKey:chromaArrays];

	    // Build results
	    for (size_t i = 0; i < rmsPerSecond.size(); i++) {
	        AudioAnalysisResult *result = [[AudioAnalysisResult alloc] init];
	        result.timeSeconds = i + 0.5;
	        result.rms = rmsPerSecond[i];
        result.loudness = i < loudnessPerSecond.size() ? loudnessPerSecond[i] : 0.0;
        result.spectralFlux = i < spectralFluxPerSecond.size() ? spectralFluxPerSecond[i] : 0.0;
        result.spectralCentroid = i < spectralCentroidPerSecond.size() ? spectralCentroidPerSecond[i] : 0.0;
        result.spectralBalance = i < spectralBalancePerSecond.size() ? spectralBalancePerSecond[i] : 0.0;
        result.zeroCrossingRate = i < zcRatePerSecond.size() ? zcRatePerSecond[i] : 0.0;
        result.bpm = i < bpmPerSecond.size() ? bpmPerSecond[i] : 0.0;
        result.tempoConfidence = i < tempoConfPerSecond.size() ? tempoConfPerSecond[i] : 0.0;
        result.rhythmEnergy = i < rhythmEnergyPerSecond.size() ? rhythmEnergyPerSecond[i] : 0.0;

        NSMutableArray<NSNumber *> *chromaArray = [NSMutableArray array];
        if (i < chromaPerSecond.size()) {
            for (double c : chromaPerSecond[i]) {
                [chromaArray addObject:@(c)];
            }
        } else {
            for (int c = 0; c < 12; c++) {
                [chromaArray addObject:@(0.0)];
            }
        }
        result.chroma = chromaArray;
        result.estimatedKey = estimatedKey;
	        result.harmonicTension = [self calculateHarmonicTension:chromaArray key:estimatedKey];
	        result.dialogueClarity = i < dialogueClarityPerSecond.size() ? dialogueClarityPerSecond[i] : 0.0;
            result.speechProbability = i < speechProbabilityPerSecond.size() ? speechProbabilityPerSecond[i] : 0.0;
            result.stereoWidth = i < stereoWidthPerSecond.size() ? stereoWidthPerSecond[i] : 0.0;
            result.reverbAmount = i < reverbAmountPerSecond.size() ? reverbAmountPerSecond[i] : 0.0;
            result.dialogueEnergy = i < dialogueEnergyPerSecond.size() ? dialogueEnergyPerSecond[i] : 0.0;
            result.musicEnergy = i < musicEnergyPerSecond.size() ? musicEnergyPerSecond[i] : 0.0;
            result.effectsEnergy = i < effectsEnergyPerSecond.size() ? effectsEnergyPerSecond[i] : 0.0;
            result.dialogueDominance = i < dialogueDominancePerSecond.size() ? dialogueDominancePerSecond[i] : 0.0;
            result.musicDominance = i < musicDominancePerSecond.size() ? musicDominancePerSecond[i] : 0.0;
            result.effectsDominance = i < effectsDominancePerSecond.size() ? effectsDominancePerSecond[i] : 0.0;

	        [results addObject:result];
	    }
    } // 结束作用域块

cleanup:
    // ✅ 统一的资源清理路径，防止早期 return 造成泄漏
    if (fftSetup) {
        vDSP_destroy_fftsetup(fftSetup);
    }
    if (splitComplex.realp) {
        free(splitComplex.realp);
    }
    if (splitComplex.imagp) {
        free(splitComplex.imagp);
    }
    if (window) {
        free(window);
    }
    if (magnitudes) {
        free(magnitudes);
    }
    if (previousMagnitudes) {
        free(previousMagnitudes);
    }

    return results;
}

#pragma mark - Chroma Extraction

+ (NSArray<NSNumber *> *)extractChroma:(const float *)magnitudes
                                length:(int)length
                            sampleRate:(double)sampleRate {

    double chroma[12] = {0};

    // Frequency resolution
    double freqPerBin = sampleRate / (2.0 * length);

    // Process frequency bins
    for (int k = 1; k < length; k++) {
        double freq = k * freqPerBin;

        // Skip DC and very low frequencies
        if (freq < 65.0) continue;  // Below C2

        // Convert frequency to MIDI note number
        double midiNote = 12.0 * log2(freq / 440.0) + 69.0;

        // Map to chroma bin (0=C, 1=C#, ..., 11=B)
        int chromaBin = ((int)round(midiNote)) % 12;
        if (chromaBin < 0) chromaBin += 12;

        // Accumulate magnitude
        chroma[chromaBin] += magnitudes[k];
    }

    // Normalize
    double sum = 0.0;
    for (int i = 0; i < 12; i++) {
        sum += chroma[i];
    }

    NSMutableArray<NSNumber *> *result = [NSMutableArray array];
    if (sum > 1e-6) {
        for (int i = 0; i < 12; i++) {
            [result addObject:@(chroma[i] / sum)];
        }
    } else {
        for (int i = 0; i < 12; i++) {
            [result addObject:@(1.0/12.0)];
        }
    }

    return result;
}

#pragma mark - Key Estimation

+ (int)estimateKey:(NSArray<NSArray<NSNumber *> *> *)chromaVectors {
    if (chromaVectors.count == 0) return -1;

    // Average chroma over time
    double avgChroma[12] = {0};
    for (NSArray<NSNumber *> *chroma in chromaVectors) {
        for (int i = 0; i < 12; i++) {
            avgChroma[i] += [chroma[i] doubleValue];
        }
    }
    for (int i = 0; i < 12; i++) {
        avgChroma[i] /= chromaVectors.count;
    }

    // Krumhansl-Schmuckler key profiles
    static const double majorProfile[12] = {6.35, 2.23, 3.48, 2.33, 4.38, 4.09, 2.52, 5.19, 2.39, 3.66, 2.29, 2.88};
    static const double minorProfile[12] = {6.33, 2.68, 3.52, 5.38, 2.60, 3.53, 2.54, 4.75, 3.98, 2.69, 3.34, 3.17};

    double bestCorr = -1.0;
    int bestKey = -1;

    // ✅ Try all 24 keys：0-11=major, 12-23=minor
    for (int root = 0; root < 12; root++) {
        // Major
        double corrMajor = 0.0;
        for (int i = 0; i < 12; i++) {
            corrMajor += avgChroma[(root + i) % 12] * majorProfile[i];
        }
        if (corrMajor > bestCorr) {
            bestCorr = corrMajor;
            bestKey = root;  // 0-11: major keys (C=0, C#=1, ..., B=11)
        }

        // Minor
        double corrMinor = 0.0;
        for (int i = 0; i < 12; i++) {
            corrMinor += avgChroma[(root + i) % 12] * minorProfile[i];
        }
        if (corrMinor > bestCorr) {
            bestCorr = corrMinor;
            bestKey = 12 + root;  // ✅ 12-23: minor keys (Cm=12, C#m=13, ..., Bm=23)
        }
    }

    return bestKey;
}

#pragma mark - Harmonic Tension

+ (double)calculateHarmonicTension:(NSArray<NSNumber *> *)chroma key:(int)key {
    if (key < 0 || key >= 24 || chroma.count != 12) {
        return 0.5;  // Neutral
    }

    // ✅ 解码 key：0-11=major, 12-23=minor
    int root = key % 12;
    bool isMinor = (key >= 12);

    // ✅ 根据大小调选择三和弦：major=root+4+7, minor=root+3+7
    int third = isMinor ? (root + 3) % 12 : (root + 4) % 12;
    int fifth = (root + 7) % 12;
    int triad[3] = {root, third, fifth};

    double triadEnergy = 0.0;
    double totalEnergy = 0.0;

    for (int i = 0; i < 12; i++) {
        double energy = [chroma[i] doubleValue];
        totalEnergy += energy;

        if (i == triad[0] || i == triad[1] || i == triad[2]) {
            triadEnergy += energy;
        }
    }

    if (totalEnergy < 1e-6) {
        return 0.5;
    }

    // Tension = 1 - consonance
    double consonance = triadEnergy / totalEnergy;
    return 1.0 - fmin(fmax(consonance, 0.0), 1.0);
}

#pragma mark - BPM Estimation

+ (void)estimateBPM:(const double *)onsetStrength
             length:(int)length
         sampleRate:(double)sampleRate
             outBPM:(double *)outBPM
      outConfidence:(double *)outConfidence {

    if (length < 10) {
        *outBPM = 0.0;
        *outConfidence = 0.0;
        return;
    }

    // Autocorrelation for periodicity detection
    // Search range: 30-180 BPM
    int minLag = (int)(60.0 / 180.0 * sampleRate);  // 180 BPM
    int maxLag = (int)(60.0 / 30.0 * sampleRate);   // 30 BPM

    // ✅ 防止奇怪 sampleRate/hopSize 导致的边界问题
    if (minLag < 1) minLag = 1;

    double maxCorr = 0.0;
    int bestLag = 0;

    for (int lag = minLag; lag < maxLag && lag < length / 2; lag++) {
        double corr = 0.0;
        int count = 0;

        for (int i = 0; i + lag < length; i++) {
            corr += onsetStrength[i] * onsetStrength[i + lag];
            count++;
        }

        if (count > 0) {
            corr /= count;
        }

        if (corr > maxCorr) {
            maxCorr = corr;
            bestLag = lag;
        }
    }

    if (bestLag > 0) {
        *outBPM = 60.0 / (bestLag / sampleRate);

        // Confidence based on autocorrelation strength
        double avgEnergy = 0.0;
        for (int i = 0; i < length; i++) {
            avgEnergy += onsetStrength[i] * onsetStrength[i];
        }
        avgEnergy /= length;

        *outConfidence = (avgEnergy > 1e-6) ? fmin(maxCorr / avgEnergy, 1.0) : 0.0;
    } else {
        *outBPM = 0.0;
        *outConfidence = 0.0;
    }
}

#pragma mark - Dialogue Clarity

+ (double)calculateDialogueClarity:(const float *)magnitudes
                            length:(int)length
                        sampleRate:(double)sampleRate {

    double freqPerBin = sampleRate / (2.0 * length);

    // Speech band: 300-3400 Hz
    int speechStartBin = (int)(300.0 / freqPerBin);
    int speechEndBin = (int)(3400.0 / freqPerBin);

    // ✅ clamp 防止 bin 索引越界
    speechStartBin = std::clamp(speechStartBin, 0, length - 1);
    speechEndBin = std::clamp(speechEndBin, 0, length - 1);

    // ✅ 使用功率（magnitude²）积累能量，更稳定
    double speechPower = 0.0;
    double lowBackgroundPower = 0.0;  // 0-300 Hz
    double highBackgroundPower = 0.0; // 3400+ Hz

    for (int k = 0; k < length; k++) {
        double power = magnitudes[k] * magnitudes[k];

        if (k >= speechStartBin && k <= speechEndBin) {
            speechPower += power;
        } else if (k < speechStartBin) {
            lowBackgroundPower += power;
        } else {
            highBackgroundPower += power;
        }
    }

    // ✅ 总背景能量（低频+高频）
    double backgroundPower = lowBackgroundPower + highBackgroundPower;

    // ✅ 防止除零并计算SNR（dB）
    if (speechPower < 1e-10 && backgroundPower < 1e-10) {
        return 0.0;  // 完全无声
    }
    if (backgroundPower < 1e-10) {
        return 1.0;  // 只有speech频段有能量
    }

    double snr = speechPower / backgroundPower;
    double snrDB = 10.0 * log10(snr + 1e-10);

    // ✅ 把 SNR(dB) 映射到 [0,1]：-6dB->0, +12dB->1
    // SNR < -6dB = 几乎全是背景噪音
    // SNR > +12dB = 语音频段明显占优
    double clarity = (snrDB + 6.0) / 18.0;  // (-6 ~ +12) -> (0 ~ 1)
    return fmin(fmax(clarity, 0.0), 1.0);
}

#pragma mark - Helper Functions

+ (double)calculateRMS:(const float *)samples length:(int)length {
    double sum = 0.0;
    for (int i = 0; i < length; i++) {
        sum += samples[i] * samples[i];
    }
    return sqrt(sum / length);
}

+ (double)calculateZeroCrossingRate:(const float *)samples length:(int)length {
    int crossings = 0;
    for (int i = 1; i < length; i++) {
        if ((samples[i-1] >= 0 && samples[i] < 0) || (samples[i-1] < 0 && samples[i] >= 0)) {
            crossings++;
        }
    }
    return (double)crossings / length;
}

+ (double)calculateSpectralCentroid:(const float *)magnitudes length:(int)length sampleRate:(double)sampleRate {
    double freqPerBin = sampleRate / (2.0 * length);

    double weightedSum = 0.0;
    double totalMag = 0.0;

    for (int k = 0; k < length; k++) {
        double freq = k * freqPerBin;
        const double mag = magnitudes[k];
        const double pwr = mag * mag;
        weightedSum += freq * pwr;
        totalMag += pwr;
    }

    return (totalMag > 1e-6) ? (weightedSum / totalMag) : 0.0;
}

+ (double)calculateSpectralBalance:(const float *)magnitudes length:(int)length {
    double lowEnergy = 0.0;
    double highEnergy = 0.0;

    int midpoint = length / 2;

    for (int k = 0; k < length; k++) {
        const double mag = magnitudes[k];
        const double pwr = mag * mag;
        if (k < midpoint) {
            lowEnergy += pwr;
        } else {
            highEnergy += pwr;
        }
    }

    if (lowEnergy + highEnergy < 1e-6) {
        return 0.5;
    }

    return highEnergy / (lowEnergy + highEnergy);
}

@end
