#import "VideoFeatureExtractor.h"
#import "FaceTrackingEngine.h"
#import "AudioAnalysisEngine.h"
#import "SubtitleNLP.h"
#import "SubtitleCue.h"

#import <Accelerate/Accelerate.h>
#import <AVFoundation/AVFoundation.h>
#import <CoreImage/CoreImage.h>
#import <Vision/Vision.h>
#import <CoreML/CoreML.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

// Import shared SubtitleCue from bridge namespace
using arcscope::bridge::SubtitleCue;

namespace {

struct FrameStats {
    double brightness{0.5};
    double saturation{0.5};
    double warmth{0.5};
    double grayness{0.5};
    double avgHue{-1.0};    // Average hue angle [0,360), -1 if invalid/gray
};

struct Cam16SecondStats {
    double J{0.0};
    double M{0.0};
    double h{-1.0};
    double Jp{0.0};
    double ap{0.0};
    double bp{0.0};
    double deltaE{0.0};
    std::array<double, 12> hueHist{};
    double harmony{0.0};
};

constexpr double kAnalysisFrameRate = 10.0;  // arcscope.md default f_analysis
constexpr size_t kAnalysisTargetHeight = 720; // arcscope.md default downsample height

static inline double clamp01(double v) {
    if (v < 0.0) return 0.0;
    if (v > 1.0) return 1.0;
    return v;
}

constexpr double kPi = 3.14159265358979323846;

double quantile_copy(std::vector<double> v, double q) {
    if (v.empty()) return 0.0;
    q = std::clamp(q, 0.0, 1.0);
    const size_t n = v.size();
    const size_t k = static_cast<size_t>(std::llround(q * static_cast<double>(n - 1)));
    std::nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
}

static inline double rad2deg(double r) { return r * 180.0 / kPi; }
static inline double deg2rad(double d) { return d * kPi / 180.0; }
static inline double wrap_deg_360(double d) {
    double r = std::fmod(d, 360.0);
    if (r < 0.0) r += 360.0;
    return r;
}

// ===== CAM16-UCS implementation (Rec.2020 linear -> XYZ -> CAM16 -> CAM16-UCS) =====
// NOTE: arcscope.md allows ACEScg or an equivalent Rec.2020 linear baseline.
struct Cam16ViewingConditions {
    // D65 reference white in XYZ with Y=100
    double Xw{95.047};
    double Yw{100.0};
    double Zw{108.883};

    // Standard "average surround" viewing conditions (not tunable model weights)
    double LA{64.0}; // adapting luminance (cd/m^2)
    double Yb{20.0}; // background relative luminance
    double F{1.0};
    double c{0.69};
    double Nc{1.0};
};

struct Cam16Ucs {
    double J{0.0};
    double M{0.0};
    double h{-1.0};
    double Jp{0.0};
    double ap{0.0};
    double bp{0.0};
};

static inline std::array<double,3> rec2020_linear_to_xyz(double r, double g, double b) {
    // Rec.2020 (D65), linear RGB -> XYZ, scaled so Y=1.0.
    constexpr double M[3][3] = {
        {0.6369580483012914, 0.1446169035862083, 0.1688809751641721},
        {0.2627002120112671, 0.6779980715188708, 0.05930171646986196},
        {0.0000000000000000, 0.02807269304908743, 1.060985057710791}
    };
    const double X = M[0][0]*r + M[0][1]*g + M[0][2]*b;
    const double Y = M[1][0]*r + M[1][1]*g + M[1][2]*b;
    const double Z = M[2][0]*r + M[2][1]*g + M[2][2]*b;
    return {X * 100.0, Y * 100.0, Z * 100.0};
}

static inline Cam16Ucs xyz_to_cam16ucs(double X, double Y, double Z, const Cam16ViewingConditions& vc) {
    constexpr double M_CAT16[3][3] = {
        { 0.401288,  0.650173, -0.051461},
        {-0.250268,  1.204414,  0.045854},
        {-0.002079,  0.048952,  0.953127}
    };
    auto mul3 = [&](double x, double y, double z) -> std::array<double,3> {
        return {
            M_CAT16[0][0]*x + M_CAT16[0][1]*y + M_CAT16[0][2]*z,
            M_CAT16[1][0]*x + M_CAT16[1][1]*y + M_CAT16[1][2]*z,
            M_CAT16[2][0]*x + M_CAT16[2][1]*y + M_CAT16[2][2]*z
        };
    };

    const auto RGBw = mul3(vc.Xw, vc.Yw, vc.Zw);
    const auto RGB  = mul3(X, Y, Z);

    const double k = 1.0 / (5.0 * vc.LA + 1.0);
    const double k4 = k*k*k*k;
    const double Fl = 0.2 * k4 * 5.0 * vc.LA + 0.1 * (1.0 - k4) * (1.0 - k4) * std::cbrt(5.0 * vc.LA);
    const double n = vc.Yb / vc.Yw;
    const double z = 1.48 + std::sqrt(n);
    const double Nbb = 0.725 * std::pow(1.0 / n, 0.2);
    const double Ncb = Nbb;

    const double D = vc.F * (1.0 - (1.0/3.6) * std::exp((-vc.LA - 42.0)/92.0));
    auto adapt = [&](double c, double cw) -> double {
        if (!(cw > 1e-12)) return c;
        return (D * (vc.Yw / cw) + 1.0 - D) * c;
    };

    const double Rc = adapt(RGB[0], RGBw[0]);
    const double Gc = adapt(RGB[1], RGBw[1]);
    const double Bc = adapt(RGB[2], RGBw[2]);

    auto f_nl = [&](double t) -> double {
        const double p = std::pow((Fl * std::abs(t) / 100.0), 0.42);
        const double s = (t < 0.0) ? -1.0 : 1.0;
        return s * (400.0 * p / (p + 27.13)) + 0.1;
    };

    const double R_a = f_nl(Rc);
    const double G_a = f_nl(Gc);
    const double B_a = f_nl(Bc);

    const double Rw_a = f_nl(adapt(RGBw[0], RGBw[0]));
    const double Gw_a = f_nl(adapt(RGBw[1], RGBw[1]));
    const double Bw_a = f_nl(adapt(RGBw[2], RGBw[2]));

    const double a = R_a - 12.0 * G_a / 11.0 + B_a / 11.0;
    const double b = (R_a + G_a - 2.0 * B_a) / 9.0;

    const double h = wrap_deg_360(rad2deg(std::atan2(b, a)));
    const double A = (2.0 * R_a + G_a + 0.05 * B_a - 0.305) * Nbb;
    const double Aw = (2.0 * Rw_a + Gw_a + 0.05 * Bw_a - 0.305) * Nbb;

    const double J = 100.0 * std::pow(std::max(0.0, A / std::max(Aw, 1e-12)), vc.c * z);

    const double et = 0.25 * (std::cos(deg2rad(h) + 2.0) + 3.8);
    const double t = (50000.0/13.0) * vc.Nc * Ncb * et * std::sqrt(a*a + b*b) /
                     (R_a + G_a + (21.0/20.0) * B_a);
    const double alpha = std::pow(t, 0.9) * std::pow(1.64 - std::pow(0.29, n), 0.73);
    const double C = alpha * std::sqrt(J / 100.0);
    const double M = C * std::pow(Fl, 0.25);

    constexpr double c1 = 0.007;
    constexpr double c2 = 0.0228;
    const double Jp = (1.0 + 100.0 * c1) * J / (1.0 + c1 * J);
    const double Mp = (1.0 / c2) * std::log(1.0 + c2 * M);
    const double hr = deg2rad(h);
    const double ap = Mp * std::cos(hr);
    const double bp = Mp * std::sin(hr);

    Cam16Ucs out;
    out.J = J;
    out.M = M;
    out.h = h;
    out.Jp = Jp;
    out.ap = ap;
    out.bp = bp;
    return out;
}

double ParseTimestamp(NSString* value) {
    NSArray<NSString*>* components = [value componentsSeparatedByCharactersInSet:[NSCharacterSet characterSetWithCharactersInString:@":,"]];
    if (components.count < 4) {
        return 0.0;
    }
    double hours = components[0].doubleValue;
    double minutes = components[1].doubleValue;
    double seconds = components[2].doubleValue;
    double millis = components[3].doubleValue;
    return hours * 3600.0 + minutes * 60.0 + seconds + millis / 1000.0;
}

/**
 * Language-aware word count estimation
 * Handles both Western (space-delimited) and CJK (character-based) languages
 */
int EstimateWordCount(NSString* line) {
    if (!line || line.length == 0) {
        return 0;
    }

    // Count spaces to detect language type
    NSUInteger spaceCount = 0;
    for (NSUInteger i = 0; i < line.length; ++i) {
        if ([line characterAtIndex:i] == ' ') {
            spaceCount++;
        }
    }

    // If very few spaces, assume CJK (Chinese/Japanese/Korean)
    if (spaceCount <= 1) {
        __block int charCount = 0;
        [line enumerateSubstringsInRange:NSMakeRange(0, line.length)
                                 options:NSStringEnumerationByComposedCharacterSequences
                              usingBlock:^(NSString* substring, NSRange range, NSRange enclosingRange, BOOL* stop) {
            unichar ch = [substring characterAtIndex:0];
            // Count alphanumeric and CJK ideographs
            if ([[NSCharacterSet alphanumericCharacterSet] characterIsMember:ch]) {
                charCount++;
            } else if (ch >= 0x4E00 && ch <= 0x9FFF) {
                // CJK Unified Ideographs
                charCount++;
            } else if (ch >= 0x3400 && ch <= 0x4DBF) {
                // CJK Extension A
                charCount++;
            } else if (ch >= 0xAC00 && ch <= 0xD7AF) {
                // Hangul syllables
                charCount++;
            } else if (ch >= 0x3040 && ch <= 0x30FF) {
                // Hiragana/Katakana
                charCount++;
            }
        }];
        return charCount;
    }

    // Western language: count space-separated words
    NSArray<NSString*>* parts = [line componentsSeparatedByCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
    int wordCount = 0;
    for (NSString* word in parts) {
        if (word.length > 0) {
            wordCount++;
        }
    }
    return wordCount;
}

std::vector<SubtitleCue> LoadSubtitleCues(NSURL* url) {
    std::vector<SubtitleCue> cues;
    NSURL* srtURL = [[url URLByDeletingPathExtension] URLByAppendingPathExtension:@"srt"];
    if (!srtURL) {
        return cues;
    }
    NSError* error = nil;
    NSString* raw = [NSString stringWithContentsOfURL:srtURL encoding:NSUTF8StringEncoding error:&error];
    if (!raw) {
        return cues;
    }
    NSString* normalized = [[raw stringByReplacingOccurrencesOfString:@"\r\n" withString:@"\n"] stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
    NSArray<NSString*>* blocks = [normalized componentsSeparatedByString:@"\n\n"];
    for (NSString* block in blocks) {
        NSArray<NSString*>* lines = [block componentsSeparatedByString:@"\n"];
        if (lines.count < 2) {
            continue;
        }
        NSUInteger index = 0;
        NSRange arrowRange = [lines[0] rangeOfString:@"-->"];
        if ([[NSCharacterSet decimalDigitCharacterSet] characterIsMember:[lines[0] characterAtIndex:0]] && arrowRange.location == NSNotFound) {
            index = 1;
        }
        if (index >= lines.count) {
            continue;
        }
        NSArray<NSString*>* timeParts = [lines[index] componentsSeparatedByString:@" --> " ];
        if (timeParts.count != 2) {
            continue;
        }
        double start = ParseTimestamp(timeParts[0]);
        double end = ParseTimestamp(timeParts[1]);
        if (end <= start) {
            continue;
        }

        // Collect subtitle text lines
        NSMutableString* textBuilder = [NSMutableString string];
        for (NSUInteger lineIdx = index + 1; lineIdx < lines.count; ++lineIdx) {
            NSString* line = lines[lineIdx];
            if (line.length > 0) {
                if (textBuilder.length > 0) {
                    [textBuilder appendString:@" "];
                }
                [textBuilder appendString:line];
            }
        }

        // Use language-aware word counting
        NSString* fullText = [textBuilder copy];
        int wordCount = EstimateWordCount(fullText);

        SubtitleCue cue;
        cue.start = start;
        cue.end = end;
        cue.wordCount = wordCount;
        cue.isExposition = wordCount >= 12;
        // Store as UTF-8 string (safe for STL containers, avoids ARC issues)
        cue.text_utf8 = fullText ? std::string([fullText UTF8String]) : std::string();
        cues.push_back(cue);
    }
    return cues;
}

void ComputeSubtitleGrids(const std::vector<SubtitleCue>& cues,
                          size_t sampleCount,
                          std::vector<double>& wordDensity,
                          std::vector<double>& expositionDensity) {
    wordDensity.assign(sampleCount, 0.0);
    expositionDensity.assign(sampleCount, 0.0);
    for (const auto& cue : cues) {
        if (cue.end <= cue.start) {
            continue;
        }
        const double duration = std::max(0.25, cue.end - cue.start);
        const int startSec = std::max(0, static_cast<int>(std::floor(cue.start)));
        const int endSec = std::min(static_cast<int>(sampleCount), static_cast<int>(std::ceil(cue.end)));
        for (int sec = startSec; sec < endSec; ++sec) {
            double secStart = static_cast<double>(sec);
            double secEnd = secStart + 1.0;
            double overlap = std::max(0.0, std::min(cue.end, secEnd) - std::max(cue.start, secStart));
            if (overlap <= 0.0) {
                continue;
            }
            double ratio = overlap / duration;
            if (sec < wordDensity.size()) {
                wordDensity[sec] += cue.wordCount * ratio;
                if (cue.isExposition) {
                    expositionDensity[sec] += ratio;
                }
            }
        }
    }
}

CVPixelBufferRef create_argb_pixel_buffer(size_t width, size_t height) {
    CVPixelBufferRef pb = nullptr;
    const OSType fmt = kCVPixelFormatType_32ARGB;
    NSDictionary* attrs = @{
        (id)kCVPixelBufferCGImageCompatibilityKey: @YES,
        (id)kCVPixelBufferCGBitmapContextCompatibilityKey: @YES,
        (id)kCVPixelBufferIOSurfacePropertiesKey: @{}
    };
    CVReturn rc = CVPixelBufferCreate(kCFAllocatorDefault, (int)width, (int)height, fmt,
                                     (__bridge CFDictionaryRef)attrs, &pb);
    return (rc == kCVReturnSuccess) ? pb : nullptr;
}

// Apply preferredTransform and scale to fixed target height (720p) for analysis.
CVPixelBufferRef create_oriented_scaled_720p(CVPixelBufferRef src,
                                             CGAffineTransform preferredTransform,
                                             CIContext* ciContext,
                                             CGColorSpaceRef colorSpace) {
    if (!src || !ciContext || !colorSpace) return nullptr;

    CIImage* img = [CIImage imageWithCVPixelBuffer:src];
    img = [img imageByApplyingTransform:preferredTransform];
    CGRect extent = img.extent;
    img = [img imageByCroppingToRect:extent];
    img = [img imageByApplyingTransform:CGAffineTransformMakeTranslation(-extent.origin.x, -extent.origin.y)];

    const double srcH = extent.size.height;
    if (!(srcH > 1.0)) return nullptr;
    const double scale = static_cast<double>(kAnalysisTargetHeight) / srcH;
    const long long dstWll = std::max(1LL, std::llround(extent.size.width * scale));
    const size_t dstW = static_cast<size_t>(dstWll);
    const size_t dstH = static_cast<size_t>(kAnalysisTargetHeight);

    CIImage* scaled = [img imageByApplyingTransform:CGAffineTransformMakeScale(scale, scale)];
    CVPixelBufferRef dst = create_argb_pixel_buffer(dstW, dstH);
    if (!dst) return nullptr;

    [ciContext render:scaled toCVPixelBuffer:dst bounds:CGRectMake(0, 0, (CGFloat)dstW, (CGFloat)dstH) colorSpace:colorSpace];
    return dst;
}

FrameStats compute_frame_stats_from_argb(CVPixelBufferRef pb) {
    FrameStats stats;
    if (!pb) return stats;

    CVPixelBufferLockBaseAddress(pb, kCVPixelBufferLock_ReadOnly);
    vImage_Buffer src;
    src.data = CVPixelBufferGetBaseAddress(pb);
    src.width = CVPixelBufferGetWidth(pb);
    src.height = CVPixelBufferGetHeight(pb);
    src.rowBytes = CVPixelBufferGetBytesPerRow(pb);

    // Channel histograms (A,R,G,B) for fast means
    vImagePixelCount histA[256] = {0};
    vImagePixelCount histR[256] = {0};
    vImagePixelCount histG[256] = {0};
    vImagePixelCount histB[256] = {0};
    vImagePixelCount* hists[4] = {histA, histR, histG, histB};
    (void)vImageHistogramCalculation_ARGB8888(&src, hists, kvImageNoFlags);

    const double N = static_cast<double>(src.width) * static_cast<double>(src.height);
    auto mean_from_hist = [&](vImagePixelCount* h) -> double {
        double acc = 0.0;
        for (int i = 0; i < 256; ++i) acc += static_cast<double>(h[i]) * static_cast<double>(i);
        return (N > 0.0) ? (acc / N) / 255.0 : 0.0;
    };

    const double avgR = mean_from_hist(histR);
    const double avgG = mean_from_hist(histG);
    const double avgB = mean_from_hist(histB);

    stats.brightness = clamp01((avgR + avgG + avgB) / 3.0);
    const double maxRGB = std::max({avgR, avgG, avgB});
    const double minRGB = std::min({avgR, avgG, avgB});
    stats.saturation = clamp01(maxRGB - minRGB);
    stats.grayness = clamp01(1.0 - stats.saturation);
    stats.warmth = clamp01((avgR - avgB + 1.0) * 0.5);

    // Avg hue from mean RGB (fallback, per-pixel hue distribution is computed in Core ColorAnalyzer)
    const double delta = maxRGB - minRGB;
    if (delta > 1e-6 && stats.saturation > 0.05) {
        double hue = 0.0;
        if (maxRGB == avgR) hue = 60.0 * std::fmod(((avgG - avgB) / delta), 6.0);
        else if (maxRGB == avgG) hue = 60.0 * (((avgB - avgR) / delta) + 2.0);
        else hue = 60.0 * (((avgR - avgG) / delta) + 4.0);
        if (hue < 0.0) hue += 360.0;
        stats.avgHue = hue;
    } else {
        stats.avgHue = -1.0;
    }

    CVPixelBufferUnlockBaseAddress(pb, kCVPixelBufferLock_ReadOnly);
    return stats;
}

// Mean absolute luma difference between two ARGB8888 frames, normalized to [0,1]
double mean_abs_luma_diff_argb(CVPixelBufferRef a, CVPixelBufferRef b,
                               std::vector<uint8_t>& lumaA,
                               std::vector<uint8_t>& lumaB,
                               std::vector<uint8_t>& lumaDiff) {
    if (!a || !b) return 0.0;
    const size_t w = CVPixelBufferGetWidth(a);
    const size_t h = CVPixelBufferGetHeight(a);
    if (w != CVPixelBufferGetWidth(b) || h != CVPixelBufferGetHeight(b)) return 0.0;

    const size_t count = w * h;
    lumaA.resize(count);
    lumaB.resize(count);
    lumaDiff.resize(count);

    CVPixelBufferLockBaseAddress(a, kCVPixelBufferLock_ReadOnly);
    CVPixelBufferLockBaseAddress(b, kCVPixelBufferLock_ReadOnly);

    vImage_Buffer srcA{CVPixelBufferGetBaseAddress(a), h, w, CVPixelBufferGetBytesPerRow(a)};
    vImage_Buffer srcB{CVPixelBufferGetBaseAddress(b), h, w, CVPixelBufferGetBytesPerRow(b)};
    vImage_Buffer dstA{lumaA.data(), h, w, w};
    vImage_Buffer dstB{lumaB.data(), h, w, w};

    // ITU-R BT.709 luma coefficients (scaled by 256): Y = 0.2126R + 0.7152G + 0.0722B
    // ARGB order in memory: A R G B
    const int16_t matrix[4] = {0, 54, 183, 19};
    const int16_t preBias[4] = {0, 0, 0, 0};
    const int32_t postBias = 0;
    (void)vImageMatrixMultiply_ARGB8888ToPlanar8(&srcA, &dstA, matrix, 256, preBias, postBias, kvImageNoFlags);
    (void)vImageMatrixMultiply_ARGB8888ToPlanar8(&srcB, &dstB, matrix, 256, preBias, postBias, kvImageNoFlags);

    // Compute mean absolute difference in float space (exact, vectorized)
    std::vector<float> fa(count);
    std::vector<float> fb(count);
    std::vector<float> fd(count);

    vDSP_vfltu8(lumaA.data(), 1, fa.data(), 1, count);
    vDSP_vfltu8(lumaB.data(), 1, fb.data(), 1, count);
    vDSP_vsub(fb.data(), 1, fa.data(), 1, fd.data(), 1, count); // fd = fa - fb
    vDSP_vabs(fd.data(), 1, fd.data(), 1, count);
    float sum = 0.0f;
    vDSP_sve(fd.data(), 1, &sum, count);
    const double meanAbs = (count > 0) ? (static_cast<double>(sum) / static_cast<double>(count)) : 0.0;

    CVPixelBufferUnlockBaseAddress(b, kCVPixelBufferLock_ReadOnly);
    CVPixelBufferUnlockBaseAddress(a, kCVPixelBufferLock_ReadOnly);

    return clamp01(meanAbs / 255.0);
}

struct OpticalFlowStats {
    double meanU{0.0};
    double meanV{0.0};
    double meanSq{0.0};   // E[u^2 + v^2]
    double zoom{0.0};     // radial component proxy (positive=expansion)
};

// Compute simple optical-flow statistics (global translation + residual + zoom proxy).
// All values are in raw flow units (pixels per frame); normalize by diag in caller.
OpticalFlowStats optical_flow_stats(VNSequenceRequestHandler* handler,
                                    CVPixelBufferRef prev,
                                    CVPixelBufferRef cur) {
    OpticalFlowStats out;
    if (!handler || !prev || !cur) return out;

    VNGenerateOpticalFlowRequest* req = [[VNGenerateOpticalFlowRequest alloc] initWithTargetedCVPixelBuffer:cur options:@{}];
    req.computationAccuracy = VNGenerateOpticalFlowRequestComputationAccuracyHigh;
    req.outputPixelFormat = kCVPixelFormatType_TwoComponent32Float;

    NSError* err = nil;
    BOOL ok = [handler performRequests:@[req] onCVPixelBuffer:prev error:&err];
    if (!ok || err) return out;

    VNPixelBufferObservation* obs = (VNPixelBufferObservation*)req.results.firstObject;
    if (!obs) return out;
    CVPixelBufferRef flow = obs.pixelBuffer;
    if (!flow) return out;

    CVPixelBufferLockBaseAddress(flow, kCVPixelBufferLock_ReadOnly);
    const size_t w = CVPixelBufferGetWidth(flow);
    const size_t h = CVPixelBufferGetHeight(flow);
    const size_t rowBytes = CVPixelBufferGetBytesPerRow(flow);
    const float* base = (const float*)CVPixelBufferGetBaseAddress(flow);
    if (!base) {
        CVPixelBufferUnlockBaseAddress(flow, kCVPixelBufferLock_ReadOnly);
        return out;
    }

    double sumU = 0.0;
    double sumV = 0.0;
    double sumUU = 0.0;
    double sumVV = 0.0;
    double sumRadial = 0.0;
    double sumRadNorm = 0.0;
    double n = 0.0;

    const double cx = 0.5 * (double)(w - 1);
    const double cy = 0.5 * (double)(h - 1);
    for (size_t y = 0; y < h; ++y) {
        const float* row = (const float*)((const uint8_t*)base + y * rowBytes);
        for (size_t x = 0; x < w; ++x) {
            const float u = row[2 * x + 0];
            const float v = row[2 * x + 1];
            sumU += (double)u;
            sumV += (double)v;
            sumUU += (double)u * (double)u;
            sumVV += (double)v * (double)v;

            const double dx = (double)x - cx;
            const double dy = (double)y - cy;
            const double r2 = dx * dx + dy * dy;
            if (r2 > 1e-6) {
                sumRadial += (double)u * dx + (double)v * dy;
                sumRadNorm += r2;
            }
            n += 1.0;
        }
    }
    CVPixelBufferUnlockBaseAddress(flow, kCVPixelBufferLock_ReadOnly);

    if (n <= 1e-9) return out;
    out.meanU = sumU / n;
    out.meanV = sumV / n;
    out.meanSq = (sumUU + sumVV) / n;
    out.zoom = (sumRadNorm > 1e-9) ? (sumRadial / sumRadNorm) : 0.0;
    return out;
}

// Mean optical flow magnitude between two frames, normalized by image diagonal (legacy helper).
double mean_optical_flow_mag(VNSequenceRequestHandler* handler,
                             CVPixelBufferRef prev,
                             CVPixelBufferRef cur) {
    const auto st = optical_flow_stats(handler, prev, cur);
    if (!prev) return 0.0;
    const size_t w = CVPixelBufferGetWidth(prev);
    const size_t h = CVPixelBufferGetHeight(prev);
    const double diag = std::hypot((double)w, (double)h);
    if (!(diag > 1e-12)) return 0.0;
    const double meanSq = std::max(0.0, st.meanSq);
    return std::sqrt(meanSq) / diag;
}

static inline float half_to_float(uint16_t hv) {
    uint32_t sign = (hv >> 15) & 1u;
    uint32_t exp  = (hv >> 10) & 0x1Fu;
    uint32_t mant = hv & 0x3FFu;
    uint32_t f;
    if (exp == 0) {
        if (mant == 0) {
            f = sign << 31;
        } else {
            exp = 127 - 15 + 1;
            while ((mant & 0x400u) == 0) { mant <<= 1; exp--; }
            mant &= 0x3FFu;
            f = (sign << 31) | (exp << 23) | (mant << 13);
        }
    } else if (exp == 31) {
        f = (sign << 31) | 0x7F800000u | (mant << 13);
    } else {
        exp = exp + (127 - 15);
        f = (sign << 31) | (exp << 23) | (mant << 13);
    }
    float outv;
    std::memcpy(&outv, &f, sizeof(float));
    return outv;
}

CVPixelBufferRef create_rgba_half_pixel_buffer(size_t width, size_t height) {
    CVPixelBufferRef pb = nullptr;
    const OSType fmt = kCVPixelFormatType_64RGBAHalf;
    NSDictionary* attrs = @{
        (id)kCVPixelBufferIOSurfacePropertiesKey: @{},
        (id)kCVPixelBufferMetalCompatibilityKey: @YES
    };
    CVReturn rc = CVPixelBufferCreate(kCFAllocatorDefault, (int)width, (int)height, fmt,
                                     (__bridge CFDictionaryRef)attrs, &pb);
    return (rc == kCVReturnSuccess) ? pb : nullptr;
}

CVPixelBufferRef render_to_linear_rec2020_half(CVPixelBufferRef argb,
                                               CIContext* ciContext,
                                               CGColorSpaceRef linearRec2020) {
    if (!argb || !ciContext || !linearRec2020) return nullptr;
    const size_t w = CVPixelBufferGetWidth(argb);
    const size_t h = CVPixelBufferGetHeight(argb);
    CVPixelBufferRef dst = create_rgba_half_pixel_buffer(w, h);
    if (!dst) return nullptr;
    CIImage* img = [CIImage imageWithCVPixelBuffer:argb];
    [ciContext render:img toCVPixelBuffer:dst bounds:CGRectMake(0, 0, (CGFloat)w, (CGFloat)h) colorSpace:linearRec2020];
    return dst;
}

Cam16SecondStats cam16_stats_from_linear_rec2020_half(CVPixelBufferRef rgbaHalf,
                                                      const Cam16ViewingConditions& vc) {
    Cam16SecondStats out;
    out.hueHist.fill(0.0);
    if (!rgbaHalf) return out;

    CVPixelBufferLockBaseAddress(rgbaHalf, kCVPixelBufferLock_ReadOnly);
    const size_t w = CVPixelBufferGetWidth(rgbaHalf);
    const size_t h = CVPixelBufferGetHeight(rgbaHalf);
    const size_t rowBytes = CVPixelBufferGetBytesPerRow(rgbaHalf);
    const uint16_t* base = (const uint16_t*)CVPixelBufferGetBaseAddress(rgbaHalf);
    if (!base) {
        CVPixelBufferUnlockBaseAddress(rgbaHalf, kCVPixelBufferLock_ReadOnly);
        return out;
    }

    double sumJ = 0.0, sumM = 0.0;
    double sumJp = 0.0, sumAp = 0.0, sumBp = 0.0;
    double n = 0.0;

    for (size_t y = 0; y < h; ++y) {
        const uint16_t* row = (const uint16_t*)((const uint8_t*)base + y * rowBytes);
        for (size_t x = 0; x < w; ++x) {
            const float r = half_to_float(row[4 * x + 0]);
            const float g = half_to_float(row[4 * x + 1]);
            const float b = half_to_float(row[4 * x + 2]);
            if (!(std::isfinite(r) && std::isfinite(g) && std::isfinite(b))) continue;

            const auto xyz = rec2020_linear_to_xyz(std::max(0.0f, r), std::max(0.0f, g), std::max(0.0f, b));
            const auto cam = xyz_to_cam16ucs(xyz[0], xyz[1], xyz[2], vc);

            sumJ += cam.J;
            sumM += cam.M;
            sumJp += cam.Jp;
            sumAp += cam.ap;
            sumBp += cam.bp;
            n += 1.0;

            const int bin = std::clamp((int)std::floor(cam.h / 30.0), 0, 11);
            out.hueHist[(size_t)bin] += std::hypot(cam.ap, cam.bp); // weight by Mp (UCS chroma)
        }
    }

    CVPixelBufferUnlockBaseAddress(rgbaHalf, kCVPixelBufferLock_ReadOnly);

    if (!(n > 0.0)) return out;

    out.J = sumJ / n;
    out.M = sumM / n;
    out.Jp = sumJp / n;
    out.ap = sumAp / n;
    out.bp = sumBp / n;
    out.h = wrap_deg_360(rad2deg(std::atan2(out.bp, out.ap)));

    double hs = 0.0;
    for (double v : out.hueHist) hs += v;
    if (hs > 1e-12) for (double& v : out.hueHist) v /= hs;
    else out.hueHist.fill(0.0);

    return out;
}

double cam16_harmony_kmedians(const std::vector<std::array<float,3>>& pts, int K) {
    if (pts.empty() || K <= 0) return 0.0;
    K = std::min<int>(K, (int)pts.size());
    if (K <= 1) return 0.0;

    auto l2 = [](const std::array<float,3>& a, const std::array<float,3>& b) -> double {
        const double dx = (double)a[0] - (double)b[0];
        const double dy = (double)a[1] - (double)b[1];
        const double dz = (double)a[2] - (double)b[2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };
    auto l1 = [](const std::array<float,3>& a, const std::array<float,3>& b) -> double {
        return std::abs((double)a[0]-b[0]) + std::abs((double)a[1]-b[1]) + std::abs((double)a[2]-b[2]);
    };

    std::vector<std::array<float,3>> centers;
    centers.reserve(K);
    // seed 0: median in J'
    {
        std::vector<float> Js;
        Js.reserve(pts.size());
        for (const auto& p : pts) Js.push_back(p[0]);
        const size_t mid = Js.size() / 2;
        std::nth_element(Js.begin(), Js.begin() + mid, Js.end());
        const float Jmed = Js[mid];
        size_t bestIdx = 0;
        double bestD = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < pts.size(); ++i) {
            const double d = std::abs((double)pts[i][0] - (double)Jmed);
            if (d < bestD) { bestD = d; bestIdx = i; }
        }
        centers.push_back(pts[bestIdx]);
    }
    // farthest-point seeding for the rest
    while ((int)centers.size() < K) {
        size_t bestIdx = 0;
        double bestMin = -1.0;
        for (size_t i = 0; i < pts.size(); ++i) {
            double mn = std::numeric_limits<double>::infinity();
            for (const auto& c : centers) mn = std::min(mn, l2(pts[i], c));
            if (mn > bestMin) { bestMin = mn; bestIdx = i; }
        }
        centers.push_back(pts[bestIdx]);
    }

    std::vector<int> assign(pts.size(), 0);
    for (int iter = 0; iter < 20; ++iter) {
        bool changed = false;
        for (size_t i = 0; i < pts.size(); ++i) {
            int bestK = 0;
            double bestD = std::numeric_limits<double>::infinity();
            for (int k = 0; k < K; ++k) {
                const double d = l1(pts[i], centers[(size_t)k]);
                if (d < bestD) { bestD = d; bestK = k; }
            }
            if (assign[i] != bestK) { assign[i] = bestK; changed = true; }
        }

        for (int k = 0; k < K; ++k) {
            std::vector<float> xs, ys, zs;
            xs.reserve(pts.size()/K + 1);
            ys.reserve(pts.size()/K + 1);
            zs.reserve(pts.size()/K + 1);
            for (size_t i = 0; i < pts.size(); ++i) {
                if (assign[i] == k) {
                    xs.push_back(pts[i][0]);
                    ys.push_back(pts[i][1]);
                    zs.push_back(pts[i][2]);
                }
            }
            if (xs.empty()) continue;
            const size_t m = xs.size()/2;
            std::nth_element(xs.begin(), xs.begin() + m, xs.end());
            std::nth_element(ys.begin(), ys.begin() + m, ys.end());
            std::nth_element(zs.begin(), zs.begin() + m, zs.end());
            centers[(size_t)k] = {xs[m], ys[m], zs[m]};
        }

        if (!changed) break;
    }

    double within = 0.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        within += l2(pts[i], centers[(size_t)assign[i]]);
    }
    within /= std::max<size_t>(1, pts.size());

    double between = 0.0;
    int pairs = 0;
    for (int i = 0; i < K; ++i) {
        for (int j = i + 1; j < K; ++j) {
            between += l2(centers[(size_t)i], centers[(size_t)j]);
            pairs++;
        }
    }
    if (pairs > 0) between /= (double)pairs;

    const double denom = between + within;
    if (!(denom > 1e-12)) return 0.0;
    return std::clamp(between / denom, 0.0, 1.0);
}

}  // namespace

bool ExtractSamplesFromVideo(NSURL* url,
                             std::vector<ArcScopeBridgeSample>& outSamples,
                             NSArray** outTracks,
                             std::vector<double>* outCutTimesSec) {
    if (!url) {
        return false;
    }
    @autoreleasepool {
        AVURLAsset* asset = [AVURLAsset URLAssetWithURL:url options:nil];
        const double duration = CMTimeGetSeconds(asset.duration);
        if (!std::isfinite(duration) || duration <= 0.0) {
            return false;
        }
        const size_t sampleCount = static_cast<size_t>(std::ceil(duration));
        const size_t frameCount = static_cast<size_t>(std::ceil(duration * kAnalysisFrameRate));
        if (outCutTimesSec) {
            outCutTimesSec->clear();
        }

        // Step 1: Extract audio features using AudioAnalysisEngine
        NSError* audioError = nil;
        NSArray<AudioAnalysisResult*>* audioResults = [AudioAnalysisEngine analyzeAudio:asset
                                                                              durationSec:duration
                                                                            progressBlock:nil
                                                                                    error:&audioError];

        // Map AudioAnalysisResult to vectors for easy indexing (1Hz, by integer second)
        std::vector<double> audioRms(sampleCount, 0.0);
        std::vector<double> audioTransient(sampleCount, 0.0);
        std::vector<double> spectral(sampleCount, 0.5);
        std::vector<double> loudness(sampleCount, 0.0);
        std::vector<double> spectralFlux(sampleCount, 0.0);
        std::vector<double> spectralCentroid(sampleCount, 0.0);
        std::vector<double> zcr(sampleCount, 0.0);
        std::vector<double> bpm(sampleCount, 0.0);
        std::vector<double> tempoConf(sampleCount, 0.0);
        std::vector<double> rhythmEnergy(sampleCount, 0.0);
        std::vector<int> estimatedKey(sampleCount, -1);
	        std::vector<double> harmonicTension(sampleCount, 0.0);
	        std::vector<double> dialogueClarity(sampleCount, 0.0);
	        std::vector<double> speechProb(sampleCount, 0.0);
	        std::vector<double> stereoWidth(sampleCount, 0.0);
	        std::vector<double> reverbAmount(sampleCount, 0.0);
	        std::vector<double> dialogueEnergy(sampleCount, 0.0);
	        std::vector<double> musicEnergy(sampleCount, 0.0);
	        std::vector<double> effectsEnergy(sampleCount, 0.0);
	        std::vector<double> dialogueDom(sampleCount, 0.0);
	        std::vector<double> musicDom(sampleCount, 0.0);
	        std::vector<double> effectsDom(sampleCount, 0.0);
	        std::vector<std::array<double,12>> chroma(sampleCount);
	        for (auto& c : chroma) c.fill(0.0);

        if (audioResults && !audioError) {
            for (AudioAnalysisResult* result in audioResults) {
                size_t i = static_cast<size_t>(result.timeSeconds);
                if (i < sampleCount) {
                    audioRms[i] = result.rms;
                    audioTransient[i] = result.spectralFlux;
                    spectral[i] = result.spectralBalance;

                    loudness[i] = result.loudness;
                    spectralFlux[i] = result.spectralFlux;
                    spectralCentroid[i] = result.spectralCentroid;
                    zcr[i] = result.zeroCrossingRate;
                    bpm[i] = result.bpm;
                    tempoConf[i] = result.tempoConfidence;
                    rhythmEnergy[i] = result.rhythmEnergy;
                    estimatedKey[i] = result.estimatedKey;
	                    harmonicTension[i] = result.harmonicTension;
	                    dialogueClarity[i] = result.dialogueClarity;
	                    speechProb[i] = result.speechProbability;
	                    stereoWidth[i] = result.stereoWidth;
	                    reverbAmount[i] = result.reverbAmount;
                        dialogueEnergy[i] = result.dialogueEnergy;
                        musicEnergy[i] = result.musicEnergy;
                        effectsEnergy[i] = result.effectsEnergy;
	                    dialogueDom[i] = result.dialogueDominance;
	                    musicDom[i] = result.musicDominance;
	                    effectsDom[i] = result.effectsDominance;
	                    if (result.chroma && result.chroma.count == 12) {
	                        for (int k = 0; k < 12; ++k) chroma[i][k] = [result.chroma[k] doubleValue];
	                    }
	                }
	            }
	        }

        // Step 2: Extract subtitle features
        auto cues = LoadSubtitleCues(url);
        std::vector<double> subtitleWords;
        std::vector<double> expositionDensity;
        ComputeSubtitleGrids(cues, sampleCount, subtitleWords, expositionDensity);

        // Step 2.5: Process NLP features from subtitle text
        using namespace arcscope::nlp;
        SubtitleNLPProcessor nlpProcessor;
        auto nlpFeatures = nlpProcessor.process(cues, duration);

        // Verify NLP features are properly aligned (should match sampleCount)
        if (nlpFeatures.size() != sampleCount) {
            // If NLP processing failed or misaligned, clear and fill with zeros
            nlpFeatures.clear();
            nlpFeatures.resize(sampleCount);
            for (size_t i = 0; i < sampleCount; ++i) {
                nlpFeatures[i].timeSeconds = static_cast<double>(i) + 0.5;
            }
        }

        // Step 3: Extract video features at fixed 10fps, downsampled to fixed 720p height.
        // Output is aggregated to 1Hz (k+0.5) samples.
        std::vector<FrameStats> secAccum(sampleCount);
        std::vector<int> secCount(sampleCount, 0);
        std::vector<Cam16SecondStats> cam16Sec(sampleCount);
        std::vector<int> cam16Count(sampleCount, 0);

	        std::vector<double> frameMotion(frameCount, 0.0);   // optical flow mag (normalized) at each analysis frame
	        std::vector<double> frameCutScore(frameCount, 0.0); // luma diff at each analysis frame

	        // Motion per second is mean optical-flow magnitude over frames in that second
	        std::vector<double> motionPerSecond(sampleCount, 0.0);
	        std::vector<int> motionCountPerSecond(sampleCount, 0);

	        // Camera vs object motion + type (doc 5.6 overlays)
	        std::vector<double> cameraMotionPerSecond(sampleCount, 0.0);
	        std::vector<double> objectMotionPerSecond(sampleCount, 0.0);
	        std::vector<double> cameraUPerSecond(sampleCount, 0.0);
	        std::vector<double> cameraVPerSecond(sampleCount, 0.0);
	        std::vector<double> zoomPerSecond(sampleCount, 0.0);
	        std::vector<int> camMotionCountPerSecond(sampleCount, 0);
	        std::vector<int> cameraMotionTypePerSecond(sampleCount, 0);

        // Video decode (streaming)
        NSArray<AVAssetTrack*>* vtracks = [asset tracksWithMediaType:AVMediaTypeVideo];
        if (vtracks.count == 0) return false;
        AVAssetTrack* vtrack = vtracks.firstObject;

        NSError* readerErr = nil;
        AVAssetReader* reader = [AVAssetReader assetReaderWithAsset:asset error:&readerErr];
        if (!reader || readerErr) return false;

        NSDictionary* outSettings = @{
            (id)kCVPixelBufferPixelFormatTypeKey: @(kCVPixelFormatType_32BGRA),
            (id)kCVPixelBufferIOSurfacePropertiesKey: @{}
        };
        AVAssetReaderTrackOutput* output = [AVAssetReaderTrackOutput assetReaderTrackOutputWithTrack:vtrack outputSettings:outSettings];
        output.alwaysCopiesSampleData = NO;
        if (![reader canAddOutput:output]) return false;
        [reader addOutput:output];
        if (![reader startReading]) return false;

        CIContext* ciContext = [CIContext contextWithOptions:nil];
        CGColorSpaceRef colorSpace = CGColorSpaceCreateDeviceRGB();
        CGColorSpaceRef linearRec2020 = CGColorSpaceCreateWithName(kCGColorSpaceLinearITUR_2020);
        const CGAffineTransform preferredTransform = vtrack.preferredTransform;
        VNSequenceRequestHandler* flowHandler = [[VNSequenceRequestHandler alloc] init];
        Cam16ViewingConditions camVC;

        CMSampleBufferRef prevSB = nullptr;
        double prevT = -1e30;
        CMSampleBufferRef curSB = [output copyNextSampleBuffer];
        double curT = curSB ? CMTimeGetSeconds(CMSampleBufferGetPresentationTimeStamp(curSB)) : 1e30;

        CMSampleBufferRef lastChosenSB = nullptr;
        CVPixelBufferRef lastScaledPB = nullptr;

        CVPixelBufferRef prevScaledPB = nullptr;
        std::vector<uint8_t> lumaA, lumaB, lumaDiff;

        for (size_t fi = 0; fi < frameCount; ++fi) {
            const double targetTime = static_cast<double>(fi) / kAnalysisFrameRate;

            while (curSB && curT < targetTime) {
                if (prevSB) CFRelease(prevSB);
                prevSB = curSB;
                prevT = curT;
                curSB = [output copyNextSampleBuffer];
                curT = curSB ? CMTimeGetSeconds(CMSampleBufferGetPresentationTimeStamp(curSB)) : 1e30;
            }

            CMSampleBufferRef chosen = nullptr;
            if (prevSB && curSB) {
                const double dPrev = std::abs(targetTime - prevT);
                const double dCur  = std::abs(curT - targetTime);
                chosen = (dPrev <= dCur) ? prevSB : curSB;
            } else if (curSB) {
                chosen = curSB;
            } else if (prevSB) {
                chosen = prevSB;
            } else {
                break;
            }

            CVPixelBufferRef srcPB = CMSampleBufferGetImageBuffer(chosen);
            if (!srcPB) continue;

            CVPixelBufferRef scaledPB = nullptr;
            if (chosen == lastChosenSB && lastScaledPB) {
                scaledPB = lastScaledPB;
                CFRetain(scaledPB);
            } else {
                scaledPB = create_oriented_scaled_720p(srcPB, preferredTransform, ciContext, colorSpace);
                if (lastScaledPB) CFRelease(lastScaledPB);
                lastScaledPB = scaledPB ? (CVPixelBufferRef)CFRetain(scaledPB) : nullptr;
                lastChosenSB = chosen;
            }
            if (!scaledPB) continue;

            // Frame stats (color observations)
            const FrameStats fs = compute_frame_stats_from_argb(scaledPB);
            const int sec = std::max(0, std::min((int)sampleCount - 1, (int)std::floor(targetTime)));
            secAccum[(size_t)sec].brightness += fs.brightness;
            secAccum[(size_t)sec].saturation += fs.saturation;
            secAccum[(size_t)sec].warmth += fs.warmth;
            secAccum[(size_t)sec].grayness += fs.grayness;
            // avgHue: accumulate only if valid
            if (fs.avgHue >= 0.0) {
                secAccum[(size_t)sec].avgHue += fs.avgHue;
            }
            secCount[(size_t)sec] += 1;

            // Color science sampling (arcscope.md §3.1.2): pick 1–2 frames per second for CAM16-UCS statistics.
            // Deterministic choice for f_analysis=10fps: frames 2 and 7 of each second.
            if (linearRec2020) {
                const int fiInSec = (int)(fi % (size_t)kAnalysisFrameRate);
                if (fiInSec == 2 || fiInSec == 7) {
                    CVPixelBufferRef half = render_to_linear_rec2020_half(scaledPB, ciContext, linearRec2020);
                    if (half) {
                        Cam16SecondStats cs = cam16_stats_from_linear_rec2020_half(half, camVC);
                        cam16Sec[(size_t)sec].J += cs.J;
                        cam16Sec[(size_t)sec].M += cs.M;
                        cam16Sec[(size_t)sec].Jp += cs.Jp;
                        cam16Sec[(size_t)sec].ap += cs.ap;
                        cam16Sec[(size_t)sec].bp += cs.bp;
                        for (size_t b = 0; b < 12; ++b) cam16Sec[(size_t)sec].hueHist[b] += cs.hueHist[b];

                        // Harmony score: K-medians in CAM16-UCS on a uniform pixel grid.
                        std::vector<std::array<float,3>> pts;
                        CVPixelBufferLockBaseAddress(half, kCVPixelBufferLock_ReadOnly);
                        const size_t w2 = CVPixelBufferGetWidth(half);
                        const size_t h2 = CVPixelBufferGetHeight(half);
                        const size_t rb2 = CVPixelBufferGetBytesPerRow(half);
                        const uint16_t* b2 = (const uint16_t*)CVPixelBufferGetBaseAddress(half);
                        if (b2) {
                            const size_t pix = w2 * h2;
                            const size_t targetPts = 20000;
                            const size_t stride = std::max<size_t>(1, (size_t)std::floor(std::sqrt((double)pix / (double)targetPts)));
                            pts.reserve(std::min(pix, targetPts));
                            for (size_t y = 0; y < h2; y += stride) {
                                const uint16_t* row = (const uint16_t*)((const uint8_t*)b2 + y * rb2);
                                for (size_t x = 0; x < w2; x += stride) {
                                    const float r = half_to_float(row[4 * x + 0]);
                                    const float g = half_to_float(row[4 * x + 1]);
                                    const float b = half_to_float(row[4 * x + 2]);
                                    if (!(std::isfinite(r) && std::isfinite(g) && std::isfinite(b))) continue;
                                    const auto xyz = rec2020_linear_to_xyz(std::max(0.0f, r), std::max(0.0f, g), std::max(0.0f, b));
                                    const auto cam = xyz_to_cam16ucs(xyz[0], xyz[1], xyz[2], camVC);
                                    pts.push_back({(float)cam.Jp, (float)cam.ap, (float)cam.bp});
                                    if (pts.size() >= targetPts) break;
                                }
                                if (pts.size() >= targetPts) break;
                            }
                        }
                        CVPixelBufferUnlockBaseAddress(half, kCVPixelBufferLock_ReadOnly);
                        if (!pts.empty()) {
                            cam16Sec[(size_t)sec].harmony += cam16_harmony_kmedians(pts, 3);
                        }

                        cam16Count[(size_t)sec] += 1;
                        CFRelease(half);
                    }
                }
            }

	            // Motion (optical flow) + cut score (luma diff)
	            if (fi > 0 && prevScaledPB) {
	                const auto flow = optical_flow_stats(flowHandler, prevScaledPB, scaledPB);
	                const size_t fw = CVPixelBufferGetWidth(scaledPB);
	                const size_t fh = CVPixelBufferGetHeight(scaledPB);
	                const double diag = std::hypot((double)fw, (double)fh);
	                const double meanSq = std::max(0.0, flow.meanSq);
	                const double camSq = flow.meanU * flow.meanU + flow.meanV * flow.meanV;
	                const double resSq = std::max(0.0, meanSq - camSq);

	                const double totalRms = (diag > 1e-12) ? (std::sqrt(meanSq) / diag) : 0.0;
	                const double camU = (diag > 1e-12) ? (flow.meanU / diag) : 0.0;
	                const double camV = (diag > 1e-12) ? (flow.meanV / diag) : 0.0;
	                const double camMag = (diag > 1e-12) ? (std::sqrt(std::max(camSq, 0.0)) / diag) : 0.0;
	                const double objMag = (diag > 1e-12) ? (std::sqrt(resSq) / diag) : 0.0;
	                const double zoom = (diag > 1e-12) ? (flow.zoom / diag) : 0.0;

	                frameMotion[fi] = totalRms;
	                frameCutScore[fi] = mean_abs_luma_diff_argb(prevScaledPB, scaledPB, lumaA, lumaB, lumaDiff);

	                motionPerSecond[(size_t)sec] += frameMotion[fi];
	                motionCountPerSecond[(size_t)sec] += 1;

	                cameraMotionPerSecond[(size_t)sec] += camMag;
	                objectMotionPerSecond[(size_t)sec] += objMag;
	                cameraUPerSecond[(size_t)sec] += camU;
	                cameraVPerSecond[(size_t)sec] += camV;
	                zoomPerSecond[(size_t)sec] += zoom;
	                camMotionCountPerSecond[(size_t)sec] += 1;
	            }

            if (prevScaledPB) CFRelease(prevScaledPB);
            prevScaledPB = scaledPB; // transfer ownership
        }

        if (prevScaledPB) CFRelease(prevScaledPB);
        if (lastScaledPB) CFRelease(lastScaledPB);
        if (prevSB) CFRelease(prevSB);
        if (curSB) CFRelease(curSB);
        CGColorSpaceRelease(colorSpace);
        if (linearRec2020) CGColorSpaceRelease(linearRec2020);

        // Finalize CAM16 per-second means and ΔE in CAM16-UCS
        for (size_t i = 0; i < sampleCount; ++i) {
            const int denom = std::max(0, cam16Count[i]);
            if (denom > 0) {
                cam16Sec[i].J /= (double)denom;
                cam16Sec[i].M /= (double)denom;
                cam16Sec[i].Jp /= (double)denom;
                cam16Sec[i].ap /= (double)denom;
                cam16Sec[i].bp /= (double)denom;
                cam16Sec[i].h = wrap_deg_360(rad2deg(std::atan2(cam16Sec[i].bp, cam16Sec[i].ap)));

                double hs = 0.0;
                for (double v : cam16Sec[i].hueHist) hs += v;
                if (hs > 1e-12) for (double& v : cam16Sec[i].hueHist) v /= hs;
                else cam16Sec[i].hueHist.fill(0.0);

                cam16Sec[i].harmony = std::clamp(cam16Sec[i].harmony / (double)denom, 0.0, 1.0);
            } else {
                cam16Sec[i].h = -1.0;
                cam16Sec[i].hueHist.fill(0.0);
                cam16Sec[i].harmony = 0.0;
            }
        }
        for (size_t i = 0; i < sampleCount; ++i) {
            if (i == 0) {
                cam16Sec[i].deltaE = 0.0;
            } else {
                const double dJ = cam16Sec[i].Jp - cam16Sec[i - 1].Jp;
                const double da = cam16Sec[i].ap - cam16Sec[i - 1].ap;
                const double db = cam16Sec[i].bp - cam16Sec[i - 1].bp;
                cam16Sec[i].deltaE = std::sqrt(dJ*dJ + da*da + db*db);
            }
        }

	        // Finalize per-second motion mean
	        for (size_t i = 0; i < sampleCount; ++i) {
	            if (motionCountPerSecond[i] > 0) motionPerSecond[i] /= (double)motionCountPerSecond[i];
	            if (camMotionCountPerSecond[i] > 0) {
	                cameraMotionPerSecond[i] /= (double)camMotionCountPerSecond[i];
	                objectMotionPerSecond[i] /= (double)camMotionCountPerSecond[i];
	                cameraUPerSecond[i] /= (double)camMotionCountPerSecond[i];
	                cameraVPerSecond[i] /= (double)camMotionCountPerSecond[i];
	                zoomPerSecond[i] /= (double)camMotionCountPerSecond[i];
	            }
	        }

	        // Camera motion type (0=none,1=pan,2=tilt,3=push,4=pull,5=track), quantile-gated to avoid noise.
	        {
	            std::vector<double> strength;
	            strength.reserve(sampleCount);
	            for (size_t i = 0; i < sampleCount; ++i) strength.push_back(cameraMotionPerSecond[i] + objectMotionPerSecond[i]);
	            const double thr = quantile_copy(strength, 0.60);
	            for (size_t i = 0; i < sampleCount; ++i) {
	                const double s = cameraMotionPerSecond[i] + objectMotionPerSecond[i];
	                if (!(s > thr)) { cameraMotionTypePerSecond[i] = 0; continue; }

	                const double absZoom = std::abs(zoomPerSecond[i]);
	                const double trans = std::hypot(cameraUPerSecond[i], cameraVPerSecond[i]);
	                if (absZoom > trans) {
	                    cameraMotionTypePerSecond[i] = (zoomPerSecond[i] >= 0.0) ? 3 : 4; // push/pull
	                    continue;
	                }
	                const double absU = std::abs(cameraUPerSecond[i]);
	                const double absV = std::abs(cameraVPerSecond[i]);
	                const double ratio = (std::max(absU, absV) > 1e-12) ? (std::min(absU, absV) / std::max(absU, absV)) : 0.0;
	                if (ratio > 0.5) cameraMotionTypePerSecond[i] = 5; // track (diagonal)
	                else cameraMotionTypePerSecond[i] = (absU >= absV) ? 1 : 2; // pan/tilt
	            }
	        }

        // Detect cuts from frameCutScore using adaptive quantile threshold (no absolute thresholds).
        // Output BOTH:
        //   - exact cut timestamps (seconds) for ShotSegments (Core)
        //   - per-second cut counts for metadata convenience / backward compatibility
        std::vector<double> cutScoresForQ;
        cutScoresForQ.reserve(frameCutScore.size());
        for (size_t i = 1; i < frameCutScore.size(); ++i) cutScoresForQ.push_back(frameCutScore[i]);
        const double cutThr = quantile_copy(std::move(cutScoresForQ), 0.98);

        std::vector<double> cutTimes;
        cutTimes.reserve(sampleCount);
        {
            std::vector<size_t> cand;
            cand.reserve(frameCutScore.size() / 16);
            for (size_t fi = 1; fi < frameCutScore.size(); ++fi) {
                if (frameCutScore[fi] > cutThr) cand.push_back(fi);
            }
            // Cluster consecutive frames above threshold; keep the argmax within each cluster.
            size_t idx = 0;
            while (idx < cand.size()) {
                size_t runStart = idx;
                size_t runEnd = idx;
                while (runEnd + 1 < cand.size() && cand[runEnd + 1] == cand[runEnd] + 1) ++runEnd;

                size_t bestFi = cand[runStart];
                double bestScore = frameCutScore[bestFi];
                for (size_t j = runStart; j <= runEnd; ++j) {
                    const size_t fi = cand[j];
                    const double s = frameCutScore[fi];
                    if (s > bestScore) { bestScore = s; bestFi = fi; }
                }

                const double t = (double)bestFi / kAnalysisFrameRate;
                if (std::isfinite(t) && t > 0.0 && t < duration) {
                    if (cutTimes.empty() || std::abs(t - cutTimes.back()) > (1.0 / kAnalysisFrameRate)) {
                        cutTimes.push_back(t);
                    }
                }
                idx = runEnd + 1;
            }
        }

        if (outCutTimesSec) {
            *outCutTimesSec = cutTimes;
        }

        std::vector<double> cutsPerSecond(sampleCount, 0.0);
        for (double t : cutTimes) {
            const int sec = std::max(0, std::min((int)sampleCount - 1, (int)std::floor(t)));
            cutsPerSecond[(size_t)sec] += 1.0;
        }

        // Step 4: Run face tracking on entire video (motion-gated)
        NSMutableArray<NSNumber*>* motionArray = [NSMutableArray arrayWithCapacity:motionPerSecond.size()];
        for (double m : motionPerSecond) {
            [motionArray addObject:@(m)];
        }

        NSError* faceError = nil;
        NSArray<FaceTrack*>* tracks = [FaceTrackingEngine analyzeFaces:asset
                                                           durationSec:duration
                                                    analysisFrameRate:3.0
                                                         motionGating:motionArray
                                                        progressBlock:nil
                                                                error:&faceError];

        if (faceError || !tracks) {
            tracks = nil;  // Ensure nil if error
        }
        if (outTracks && tracks) {
            *outTracks = tracks;
        }

        // Step 5: Build per-second samples
        outSamples.clear();
        outSamples.reserve(sampleCount);
        for (size_t i = 0; i < sampleCount; ++i) {
            double audioValue = i < audioRms.size() ? audioRms[i] : 0.0;
            double expositionScore = i < expositionDensity.size() ? expositionDensity[i] : 0.0;

            // ARCHITECTURE CONTRACT (arcscope.md § 4.5.2):
            // Bridge 层只输出原始观测（tracks），不做解释
            // Face metrics 由 Core 的 FaceAnalyzer 从 tracks 生成
            // 这样确保单一数据源，避免 Bridge 和 Core 的重复推断

            ArcScopeBridgeSample sample{};
            // TIME INVARIANT CONTRACT ENFORCEMENT (§ I.R3):
            // sample.time_seconds MUST be center-aligned at i + 0.5
            // This ensures consistency with:
            //   - Core FilmEngine calculations (t_center = k + 0.5)
            //   - UI curve display (times[i] = i + 0.5)
            //   - Scene/Sequence boundary detection
            sample.time_seconds = static_cast<double>(i) + 0.5;
            const int cnt = secCount[i];
	            sample.cut_density = (i < cutsPerSecond.size()) ? cutsPerSecond[i] : 0.0;
	            sample.motion_amplitude = (i < motionPerSecond.size()) ? motionPerSecond[i] : 0.0;
	            sample.camera_motion = (i < cameraMotionPerSecond.size()) ? cameraMotionPerSecond[i] : 0.0;
	            sample.object_motion = (i < objectMotionPerSecond.size()) ? objectMotionPerSecond[i] : 0.0;
	            sample.camera_motion_type = (i < cameraMotionTypePerSecond.size()) ? cameraMotionTypePerSecond[i] : 0;
	            sample.audio_rms = audioValue;
	            sample.audio_transient = i < audioTransient.size() ? audioTransient[i] : 0.0;
	            sample.spectral_balance = i < spectral.size() ? spectral[i] : 0.5;

            // AudioAnalysisEngine pass-through (observation-level)
            sample.audio_loudness = i < loudness.size() ? loudness[i] : 0.0;
            sample.audio_spectral_flux = i < spectralFlux.size() ? spectralFlux[i] : 0.0;
            sample.audio_spectral_centroid = i < spectralCentroid.size() ? spectralCentroid[i] : 0.0;
            sample.audio_zero_crossing_rate = i < zcr.size() ? zcr[i] : 0.0;
            sample.audio_bpm = i < bpm.size() ? bpm[i] : 0.0;
            sample.audio_tempo_confidence = i < tempoConf.size() ? tempoConf[i] : 0.0;
            sample.audio_rhythm_energy = i < rhythmEnergy.size() ? rhythmEnergy[i] : 0.0;
            sample.audio_estimated_key = i < estimatedKey.size() ? estimatedKey[i] : -1;
	            sample.audio_harmonic_tension = i < harmonicTension.size() ? harmonicTension[i] : 0.0;
	            sample.audio_dialogue_clarity = i < dialogueClarity.size() ? dialogueClarity[i] : 0.0;
	            sample.audio_speech_probability = i < speechProb.size() ? speechProb[i] : 0.0;
	            sample.audio_stereo_width = i < stereoWidth.size() ? stereoWidth[i] : 0.0;
	            sample.audio_reverb_amount = i < reverbAmount.size() ? reverbAmount[i] : 0.0;
	            sample.audio_dialogue_energy = i < dialogueEnergy.size() ? dialogueEnergy[i] : 0.0;
	            sample.audio_music_energy = i < musicEnergy.size() ? musicEnergy[i] : 0.0;
	            sample.audio_effects_energy = i < effectsEnergy.size() ? effectsEnergy[i] : 0.0;
	            sample.audio_dialogue_dominance = i < dialogueDom.size() ? dialogueDom[i] : 0.0;
	            sample.audio_music_dominance = i < musicDom.size() ? musicDom[i] : 0.0;
	            sample.audio_effects_dominance = i < effectsDom.size() ? effectsDom[i] : 0.0;
	            for (int k = 0; k < 12; ++k) {
	                sample.audio_chroma[k] = (i < chroma.size()) ? chroma[i][k] : 0.0;
	            }

            if (cnt > 0) {
                sample.brightness = secAccum[i].brightness / (double)cnt;
                sample.saturation = secAccum[i].saturation / (double)cnt;
                sample.warmth = secAccum[i].warmth / (double)cnt;
                sample.grayness = secAccum[i].grayness / (double)cnt;
                // avgHue is only meaningful if a lot of frames are colorful; keep -1 if gray.
                sample.avg_hue = (secAccum[i].avgHue > 0.0) ? (secAccum[i].avgHue / (double)cnt) : -1.0;
            } else {
                sample.brightness = 0.0;
                sample.saturation = 0.0;
                sample.warmth = 0.0;
                sample.grayness = 0.0;
                sample.avg_hue = -1.0;
            }

            sample.cam16_j = cam16Sec[i].J;
            sample.cam16_m = cam16Sec[i].M;
            sample.cam16_h = cam16Sec[i].h;
            sample.cam16_jp = cam16Sec[i].Jp;
            sample.cam16_ap = cam16Sec[i].ap;
            sample.cam16_bp = cam16Sec[i].bp;
            sample.cam16_delta_e = cam16Sec[i].deltaE;
            sample.cam16_harmony = cam16Sec[i].harmony;
            for (int b = 0; b < 12; ++b) sample.cam16_hue_hist[b] = cam16Sec[i].hueHist[(size_t)b];

            // NLP features (from SubtitleNLPProcessor)
            // If NLP processing succeeded, these will be populated; otherwise they're 0
            sample.nlp_word_count = (i < nlpFeatures.size()) ? nlpFeatures[i].wordCount : 0;
            sample.nlp_sentence_complexity = (i < nlpFeatures.size()) ? nlpFeatures[i].sentenceComplexity : 0.0;
            sample.nlp_new_concept_ratio = (i < nlpFeatures.size()) ? nlpFeatures[i].newConceptRatio : 0.0;
            sample.nlp_entity_density = (i < nlpFeatures.size()) ? nlpFeatures[i].entityDensity : 0.0;

            // UNIFIED SUBTITLE DENSITY: Use NLP tokenizer word count for consistency
            sample.subtitle_density = clamp01(static_cast<double>(sample.nlp_word_count) / 12.0);

            // Exposition probability from overlap-weighted computation
            sample.exposition_probability = clamp01(expositionScore);

            // Bridge-level coarse info_density estimate (subtitle_density * exposition weight)
            sample.info_density = clamp01(sample.subtitle_density * 0.5 + sample.exposition_probability * 0.5);

            // Face metrics 置零：Core 的 FaceAnalyzer 会从 tracks 生成完整曲线
            // 这些字段保留是为了向后兼容，但不应被 Core 的 PCA/分析使用
            sample.face_presence = 0.0;
            sample.face_arousal = 0.0;
            sample.face_valence = 0.0;
            sample.face_confidence = 0.0;
            outSamples.push_back(sample);
        }
        return !outSamples.empty();
    }
}
