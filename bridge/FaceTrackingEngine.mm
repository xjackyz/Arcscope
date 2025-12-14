//
//  FaceTrackingEngine.mm
//  ArcScope Bridge - Face Tracking Implementation
//

#import "FaceTrackingEngine.h"
#import <CoreImage/CoreImage.h>
#include <optional>
#include <unordered_map>
#import <vector>
#import <algorithm>
#import <cmath>

@implementation FaceDetection
@end

@implementation FaceTrack

- (instancetype)init {
    self = [super init];
    if (self) {
        _trackId = -1;
        _timePoints = [NSMutableArray array];
        _boundingBoxes = [NSMutableArray array];
        _confidences = [NSMutableArray array];
        _faceAreas = [NSMutableArray array];
        _valences = [NSMutableArray array];
        _arousals = [NSMutableArray array];
        _expressions = [NSMutableArray array];
        _yaws = [NSMutableArray array];
        _rolls = [NSMutableArray array];
        _downcastRatios = [NSMutableArray array];
    }
    return self;
}

- (void)addDetection:(FaceDetection *)detection atTime:(double)timeSeconds {
    [self.timePoints addObject:@(timeSeconds)];
    [self.boundingBoxes addObject:[NSValue valueWithRect:detection.boundingBox]];
    [self.confidences addObject:@(detection.confidence)];
    [self.faceAreas addObject:@(detection.faceArea)];
    [self.valences addObject:@(detection.valence)];
    [self.arousals addObject:@(detection.arousal)];
    [self.expressions addObject:@(detection.expressionIntensity)];
    [self.yaws addObject:@(detection.yaw)];
    [self.rolls addObject:@(detection.roll)];
    [self.downcastRatios addObject:@(detection.downcastRatio)];
}

- (void)computeStatistics:(double)totalDuration {
    if (self.timePoints.count == 0) {
        self.appearanceDuration = 0.0;
        self.avgFaceArea = 0.0;
        self.avgConfidence = 0.0;
        self.score = 0.0;
        return;
    }

    // Appearance duration: count of unique integer seconds (not percentage)
    // Convert each timestamp to integer second via floor, then de-duplicate
    NSMutableSet<NSNumber *> *uniqueSecondsSet = [NSMutableSet set];
    for (NSNumber *timeNum in self.timePoints) {
        double t = [timeNum doubleValue];
        int sec = (int)floor(t);
        [uniqueSecondsSet addObject:@(sec)];
    }
    // Store as count of seconds (Core will divide by filmDurationSec to get Ti)
    self.appearanceDuration = (double)uniqueSecondsSet.count;

    // Average face area
    double sumArea = 0.0;
    for (NSNumber *area in self.faceAreas) {
        sumArea += [area doubleValue];
    }
    self.avgFaceArea = sumArea / self.faceAreas.count;

    // Average confidence
    double sumConf = 0.0;
    for (NSNumber *conf in self.confidences) {
        sumConf += [conf doubleValue];
    }
    self.avgConfidence = sumConf / self.confidences.count;

    // Score: (seconds / totalDuration) * S_i * C_i
    // Note: This score is for Bridge-internal use only; Core computes its own score
    double Ti = self.appearanceDuration / totalDuration;
    self.score = Ti * self.avgFaceArea * self.avgConfidence;
}

@end

@implementation FaceTrackingEngine

#pragma mark - CoreML Emotion Model (required, no heuristics)

namespace {

struct EmotionPrediction {
    double valence{0.0};     // [-1,1]
    double arousal{0.0};     // [0,1]
    double intensity{0.0};   // [0,1]
};

static inline double clamp(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static inline double clamp01(double v) { return clamp(v, 0.0, 1.0); }

std::optional<double> face_observation_angle(VNFaceObservation* obs, NSString* key) {
    if (!obs || !key) return std::nullopt;
    @try {
        id v = [obs valueForKey:key];
        if ([v isKindOfClass:[NSNumber class]]) {
            return ((NSNumber*)v).doubleValue;
        }
    } @catch (__unused NSException* ex) {
    }
    return std::nullopt;
}

// Landmark proxy in face coords: (nose-eye)/(mouth-eye).
// Smaller ratios imply the nose appears closer to eyes (head pitched down / downcast),
// larger ratios imply closer to mouth (upcast).
double downcast_ratio_from_landmarks(VNFaceObservation* obs) {
    if (!obs) return 0.5;
    VNFaceLandmarks2D* lm = obs.landmarks;
    if (!lm) return 0.5;

    auto mean_y = [](VNFaceLandmarkRegion2D* r) -> std::optional<double> {
        if (!r || r.pointCount == 0) return std::nullopt;
        const CGPoint* pts = r.normalizedPoints;
        if (!pts) return std::nullopt;
        double acc = 0.0;
        for (NSUInteger i = 0; i < r.pointCount; ++i) acc += pts[i].y;
        return acc / (double)r.pointCount;
    };

    auto le = mean_y(lm.leftEye);
    auto re = mean_y(lm.rightEye);
    auto no = mean_y(lm.nose);
    auto mo = mean_y(lm.outerLips);
    if (!no || !mo || (!le && !re)) return 0.5;

    const double eyeY = (le && re) ? (0.5 * (*le + *re)) : (le ? *le : *re);
    const double noseY = *no;
    const double mouthY = *mo;
    const double denom = std::max(1e-6, mouthY - eyeY);
    const double ratio = (noseY - eyeY) / denom;
    return clamp(ratio, 0.0, 1.0);
}

// Load a required CoreML model from app bundle resources:
//   app/Sources/ArcScopeApp/Resources/ArcScopeFaceEmotion.mlmodelc
MLModel* loadEmotionModelOrThrow() {
    static MLModel* cached = nil;
    static dispatch_once_t onceToken;
    static NSError* cachedError = nil;

    dispatch_once(&onceToken, ^{
        NSURL* url = [[NSBundle mainBundle] URLForResource:@"ArcScopeFaceEmotion" withExtension:@"mlmodelc"];
        if (!url) {
            cachedError = [NSError errorWithDomain:@"FaceTrackingEngine"
                                             code:-1001
                                         userInfo:@{NSLocalizedDescriptionKey:
                                                        @"Missing CoreML model: ArcScopeFaceEmotion.mlmodelc (place it in app Resources)."}];
            return;
        }
        MLModelConfiguration* cfg = [[MLModelConfiguration alloc] init];
        NSError* err = nil;
        cached = [MLModel modelWithContentsOfURL:url configuration:cfg error:&err];
        cachedError = err;
    });

    if (!cached) {
        @throw [NSException exceptionWithName:@"ArcScopeFaceModelError"
                                       reason:cachedError.localizedDescription ?: @"Failed to load ArcScopeFaceEmotion model"
                                     userInfo:@{@"error": cachedError ?: [NSNull null]}];
    }
    return cached;
}

// Discover the first image input feature in the model.
std::pair<NSString*, MLImageConstraint*> discoverImageInputOrThrow(MLModel* model) {
    NSDictionary<NSString*, MLFeatureDescription*>* inputs = model.modelDescription.inputDescriptionsByName;
    for (NSString* name in inputs) {
        MLFeatureDescription* desc = inputs[name];
        if (desc.type == MLFeatureTypeImage && desc.imageConstraint) {
            return {name, desc.imageConstraint};
        }
    }
    @throw [NSException exceptionWithName:@"ArcScopeFaceModelError"
                                   reason:@"ArcScopeFaceEmotion model must have an image input (MLFeatureTypeImage)."
                                 userInfo:nil];
}

// Parse model outputs into (valence, arousal, intensity).
EmotionPrediction parseEmotionOutputsOrThrow(id<MLFeatureProvider> out) {
    // Preferred: separate scalar outputs.
    auto getScalar = [&](NSString* key) -> std::optional<double> {
        if (!out.featureNames || ![out.featureNames containsObject:key]) return std::nullopt;
        MLFeatureValue* v = [out featureValueForName:key];
        if (!v) return std::nullopt;
        if (v.type == MLFeatureTypeDouble) return v.doubleValue;
        if (v.type == MLFeatureTypeInt64) return (double)v.int64Value;
        return std::nullopt;
    };

    auto val = getScalar(@"valence");
    auto aro = getScalar(@"arousal");
    auto inten = getScalar(@"intensity");
    if (val && aro && inten) {
        EmotionPrediction p;
        p.valence = clamp(*val, -1.0, 1.0);
        p.arousal = clamp(*aro, 0.0, 1.0);
        p.intensity = clamp(*inten, 0.0, 1.0);
        return p;
    }

    // Alternative: a single multiarray output (length >= 3).
    for (NSString* name in out.featureNames) {
        MLFeatureValue* v = [out featureValueForName:name];
        if (!v || v.type != MLFeatureTypeMultiArray) continue;
        MLMultiArray* a = v.multiArrayValue;
        if (!a || a.count < 3) continue;

        // Flatten in row-major order; assume first 3 values are (valence, arousal, intensity).
        double v0 = ((NSNumber*)a[0]).doubleValue;
        double v1 = ((NSNumber*)a[1]).doubleValue;
        double v2 = ((NSNumber*)a[2]).doubleValue;

        EmotionPrediction p;
        p.valence = clamp(v0, -1.0, 1.0);
        p.arousal = clamp(v1, 0.0, 1.0);
        p.intensity = clamp(v2, 0.0, 1.0);
        return p;
    }

    @throw [NSException exceptionWithName:@"ArcScopeFaceModelError"
                                   reason:@"ArcScopeFaceEmotion outputs must provide valence/arousal/intensity (scalars or a 3+ MLMultiArray)."
                                 userInfo:nil];
}

// Crop face region from frame and resize to model input constraint.
// Uses scale-fill semantics to match Vision's typical CoreML preprocessing.
CVPixelBufferRef createFaceInputPixelBuffer(CVPixelBufferRef frame,
                                            VNFaceObservation* face,
                                            MLImageConstraint* constraint) {
    if (!frame || !face || !constraint) return nil;

    const size_t frameW = CVPixelBufferGetWidth(frame);
    const size_t frameH = CVPixelBufferGetHeight(frame);
    if (frameW == 0 || frameH == 0) return nil;

    const CGRect faceRect = VNImageRectForNormalizedRect(face.boundingBox, (int)frameW, (int)frameH);
    if (faceRect.size.width < 2 || faceRect.size.height < 2) return nil;

    CIImage* ci = [CIImage imageWithCVPixelBuffer:frame];
    CIImage* cropped = [ci imageByCroppingToRect:faceRect];

    const int targetW = (int)constraint.pixelsWide;
    const int targetH = (int)constraint.pixelsHigh;
    if (targetW <= 0 || targetH <= 0) return nil;

    const double sx = (double)targetW / std::max(1.0, faceRect.size.width);
    const double sy = (double)targetH / std::max(1.0, faceRect.size.height);
    const double scale = std::max(sx, sy); // scale fill

    CIImage* scaled = [cropped imageByApplyingTransform:CGAffineTransformMakeScale((CGFloat)scale, (CGFloat)scale)];
    const CGRect scaledExtent = scaled.extent;
    const double cx = CGRectGetMidX(scaledExtent);
    const double cy = CGRectGetMidY(scaledExtent);
    const CGRect outRect = CGRectMake(cx - targetW * 0.5, cy - targetH * 0.5, targetW, targetH);
    CIImage* final = [scaled imageByCroppingToRect:outRect];

    const OSType fmt = constraint.pixelFormatType != 0 ? constraint.pixelFormatType : kCVPixelFormatType_32BGRA;
    NSDictionary* attrs = @{
        (id)kCVPixelBufferCGImageCompatibilityKey: @YES,
        (id)kCVPixelBufferCGBitmapContextCompatibilityKey: @YES
    };

    CVPixelBufferRef outPB = nil;
    CVReturn rc = CVPixelBufferCreate(kCFAllocatorDefault, targetW, targetH, fmt,
                                      (__bridge CFDictionaryRef)attrs, &outPB);
    if (rc != kCVReturnSuccess || !outPB) return nil;

    static CIContext* ctx = nil;
    static dispatch_once_t once;
    dispatch_once(&once, ^{
        ctx = [CIContext contextWithOptions:@{kCIContextUseSoftwareRenderer: @NO}];
    });

    CGColorSpaceRef cs = CGColorSpaceCreateDeviceRGB();
    [ctx render:final toCVPixelBuffer:outPB bounds:CGRectMake(0, 0, targetW, targetH) colorSpace:cs];
    CGColorSpaceRelease(cs);
    return outPB;
}

EmotionPrediction inferEmotionOrThrow(MLModel* model,
                                     NSString* inputName,
                                     MLImageConstraint* constraint,
                                     CVPixelBufferRef frame,
                                     VNFaceObservation* face) {
    CVPixelBufferRef inputPB = createFaceInputPixelBuffer(frame, face, constraint);
    if (!inputPB) {
        @throw [NSException exceptionWithName:@"ArcScopeFaceModelError"
                                       reason:@"Failed to build face crop input for CoreML."
                                     userInfo:nil];
    }

    NSError* err = nil;
    MLFeatureValue* imageValue = [MLFeatureValue featureValueWithPixelBuffer:inputPB];
    NSDictionary<NSString*, MLFeatureValue*>* dict = @{ inputName: imageValue };
    id<MLFeatureProvider> provider = [[MLDictionaryFeatureProvider alloc] initWithDictionary:dict error:&err];
    if (!provider || err) {
        CVPixelBufferRelease(inputPB);
        @throw [NSException exceptionWithName:@"ArcScopeFaceModelError"
                                       reason:(err.localizedDescription ?: @"Failed to create ML feature provider.")
                                     userInfo:nil];
    }

    id<MLFeatureProvider> out = [model predictionFromFeatures:provider error:&err];
    CVPixelBufferRelease(inputPB);
    if (!out || err) {
        @throw [NSException exceptionWithName:@"ArcScopeFaceModelError"
                                       reason:(err.localizedDescription ?: @"CoreML prediction failed.")
                                     userInfo:nil];
    }

    return parseEmotionOutputsOrThrow(out);
}

} // namespace

#pragma mark - Main Analysis Pipeline

+ (nullable NSArray<FaceTrack *> *)analyzeFaces:(AVAsset *)asset
                                     durationSec:(double)durationSec
                              analysisFrameRate:(double)analysisFrameRate
                                   motionGating:(nullable NSArray<NSNumber *> *)motionPerSecond
                                  progressBlock:(nullable void (^)(double progress))progressBlock
                                          error:(NSError **)error {

    // Get first video track
    NSArray<AVAssetTrack *> *videoTracks = [asset tracksWithMediaType:AVMediaTypeVideo];
    if (videoTracks.count == 0) {
        if (error) {
            *error = [NSError errorWithDomain:@"FaceTrackingEngine"
                                         code:-1
                                     userInfo:@{NSLocalizedDescriptionKey: @"No video track found"}];
        }
        return nil;
    }

    // Setup reader
    AVAssetReader *reader = [[AVAssetReader alloc] initWithAsset:asset error:error];
    if (!reader) return nil;

    NSDictionary *outputSettings = @{
        (id)kCVPixelBufferPixelFormatTypeKey: @(kCVPixelFormatType_32BGRA)
    };

    AVAssetReaderTrackOutput *output = [[AVAssetReaderTrackOutput alloc] initWithTrack:videoTracks[0]
                                                                          outputSettings:outputSettings];
    [reader addOutput:output];

    if (![reader startReading]) {
        if (error) *error = reader.error;
        return nil;
    }

    // Setup face detection request (need landmarks for downcast proxy).
    VNDetectFaceLandmarksRequest *faceRequest = [[VNDetectFaceLandmarksRequest alloc] init];

    // Required CoreML emotion model (no heuristics)
    MLModel* emotionModel = nil;
    NSString* emotionInputName = nil;
    MLImageConstraint* emotionConstraint = nil;
    @try {
        emotionModel = loadEmotionModelOrThrow();
        auto in = discoverImageInputOrThrow(emotionModel);
        emotionInputName = in.first;
        emotionConstraint = in.second;
    } @catch (NSException* ex) {
        if (error) {
            *error = [NSError errorWithDomain:@"FaceTrackingEngine"
                                         code:-1002
                                     userInfo:@{NSLocalizedDescriptionKey: ex.reason ?: @"Face emotion model error"}];
        }
        return nil;
    }

    // Tracks
    NSMutableArray<FaceTrack *> *tracks = [NSMutableArray array];
    int nextTrackId = 0;

    // Frame processing
    double currentTime = 0.0;
    double frameDuration = 1.0 / analysisFrameRate;
    int frameCount = 0;
    int totalFrames = (int)(durationSec * analysisFrameRate);

    // Motion gating threshold (Q_M(0.6) from arcscope.md 3.1.2)
    double motionThreshold = 0.6;
    if (motionPerSecond && motionPerSecond.count > 0) {
        NSArray *sorted = [motionPerSecond sortedArrayUsingSelector:@selector(compare:)];
        int idx = (int)llround(0.6 * (double)std::max<NSInteger>(1, (NSInteger)sorted.count - 1));
        idx = std::max(0, std::min((int)sorted.count - 1, idx));
        motionThreshold = [sorted[(NSUInteger)idx] doubleValue];
    }

    // Cache per-track emotion when motion-gating allows reuse (doc 3.1.2).
    std::unordered_map<int, EmotionPrediction> cachedByTrack;
    int lastSecondProcessed = -1;
    bool prevSecondHadFace = false;       // presence over the previous second
    bool thisSecondHadFace = false;       // presence accumulated for the current second

    while (reader.status == AVAssetReaderStatusReading) {
        CMSampleBufferRef sampleBuffer = [output copyNextSampleBuffer];
        if (!sampleBuffer) break;

        CVImageBufferRef imageBuffer = CMSampleBufferGetImageBuffer(sampleBuffer);
        if (!imageBuffer) {
            CFRelease(sampleBuffer);
            continue;
        }

        // Skip frames to achieve desired frame rate
        CMTime presentationTime = CMSampleBufferGetPresentationTimeStamp(sampleBuffer);
        double timestamp = CMTimeGetSeconds(presentationTime);

        if (timestamp >= currentTime) {
            currentTime += frameDuration;
            int currentSecond = (int)floor(timestamp);

            // Detect faces
            VNImageRequestHandler *handler = [[VNImageRequestHandler alloc]
                                               initWithCVPixelBuffer:imageBuffer
                                               options:@{}];

            [handler performRequests:@[faceRequest] error:nil];

            NSArray<VNFaceObservation *> *faceObservations = faceRequest.results;
            NSMutableArray<FaceDetection *> *detections = [NSMutableArray array];

	            for (VNFaceObservation *obs in faceObservations) {
	                FaceDetection *det = [[FaceDetection alloc] init];
	                det.boundingBox = obs.boundingBox;
	                det.confidence = obs.confidence;

	                // Calculate face area (bbox area / frame area)
	                det.faceArea = obs.boundingBox.size.width * obs.boundingBox.size.height;

                    // Head pose (radians, observation-level)
                    det.yaw = face_observation_angle(obs, @"yaw").value_or(0.0);
                    det.roll = face_observation_angle(obs, @"roll").value_or(0.0);

                    // Downcast proxy from landmarks (observation-level)
                    det.downcastRatio = downcast_ratio_from_landmarks(obs);

	                [detections addObject:det];
	            }

            const bool frameHadFace = (detections.count > 0);
            if (currentSecond != lastSecondProcessed) {
                // Advance second boundary: snapshot previous second's presence.
                if (lastSecondProcessed >= 0) prevSecondHadFace = thisSecondHadFace;
                lastSecondProcessed = currentSecond;
                thisSecondHadFace = false;
            }
            thisSecondHadFace = thisSecondHadFace || frameHadFace;

            // Associate detections with tracks
            NSArray *associations = [self associateDetectionsWithTracks:detections tracks:tracks threshold:0.3];

            // Motion gating: only re-analyze expression if motion above threshold OR face presence toggles 0<->1 (doc 3.1.2).
            bool shouldAnalyzeExpression = true;
            if (motionPerSecond && currentSecond >= 0 && currentSecond < (int)motionPerSecond.count) {
                double motion = [motionPerSecond[(NSUInteger)currentSecond] doubleValue];
                shouldAnalyzeExpression = (motion >= motionThreshold);
            }
            const bool toggledPresence = (thisSecondHadFace != prevSecondHadFace);
            if (toggledPresence) shouldAnalyzeExpression = true;

            // Build a mapping detIdx -> trackIdx (track array index).
            std::vector<int> detToTrack(detections.count, -1);
            for (NSArray* pair in associations) {
                int detIdx = [pair[0] intValue];
                int trackIdx = [pair[1] intValue];
                if (detIdx >= 0 && detIdx < (int)detections.count) detToTrack[(size_t)detIdx] = trackIdx;
            }

            // Emotion inference / reuse per detection (never heuristic).
            for (NSUInteger detIdx = 0; detIdx < detections.count; ++detIdx) {
                FaceDetection* det = detections[detIdx];
                VNFaceObservation* obs = faceObservations[detIdx];
                const int trackIdx = (detIdx < detToTrack.size()) ? detToTrack[(size_t)detIdx] : -1;

                bool inferred = false;
                if (!shouldAnalyzeExpression && trackIdx >= 0) {
                    auto it = cachedByTrack.find(trackIdx);
                    if (it != cachedByTrack.end()) {
                        det.valence = it->second.valence;
                        det.arousal = it->second.arousal;
                        det.expressionIntensity = it->second.intensity;
                        inferred = true;
                    }
                }

                if (!inferred && frameHadFace) {
                    @try {
                        EmotionPrediction p = inferEmotionOrThrow(emotionModel, emotionInputName, emotionConstraint, (CVPixelBufferRef)imageBuffer, obs);
                        det.valence = p.valence;
                        det.arousal = p.arousal;
                        det.expressionIntensity = p.intensity;
                        if (trackIdx >= 0) cachedByTrack[trackIdx] = p;
                    } @catch (NSException* ex) {
                        if (error) {
                            *error = [NSError errorWithDomain:@"FaceTrackingEngine"
                                                         code:-1003
                                                     userInfo:@{NSLocalizedDescriptionKey: ex.reason ?: @"Face emotion inference failed"}];
                        }
                        CFRelease(sampleBuffer);
                        return nil;
                    }
                }
            }

            // Update tracks
            NSMutableSet *updatedTrackIds = [NSMutableSet set];
            for (NSArray *pair in associations) {
                int detIdx = [pair[0] intValue];
                int trackId = [pair[1] intValue];

                FaceDetection *det = detections[detIdx];
                FaceTrack *track = tracks[trackId];
                [track addDetection:det atTime:timestamp];
                [updatedTrackIds addObject:@(trackId)];
            }

            // Create new tracks for unmatched detections
            for (int i = 0; i < detections.count; i++) {
                BOOL matched = NO;
                for (NSArray *pair in associations) {
                    if ([pair[0] intValue] == i) {
                        matched = YES;
                        break;
                    }
                }

                if (!matched) {
                    FaceTrack *newTrack = [[FaceTrack alloc] init];
                    newTrack.trackId = nextTrackId++;
                    [newTrack addDetection:detections[i] atTime:timestamp];
                    [tracks addObject:newTrack];

                    // Cache emotion for this new track (so motion gating reuse is well-defined).
                    FaceDetection* det = detections[i];
                    cachedByTrack[(int)tracks.count - 1] = EmotionPrediction{det.valence, det.arousal, det.expressionIntensity};
                }
            }

            frameCount++;
            if (progressBlock && frameCount % 30 == 0) {
                progressBlock((double)frameCount / totalFrames);
            }
        }

        CFRelease(sampleBuffer);
    }

    // Compute track statistics
    for (FaceTrack *track in tracks) {
        [track computeStatistics:durationSec];
    }

    // ARCHITECTURE CONTRACT:
    // Bridge only outputs raw tracks. Core's FaceAnalyzer will:
    // - Identify dominant track using adaptive quantiles
    // - Calculate close-up weights
    // - Build per-second FaceFeatures with presence/valence/arousal
    return tracks;
}

#pragma mark - Detection-Track Association

+ (NSArray *)associateDetectionsWithTracks:(NSArray<FaceDetection *> *)detections
                                    tracks:(NSArray<FaceTrack *> *)tracks
                                 threshold:(double)threshold {

    NSMutableArray *associations = [NSMutableArray array];

    if (detections.count == 0 || tracks.count == 0) {
        return associations;
    }

    // Build cost matrix: IoU between each detection and last bbox of each track
    std::vector<std::vector<double>> iouMatrix(detections.count, std::vector<double>(tracks.count, 0.0));

    for (NSUInteger i = 0; i < detections.count; i++) {
        FaceDetection *det = detections[i];

        for (NSUInteger j = 0; j < tracks.count; j++) {
            FaceTrack *track = tracks[j];

            if (track.boundingBoxes.count > 0) {
                CGRect lastBox = [track.boundingBoxes.lastObject rectValue];
                double iou = [self calculateIoU:det.boundingBox box2:lastBox];
                iouMatrix[i][j] = iou;
            }
        }
    }

    // Greedy matching: assign each detection to best track above threshold
    NSMutableSet *matchedDetections = [NSMutableSet set];
    NSMutableSet *matchedTracks = [NSMutableSet set];

    for (NSUInteger i = 0; i < detections.count; i++) {
        double maxIoU = threshold;
        int bestTrack = -1;

        for (NSUInteger j = 0; j < tracks.count; j++) {
            if (![matchedTracks containsObject:@(j)] && iouMatrix[i][j] > maxIoU) {
                maxIoU = iouMatrix[i][j];
                bestTrack = (int)j;
            }
        }

        if (bestTrack >= 0) {
            [associations addObject:@[@(i), @(bestTrack)]];
            [matchedDetections addObject:@(i)];
            [matchedTracks addObject:@(bestTrack)];
        }
    }

    return associations;
}

+ (double)calculateIoU:(CGRect)box1 box2:(CGRect)box2 {
    CGRect intersection = CGRectIntersection(box1, box2);

    if (CGRectIsNull(intersection)) {
        return 0.0;
    }

    double intersectionArea = intersection.size.width * intersection.size.height;
    double area1 = box1.size.width * box1.size.height;
    double area2 = box2.size.width * box2.size.height;
    double unionArea = area1 + area2 - intersectionArea;

    if (unionArea < 1e-6) {
        return 0.0;
    }

    return intersectionArea / unionArea;
}

@end
