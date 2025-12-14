//
//  FaceTrackingEngine.h
//  ArcScope Bridge - Multi-Object Tracking for Faces
//
//  Implements arcscope.md section 4.5.2
//  - Vision-based face detection
//  - Multi-object tracking (MOT) across frames
//  - Dominant face trajectory identification
//  - Expression/emotion analysis
//

#import <Foundation/Foundation.h>
#import <Vision/Vision.h>
#import <CoreML/CoreML.h>
#import <AVFoundation/AVFoundation.h>

NS_ASSUME_NONNULL_BEGIN

/**
 * Single face detection in a frame
 */
@interface FaceDetection : NSObject
@property (nonatomic, assign) CGRect boundingBox;      // Normalized [0,1]
@property (nonatomic, assign) double confidence;       // Detection confidence
@property (nonatomic, assign) double faceArea;         // Bbox area / frame area
@property (nonatomic, assign) double valence;          // Emotion valence [-1,1]
@property (nonatomic, assign) double arousal;          // Emotion arousal [0,1]
@property (nonatomic, assign) double expressionIntensity;  // Expression strength
@property (nonatomic, assign) double yaw;              // Head yaw (radians)
@property (nonatomic, assign) double roll;             // Head roll (radians)
@property (nonatomic, assign) double downcastRatio;    // (nose-eye)/(mouth-eye) in face coords
@end

/**
 * Face track across time
 */
@interface FaceTrack : NSObject
@property (nonatomic, assign) int trackId;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *timePoints;     // Seconds
@property (nonatomic, strong) NSMutableArray<NSValue *> *boundingBoxes;   // CGRect
@property (nonatomic, strong) NSMutableArray<NSNumber *> *confidences;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *faceAreas;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *valences;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *arousals;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *expressions;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *yaws;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *rolls;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *downcastRatios;

// Statistics (computed after tracking)
@property (nonatomic, assign) double appearanceDuration;   // Total seconds on screen
@property (nonatomic, assign) double avgFaceArea;
@property (nonatomic, assign) double avgConfidence;
@property (nonatomic, assign) double score;                // T_i * S_i * C_i

- (void)addDetection:(FaceDetection *)detection atTime:(double)timeSeconds;
- (void)computeStatistics:(double)totalDuration;
@end

/**
 * Face Tracking Engine
 *
 * Performs multi-object tracking using Vision framework
 */
@interface FaceTrackingEngine : NSObject

/**
 * Analyze faces in video and track all faces
 *
 * @param asset Video asset
 * @param durationSec Total duration
 * @param analysisFrameRate Frame rate for analysis (default 3 fps per arcscope.md 3.1.2)
 * @param motionGating Motion threshold for expression re-analysis
 * @param progressBlock Progress callback
 * @return Array of FaceTrack objects
 *
 * NOTE: This method ONLY extracts raw face tracks.
 * All interpretation (dominant track identification, close-up weight, etc.)
 * is performed by Core's FaceAnalyzer.
 */
+ (nullable NSArray<FaceTrack *> *)analyzeFaces:(AVAsset *)asset
                                     durationSec:(double)durationSec
                              analysisFrameRate:(double)analysisFrameRate
                                   motionGating:(nullable NSArray<NSNumber *> *)motionPerSecond
                                  progressBlock:(nullable void (^)(double progress))progressBlock
                                          error:(NSError **)error;

/**
 * Associate detections with existing tracks using IoU matching
 *
 * @param detections Current frame detections
 * @param tracks Existing tracks
 * @param threshold IoU threshold for matching (default 0.3)
 * @return Array of (detection index, track ID) pairs
 */
+ (NSArray *)associateDetectionsWithTracks:(NSArray<FaceDetection *> *)detections
                                    tracks:(NSArray<FaceTrack *> *)tracks
                                 threshold:(double)threshold;

/**
 * Calculate IoU (Intersection over Union) between two bounding boxes
 */
+ (double)calculateIoU:(CGRect)box1 box2:(CGRect)box2;

@end

NS_ASSUME_NONNULL_END
