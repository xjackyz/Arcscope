#pragma once

#include <vector>

#include "arcscope_bridge.h"

#ifdef __OBJC__
@class NSURL;
@class NSArray;
@class FaceTrack;
#else
typedef void NSURL;
typedef void NSArray;
typedef void FaceTrack;
#endif

bool ExtractSamplesFromVideo(NSURL* url,
                             std::vector<ArcScopeBridgeSample>& outSamples,
                             NSArray** outTracks,
                             std::vector<double>* outCutTimesSec);
