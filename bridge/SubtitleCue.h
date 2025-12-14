#pragma once

#import <Foundation/Foundation.h>
#include <string>

/**
 * Shared SubtitleCue structure for Bridge layer
 *
 * IMPORTANT:
 * - Uses std::string instead of NSString* to avoid ARC issues in STL containers
 * - Shared across VideoFeatureExtractor and SubtitleNLPProcessor
 * - Do NOT place in anonymous namespace
 */
namespace arcscope {
namespace bridge {

struct SubtitleCue {
    double start{0.0};        // Start time in seconds
    double end{0.0};          // End time in seconds
    int wordCount{0};         // Word count (tokenizer-based, language-aware)
    bool isExposition{false}; // High word density flag
    std::string text_utf8;    // UTF-8 encoded text (safe for STL containers)
};

}  // namespace bridge
}  // namespace arcscope
