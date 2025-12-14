#import "SubtitleNLP.h"

#import <Foundation/Foundation.h>
#import <NaturalLanguage/NaturalLanguage.h>

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

namespace arcscope {
namespace nlp {

namespace {

/**
 * Extract tokens from text using NLTokenizer
 */
std::vector<NSString*> tokenizeText(NSString* text) {
    __block std::vector<NSString*> tokens;
    if (!text || text.length == 0) {
        return tokens;
    }

    NLTokenizer* tokenizer = [[NLTokenizer alloc] initWithUnit:NLTokenUnitWord];
    [tokenizer setString:text];

    [tokenizer enumerateTokensInRange:NSMakeRange(0, text.length)
                           usingBlock:^(NSRange tokenRange, NLTokenizerAttributes attributes, BOOL* stop) {
        NSString* token = [text substringWithRange:tokenRange];
        if (token.length > 0) {
            tokens.push_back([token lowercaseString]);
        }
    }];

    return tokens;
}

/**
 * Calculate sentence complexity from text
 * Uses: average sentence length, punctuation density, connector word ratio
 * Supports both English and Chinese text
 */
double calculateSentenceComplexity(NSString* text, const std::vector<NSString*>& tokens) {
    if (!text || text.length == 0 || tokens.empty()) {
        return 0.0;
    }

    // Detect dominant language
    NLLanguageRecognizer* recognizer = [[NLLanguageRecognizer alloc] init];
    [recognizer processString:text];
    NLLanguage dominantLang = recognizer.dominantLanguage;
    BOOL isChinese = dominantLang && [dominantLang hasPrefix:@"zh"];

    // Count sentences (rough heuristic)
    NSUInteger sentenceCount = 1;
    for (NSUInteger i = 0; i < text.length; ++i) {
        unichar ch = [text characterAtIndex:i];
        if (ch == '.' || ch == '!' || ch == '?' || ch == 0x3002 || ch == 0xFF01 || ch == 0xFF1F) {  // Chinese punctuation
            sentenceCount++;
        }
    }
    sentenceCount = std::max<NSUInteger>(1, sentenceCount);

    // Average sentence length (in words)
    double avgSentenceLength = static_cast<double>(tokens.size()) / static_cast<double>(sentenceCount);

    // Connector word density (proxy for clause complexity)
    static NSSet<NSString*>* englishConnectors = [NSSet setWithArray:@[
        @"because", @"although", @"however", @"therefore", @"moreover",
        @"furthermore", @"nevertheless", @"consequently", @"meanwhile",
        @"whereas", @"while", @"since", @"unless", @"until", @"after",
        @"before", @"when", @"where", @"which", @"who", @"that"
    ]];

    static NSSet<NSString*>* chineseConnectors = [NSSet setWithArray:@[
        @"因为", @"所以", @"但是", @"然而", @"而且", @"并且", @"或者",
        @"虽然", @"如果", @"假如", @"由于", @"于是", @"然后", @"接着",
        @"不过", @"可是", @"而", @"却", @"况且", @"何况", @"既然",
        @"即使", @"尽管", @"无论", @"不管", @"只要", @"只有", @"除非",
        @"的", @"了", @"着", @"把", @"被", @"将", @"对", @"向"
    ]];

    NSSet<NSString*>* connectorWords = isChinese ? chineseConnectors : englishConnectors;

    int connectorCount = 0;
    for (NSString* token : tokens) {
        if ([connectorWords containsObject:token]) {
            connectorCount++;
        }
    }
    double connectorRatio = static_cast<double>(connectorCount) / static_cast<double>(tokens.size());

    // Combine metrics (normalized to [0, 1])
    // Long sentences (>15 words) and high connector density → high complexity
    double lengthScore = std::min(1.0, avgSentenceLength / 15.0);
    double connectorScore = std::min(1.0, connectorRatio * 5.0);  // Scale up connector influence

    return std::min(1.0, 0.7 * lengthScore + 0.3 * connectorScore);
}

/**
 * Extract concept candidates (nouns, named entities)
 * Supports both English and Chinese text with language-specific filtering
 */
std::vector<NSString*> extractConcepts(NSString* text, NLTagger* tagger) {
    __block std::vector<NSString*> concepts;
    if (!text || text.length == 0) {
        return concepts;
    }

    // Detect dominant language
    NLLanguageRecognizer* recognizer = [[NLLanguageRecognizer alloc] init];
    [recognizer processString:text];
    NLLanguage dominantLang = recognizer.dominantLanguage;
    BOOL isChinese = dominantLang && [dominantLang hasPrefix:@"zh"];

    [tagger setString:text];

    // Extract nouns and named entities with proper options
    NLTaggerOptions options = NLTaggerOmitWhitespace | NLTaggerOmitPunctuation | NLTaggerJoinNames;
    [tagger enumerateTagsInRange:NSMakeRange(0, text.length)
                            unit:NLTokenUnitWord
                          scheme:NLTagSchemeNameTypeOrLexicalClass
                         options:options
                      usingBlock:^(NLTag tag, NSRange tokenRange, BOOL* stop) {
        if (!tag) return;

        // Accept nouns and named entities
        BOOL isNoun = [tag isEqualToString:NLTagNoun];
        BOOL isEntity = [tag isEqualToString:NLTagPersonalName] ||
                        [tag isEqualToString:NLTagPlaceName] ||
                        [tag isEqualToString:NLTagOrganizationName];

        if (isNoun || isEntity) {
            NSString* token = [[text substringWithRange:tokenRange] lowercaseString];

            // Language-specific filtering
            NSUInteger minLength = isChinese ? 1 : 3;  // Chinese: 1+ chars, English: 3+ chars

            if (token.length >= minLength) {
                concepts.push_back(token);
            }
        }
    }];

    return concepts;
}

/**
 * Calculate new concept ratio using sliding window
 * Concepts seen within last 30 seconds are considered "known"
 */
double calculateNewConceptRatio(
    const std::vector<NSString*>& concepts,
    std::unordered_map<std::string, double>& conceptLastSeen,
    double currentTime,
    double windowSeconds = 30.0
) {
    if (concepts.empty()) {
        return 0.0;
    }

    int newCount = 0;
    for (NSString* concept : concepts) {
        std::string key = [concept UTF8String];
        auto it = conceptLastSeen.find(key);

        bool isNew = (it == conceptLastSeen.end()) ||
                     ((currentTime - it->second) > windowSeconds);

        if (isNew) {
            newCount++;
        }

        // Update last seen time
        conceptLastSeen[key] = currentTime;
    }

    return static_cast<double>(newCount) / static_cast<double>(concepts.size());
}

/**
 * Calculate entity density (named entities per word)
 */
double calculateEntityDensity(NSString* text, NLTagger* tagger, int wordCount) {
    if (!text || text.length == 0 || wordCount == 0) {
        return 0.0;
    }

    [tagger setString:text];

    __block int entityCount = 0;
    NLTaggerOptions options = NLTaggerOmitWhitespace | NLTaggerOmitPunctuation | NLTaggerJoinNames;
    [tagger enumerateTagsInRange:NSMakeRange(0, text.length)
                            unit:NLTokenUnitWord
                          scheme:NLTagSchemeNameType
                         options:options
                      usingBlock:^(NLTag tag, NSRange tokenRange, BOOL* stop) {
        if ([tag isEqualToString:NLTagPersonalName] ||
            [tag isEqualToString:NLTagPlaceName] ||
            [tag isEqualToString:NLTagOrganizationName]) {
            entityCount++;
        }
    }];

    return std::min(1.0, static_cast<double>(entityCount) / static_cast<double>(wordCount));
}

}  // anonymous namespace

// ============================================================================
// SubtitleNLPProcessor Implementation
// ============================================================================

class SubtitleNLPProcessor::Impl {
public:
    NLTagger* tagger_{nil};
    std::unordered_map<std::string, double> conceptMemory_;

    Impl() {
        // Initialize NLTagger with required schemes
        NSArray<NLTagScheme>* schemes = @[
            NLTagSchemeNameTypeOrLexicalClass,
            NLTagSchemeNameType
        ];
        tagger_ = [[NLTagger alloc] initWithTagSchemes:schemes];
    }

    ~Impl() {
        tagger_ = nil;
    }
};

SubtitleNLPProcessor::SubtitleNLPProcessor() {
    impl_ = new Impl();
}

SubtitleNLPProcessor::~SubtitleNLPProcessor() {
    delete impl_;
}

std::vector<NLPSecondFeatures> SubtitleNLPProcessor::process(
    const std::vector<arcscope::bridge::SubtitleCue>& cues,
    double durationSeconds
) {
    const size_t sampleCount = static_cast<size_t>(std::ceil(durationSeconds));
    std::vector<NLPSecondFeatures> result(sampleCount);
    std::vector<double> weightSum(sampleCount, 0.0);  // Track overlap weights for averaging

    // Initialize all samples with timeSeconds
    for (size_t i = 0; i < sampleCount; ++i) {
        result[i].timeSeconds = static_cast<double>(i) + 0.5;
    }

    // If no cues, return empty features (caller will use fallback)
    if (cues.empty()) {
        return result;
    }

    // Reset concept memory for this video
    impl_->conceptMemory_.clear();

    // Process each subtitle cue and distribute to second buckets
    for (const auto& cue : cues) {
        // Convert std::string back to NSString for NLP processing
        if (cue.text_utf8.empty()) {
            continue;
        }
        NSString* text = [NSString stringWithUTF8String:cue.text_utf8.c_str()];
        if (!text || text.length == 0) {
            continue;
        }

        // Extract NLP features from this cue
        auto tokens = tokenizeText(text);
        if (tokens.empty()) {
            continue;
        }

        int trueWordCount = static_cast<int>(tokens.size());
        double complexity = calculateSentenceComplexity(text, tokens);
        auto concepts = extractConcepts(text, impl_->tagger_);
        double newConceptRatio = calculateNewConceptRatio(
            concepts,
            impl_->conceptMemory_,
            (cue.start + cue.end) / 2.0  // Use cue midpoint time
        );
        double entityDensity = calculateEntityDensity(text, impl_->tagger_, trueWordCount);

        // Distribute this cue's features to overlapping second buckets
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

            // Accumulate features (weighted by overlap)
            // wordCount is a quantity - accumulate (sum)
            result[sec].wordCount += static_cast<int>(std::round(trueWordCount * ratio));

            // Complexity/ratio metrics - accumulate for weighted average
            weightSum[sec] += ratio;
            result[sec].sentenceComplexity += complexity * ratio;
            result[sec].newConceptRatio += newConceptRatio * ratio;
            result[sec].entityDensity += entityDensity * ratio;
        }
    }

    // WEIGHTED AVERAGING: Normalize ratio metrics by dividing by weight sum
    // This ensures proper averaging when multiple cues overlap the same second
    for (size_t sec = 0; sec < sampleCount; ++sec) {
        if (weightSum[sec] > 1e-9) {
            result[sec].sentenceComplexity /= weightSum[sec];
            result[sec].newConceptRatio /= weightSum[sec];
            result[sec].entityDensity /= weightSum[sec];
        }

        // Clamp to [0, 1] after normalization
        result[sec].sentenceComplexity = std::clamp(result[sec].sentenceComplexity, 0.0, 1.0);
        result[sec].newConceptRatio = std::clamp(result[sec].newConceptRatio, 0.0, 1.0);
        result[sec].entityDensity = std::clamp(result[sec].entityDensity, 0.0, 1.0);
    }

    return result;
}

}  // namespace nlp
}  // namespace arcscope
