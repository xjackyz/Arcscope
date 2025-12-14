#include "arcscope/Utility.h"

#include <stdexcept>

namespace arcscope {

std::string curve_kind_to_string(CurveKind kind) {
    switch (kind) {
        case CurveKind::Pace:
            return "pace";
        case CurveKind::Sound:
            return "sound";
        case CurveKind::Color:
            return "color";
        case CurveKind::Info:
            return "info";
        case CurveKind::Arousal:
            return "arousal";
        case CurveKind::FaceAffect:
            return "face_affect";
        case CurveKind::FacePresence:
            return "face_presence";
    }
    return "pace";
}

CurveKind curve_kind_from_string(const std::string& value) {
    if (value == "pace") {
        return CurveKind::Pace;
    }
    if (value == "sound") {
        return CurveKind::Sound;
    }
    if (value == "color") {
        return CurveKind::Color;
    }
    if (value == "info") {
        return CurveKind::Info;
    }
    if (value == "arousal") {
        return CurveKind::Arousal;
    }
    if (value == "face_affect") {
        return CurveKind::FaceAffect;
    }
    if (value == "face_presence") {
        return CurveKind::FacePresence;
    }
    throw std::runtime_error("Unknown curve kind: " + value);
}

}  // namespace arcscope
