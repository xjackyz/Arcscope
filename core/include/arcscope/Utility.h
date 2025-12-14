#pragma once

#include "arcscope/FilmTypes.h"

#include <string>

namespace arcscope {

std::string curve_kind_to_string(CurveKind kind);
CurveKind curve_kind_from_string(const std::string& value);

}

