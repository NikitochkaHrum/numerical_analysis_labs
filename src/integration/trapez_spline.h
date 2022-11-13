#pragma once

#include "function_tool.h"

ld trapezoidal_rule_spline(FunctionTool & f, const function<ld(ld)> & der_f, ld a, ld b, ld eps, ld exact_value);
