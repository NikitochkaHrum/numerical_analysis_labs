#pragma once

#include "function_tool.h"

// Трехточечная квадратурная формула Гаусса
ld gauss(FunctionTool& f, ld a, ld b, ld eps, ld exact_value);
