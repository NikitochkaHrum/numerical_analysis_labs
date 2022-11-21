#pragma once

#include "../matrices.h"

ld solve_nle(const std::function<ld(ld)>& f, ld a, ld b, ld eps); // Метод половинного деления
