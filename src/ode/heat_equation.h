#pragma once

#include "../matrices.h"

void heat_eq_explicit_schema(function<ld(ld, ld)> &f, ld a, ld b, function<ld(ld)> &phi, ld hi, int n, ld end_t, function<ld(ld, ld)> &exact_solution);

void heat_eq_implicit_schema(function<ld(ld, ld)> &f, ld a, ld b, function<ld(ld)> &phi, ld hi, int n, ld end_t, function<ld(ld, ld)> &exact_solution);
