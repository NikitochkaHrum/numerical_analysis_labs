#pragma once

#include "../matrices.h"

tuple<ld, ld, Vec> calc_params(Vec* x, Vec* x_next, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2);
tuple<ld, ld, ld, Vec> calc_params(Vec* x, Vec* x_next, Vec* solution, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2);