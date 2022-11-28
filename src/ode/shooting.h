#pragma once

#include "runge_kutt.h"

//метод стрельб
ld shooting_method(function<ld(ld, ld, ld)>& Z, function<ld(ld)>& Y, ld A, ld B, ld x_left, ld x_right, ld h, ld h_eps, ld y_eps);