#pragma once
#include "../matrices.h"

Vec simple_iteration(Vec init, Vec solution, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2,
                        function<ld(ld, ld)> fi1, function<ld(ld, ld)> fi2,
                        function<ld(ld, ld)> der_fi1_x, function<ld(ld, ld)> der_fi1_y,
                        function<ld(ld, ld)> der_fi2_x, function<ld(ld, ld)> der_fi2_y,
                        ld eps);