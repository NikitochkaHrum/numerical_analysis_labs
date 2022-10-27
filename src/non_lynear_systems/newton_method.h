#pragma once

#include "../matrices.h"

Vec newton_method(Vec init, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2,
                        function<ld(ld, ld)> der_f1_x, function<ld(ld, ld)> der_f1_y, function<ld(ld, ld)> der_f2_x, function<ld(ld, ld)> der_f2_y,
                        ld eps);

Mat get_matrix_of_inversed_jacobian(function<ld(ld, ld)> der_f1_x, function<ld(ld, ld)> der_f1_y, function<ld(ld, ld)> der_f2_x, function<ld(ld, ld)> der_f2_y,
                                        Vec vals);
