#pragma once

#include "../function_tool.h"

const ld alpha = 0.5; // Во всех наших задачах h = h / 2 => alpha = 1 / 2
const ld d = 1 / log(alpha); // Константа 1 / ln(alpha)

// Оценка погрешности методом Рунге
ld error_estimation(ld s1, ld s2, ld h, int k); // S(f, h1), S(f, h2), h1, h2, k

// Эмпирическая оценка порядка апроксимации (k)
ld degree_estimation(ld s1, ld s2, ld s3);

void print_header();

void print_iteration(int n, ld h, ld s, ld s1, ld s2, ld s3, int k);
