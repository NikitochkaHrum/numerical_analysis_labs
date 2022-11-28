#include "shooting.h"
#include "method_log/method_info.h"

ld shooting_method(function<ld(ld, ld, ld)>& Z, function<ld(ld)>& Y, ld A, ld B, ld x_left, ld x_right, ld h, ld h_eps, ld y_eps){
    log_shooting_method_header();
	ld z_l = 0, z_r = 3, z_m;
	ld y_right;
	int i = 1;

	do
	{
		// Метод половинного деления
		z_m = (z_l + z_r) / 2;
		
		// Получим значения с помощью метода РК 4-го порядка
		auto values = runge_kutta_method(Z, Y, x_left, x_right, A, z_m, h, h_eps, false);
		y_right = get<1>(values.back());

		if (y_right < B)
			z_l = z_m;
		else
			z_r = z_m;

		log_shooting_method_iteration(i, z_m, y_right, B);
		i++;
	} while (abs(y_right - B) > y_eps); // Продолжаем сужать отрезок пока не добьемся нужно точности

	return z_m;
}