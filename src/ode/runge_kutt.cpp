#include "runge_kutt.h"
#include "method_log/method_info.h"

// Метод Рунге-Кутта 4-го порядка
pair<ld, ld> runge_kutta_step(function<ld(ld, ld, ld)> &Z, ld h, ld x, ld y, ld z){
    vector<ld> k_y(4, 0), k_z(4, 0);
	// для y: y'(x) = z(x) ->
	// k0_y(h) = h * z
	// k1_y(h) = h * (z + k0_z / 2)
	// k2_y(h) = h * (z + k1_z / 2)
	// k3_y(h) = h * (z + k2_z)

	// для z: z'(x) = Z(x, y(x), z(x))
	// k0_z(h) = h * Z(x, y, z)
	// k1_z(h) = h * Z(x + h/2, y + k0_y / 2, z + k0_z / 2)
	// k2_z(h) = h * Z(x + h/2, y + k1_y / 2, z + k1_z / 2)
	// k3_z(h) = h * Z(x + h, y + k2_y, z + k2_z)

	k_y[0] = h * z;
	k_z[0] = h * Z(x, y, z);

	k_y[1] = h * (z + k_z[0] / 2.);
	k_z[1] = h * Z(x + h / 2., y + k_y[0] / 2., z + k_z[0] / 2.);

	k_y[2] = h * (z + k_z[1] / 2.);
	k_z[2] = h * Z(x + h / 2., y + k_y[1] / 2., z + k_z[1] / 2.);

	k_y[3] = h * (z + k_z[2]);
	k_z[3] = h * Z(x + h, y + k_y[2], z + k_z[2]);

	ld y_next = y + (k_y[0] + 2. * k_y[1] + 2. * k_y[2] + k_y[3]) / 6.;
	ld z_next = z + (k_z[0] + 2. * k_z[1] + 2. * k_z[2] + k_z[3]) / 6.;

	return { y_next, z_next };
}
vector<triplet> runge_kutta_method(function<ld(ld, ld, ld)> &Z, function<ld(ld)> &Y_r, ld x_first, ld x_last, ld y, ld z, ld h, ld EPS, bool log){
    // x = 0
	// y = 1
	// z = alpha

    if(log)
		log_runge_kutta_header();
	vector<triplet> result;
	pair<ld, ld> left_values, center_values, right_values, next_values;
	ld _h = h;

	for (ld x = x_first; x < x_last; x += h) {
		h = _h;
		// Автоматический контроль точности
		do
		{
			left_values = runge_kutta_step(Z, h, x, y, z);
			center_values = runge_kutta_step(Z, h / 2., x, y, z);
			right_values = runge_kutta_step(Z, h / 2., x, center_values.first, center_values.second);
			
			h /= 2.;
		} while (abs(left_values.first - right_values.first) > EPS); // продолжаем уменьшать отрезок пока не добьемс¤ нужной точности
		
		h *= 2;
		result.push_back({ x, y, z });
		if(log)
			log_runge_kutta_iteration(x, y, Y_r(x), z); // выведем информацию о итерации
		// обновим новыми значени¤ми

		y = left_values.first;
		z = left_values.second;
	}

	return result;
}