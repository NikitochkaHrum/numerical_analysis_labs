#pragma once
#include "../matrices.h"

// Вспомогательный класс, для подсчета количества обращений к функции
class FunctionTool
{
private:
	unsigned int count = 0;
	function<ld(ld)>  func;
public:
	FunctionTool(function<ld(ld)>  _func);
	ld operator() (ld x);
	void reset();
	unsigned int get_calls_count();
};