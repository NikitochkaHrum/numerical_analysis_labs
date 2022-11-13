#include "function_tool.h"

FunctionTool::FunctionTool(function<ld(ld)> _func) {
	func = _func;
}
ld FunctionTool::operator() (ld x)
{
	count++;
	return func(x);
}
void FunctionTool::reset() {
	count = 0;
}
unsigned int FunctionTool::get_calls_count()
{
	return count;
}