#ifndef _UNIFORM_FUNC_H_
#define _UNIFORM_FUNC_H_
#include <iostream>
#include <sstream>

template<typename T>
inline std::string ToString(T value) 
{
	std::stringstream ss;
	ss << value;
	return ss.str();
}

#endif