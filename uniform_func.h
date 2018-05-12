#ifndef _UNIFORM_FUNC_H_
#define _UNIFORM_FUNC_H_
#include <iostream>
#include <sstream>

namespace mom {
#define _DEBUG_
    template<typename T>
    inline std::string ToString(T value)
    {
        std::stringstream ss;
        ss << value;
        return ss.str();
    }
    inline bool SkipComment(std::fstream &in, std::stringstream &ss)
    {
        std::string line;
        ss.clear();
        std::getline(in, line);
        ss << line.substr(0, line.find("#"));
        if (ss.bad() || in.bad() || ss.fail()) {
            return false;
        }
        return true;
    }

    inline int GetInt(std::stringstream &ss)
    {
        int a;
        ss >> a;
        return a;
    }
    inline int MSG(const char* format, ...) {
        int _Result;
        va_list _ArgList;
        char _Format[1000];
        strcpy_s(_Format, format);
        strcat_s(_Format, "\n");
        __crt_va_start(_ArgList, format);
        _Result = _vfprintf_l(stdout, _Format, NULL, _ArgList);
        __crt_va_end(_ArgList);
        return _Result;
    }
}
#endif
