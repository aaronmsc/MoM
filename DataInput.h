#ifndef _DATA_INPUT_H_
#define _DATA_INPUT_H_
//#include "PortMatConv.h"
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "Port.h"
#include "Basis.h"
#include "EMOption.h"
#include "Element.h"


namespace mom
{
    class DataInput {
    public:
        DataInput() {}
        bool ReadData(const std::string prefix_path, const std::vector<std::string> &files_name,
            std::vector<Basis> &basis_list, std::vector<Element> &elem_list,
            std::vector<Port *> &port_list);
        //bool PythonFormat(std::vector<OutputParamInfor> infor, int type, int port_num);
    };
}
#endif // _DATA_INPUT_H_
