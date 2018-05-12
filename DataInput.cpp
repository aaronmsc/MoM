#include "DataInput.h"

namespace mom
{
#define _DEBUG_RWG_
#define _DEBUG_VOLTAGE_



    //bool PythonFormat(std::vector<OutputParamInfor> infor, int type, int port_num, std::string out_file)
    //{
    //    int num_param = 0;
    //    // int port_num_square = port_num * port_num;
    //    // int group_num = infor.size() / num_param / port_num_square;
    //    // if (num_param == 0) {
    //    // 	return false;
    //    // }

    //    auto ConvertParamType = [](PortMatType type) {
    //        switch (type) {
    //        case PortMatType::S: {
    //            return "S";
    //        }
    //        case PortMatType::Y: {
    //            return "Y";
    //        }
    //        case PortMatType::Z: {
    //            return "Z";
    //        }
    //        case PortMatType::UNKOWN_PORT_TYPE: {
    //            throw runtime_error("Unknow PortMatType!");
    //        }
    //        default: {
    //            throw runtime_error("Unknow PortMatType!");
    //        }
    //        }
    //    };

    //    ofstream of(out_file);
    //    if (!of.is_open()) {
    //        MSG("%s doesn't exist", out_file.c_str());
    //        return false;
    //    }
    //    of << port_num << endl;
    //    for (int type_flag = 0x1; type_flag != 0x8; type_flag = type_flag << 1) { // mask
    //        PortMatType type_tmp = static_cast<PortMatType>(type & type_flag);
    //        if (type_tmp != UNKOWN_PORT_TYPE) {
    //            for (int m = 1; m <= port_num; m++) {
    //                for (int n = 1; n <= port_num; n++) {
    //                    char param_port[256];
    //                    sprintf_s(param_port, "%s(%d,%d)", ConvertParamType(type_tmp), m, n);
    //                    of << param_port << endl;

    //                    // x
    //                    for (size_t i=0; i < infor.size(); i++) {
    //                        if (infor[i].param_port == param_port)
    //                            of << infor[i].freq << " ";
    //                    }
    //                    of << endl;
    //                    // of << "S_phase=[ ";
    //                    for (size_t i=0; i < infor.size(); i++) {
    //                        if (infor[i].param_port == param_port)
    //                            of << infor[i].phase << " ";
    //                    }
    //                    of << endl;

    //                    // of << "S_real=[ ";
    //                    for (size_t i=0; i < infor.size(); i++) {
    //                        if (infor[i].param_port == param_port)
    //                            of << infor[i].real_part << " ";
    //                    }
    //                    of << endl;

    //                    // of << "S_image=[ ";
    //                    for (size_t i=0; i < infor.size(); i++) {
    //                        if (infor[i].param_port == param_port)
    //                            of << infor[i].image_part << " ";
    //                    }
    //                    of << endl;

    //                    // of << "S_mag=[ ";
    //                    for (size_t i=0; i < infor.size(); i++) {
    //                        if (infor[i].param_port == param_port)
    //                            of << infor[i].magnitude << " ";
    //                    }
    //                    of << endl;

    //                    // of << "S_mag_DB=[ ";
    //                    for (size_t i=0; i < infor.size(); i++) {
    //                        if (infor[i].param_port == param_port)
    //                            of << infor[i].magnitude_dB << " ";
    //                    }
    //                    of << endl;
    //                    of << endl;
    //                }
    //            }
    //        }
    //    }

    //    of.close();
    //    return true;
    //}


    std::string& Trim(std::string& s) {
        if (s.empty()) {
            return s;
        }

        s.erase(0, s.find_first_of(" "));
        s.erase(s.find_last_not_of(" ") + 1);
        return s;
    }




    bool DataInput::ReadData(const std::string prefix_path, const std::vector<std::string> &files_name,
        std::vector<Basis> &basis_list, std::vector<Element> &elem_list,
        std::vector<Port *> &port_list)
    {
        std::stringstream ss;
        std::string ignore;
        auto DotInPolygon = [](const std::vector<arma::vec3> &figure, const arma::vec3 &point, const CoordinatePlane &plane) {
            bool oddNode = false;
            double xp;
            double eps = 1e-8;
            double x = 0, y = 0;
            switch (plane) {
            case XY_PLANE:
                x = point[0];
                y = point[1];
                break;
            case YZ_PLANE:
                x = point[1];
                y = point[2];
                break;
            case XZ_PLANE:
                x = point[0];
                y = point[2];
                break;
            }

            for (size_t i = 0, j = figure.size() - 1; i < figure.size(); j = i, i++) {
                if ((figure[i][1] < y && figure[j][1] >= y || figure[j][1] < y && figure[i][1] >= y) &&
                    (figure[i][0] <= x || figure[j][0] <= x)) {
                    xp = (y - figure[i][1]) * (figure[j][0] - figure[i][0]) / (figure[j][1] - figure[i][1]) + figure[i][0];
                    if (std::abs(xp - x) < eps) {
                        return true;
                    }
                    else if (xp < x) {
                        oddNode = !oddNode;
                    }
                }
                else if ((figure[i][0] <= x && figure[j][0] >= x || figure[j][0] <= x && figure[i][0] >= x) &&
                    (std::abs(figure[i][1] - y) < eps || std::abs(figure[j][1] - y) < eps)) {
                    return true;
                }
            }
            return oddNode;
        };


        for (auto file_name : files_name) {
            if (file_name == "rwg.txt") {
                int num = 0, basis_left = 0, basis_right = 0;
                std::fstream basis_in(prefix_path + "rwg.txt");
                if (!basis_in.is_open()) {
                    MSG("file doesn't exist");
                    system("pause");
                    exit(1);
                }

                SkipComment(basis_in, ss);
                ss >> num;
                std::string k, l;
                while (num--) {
                    SkipComment(basis_in, ss);
                    ss >> basis_left >> basis_right >> ignore >> ignore;
                    basis_list.push_back(Basis(basis_left, basis_right));
                }
                basis_in.close();
            }
            else if (file_name == "tri.txt") {
                int num = 0;
                double x, y, z;
                std::vector<arma::vec3> vertex_list;
                std::fstream elem_in(prefix_path + "tri.txt");
                if (!elem_in.is_open()) {
                    MSG("file doesn't exist");
                    system("pause");
                    exit(1);
                }

                SkipComment(elem_in, ss);
                ss >> num;
                for (int i = 0; i < num; i++) {
                    SkipComment(elem_in, ss);
                    ss >> ignore >> ignore;

                    for (int j = 0; j < 3; j++) {
                        SkipComment(elem_in, ss);
                        ss >> x >> y >> z;
                        vertex_list.push_back(arma::vec3({ x, y, z }));
                    }

                    for (int j = 0; j < 3; j++) {
                        SkipComment(elem_in, ss);
                        ss >> ignore >> ignore >> ignore;
                    }
                    arma::vec3 mid = (vertex_list[0] + vertex_list[1] + vertex_list[2]) / 3;
                    elem_list.push_back(Element(i, vertex_list));
                    vertex_list.clear();
                }
                elem_in.close();
            }
            else if (file_name == "feed.txt") {
                int num = 0;
                int index;
                std::string feed_direction;
                std::string port_name;
                std::string projection_plane;
                real x, y;
                arma::vec3 bas_point, direction_point;
                std::fstream feed_in(prefix_path + "feed.txt");
                if (!feed_in.is_open()) {
                    MSG("file doesn't exist");
                    system("pause");
                    exit(1);
                }
                SkipComment(feed_in, ss);
                ss >> num;
                SkipComment(feed_in, ss);
                ss >> projection_plane;

                CoordinatePlane plane = XY_PLANE;
                if (projection_plane == "XY")
                    plane = XY_PLANE;
                else if (projection_plane == "YZ")
                    plane = YZ_PLANE;
                else if (projection_plane == "XZ")
                    plane = XZ_PLANE;
                while (num--) {

                    std::vector<arma::vec3> port_rect;
                    SkipComment(feed_in, ss);
                    ss >> port_name;
                    SkipComment(feed_in, ss);
                    ss >> index;
                    for (size_t i = 0; i < 4; i++) {
                        SkipComment(feed_in, ss);
                        ss >> x >> y;
                        arma::vec3 p({ x, y ,0 });
                        port_rect.push_back(arma::vec3({ x, y ,0 }));
                    }

                    Port *port_tmp = new Port();
                    // if (feed_direction == "+")
                    // 	port_tmp->SetPortPolarity(true);
                    // if (feed_direction == "-")
                    // 	port_tmp->SetPortPolarity(false);
                    port_tmp->InitPortByRect(port_rect, index, port_name);

                    port_list.push_back(port_tmp);
                }
                /*for(auto port:port_list){
                std::vector<arma::vec3> port_rect;
                port->GetVecPoints(&port_rect);
                for(auto point:port_rect){
                point.Show();
                }
                }*/
                int count = 0;
                for (auto &elem : elem_list) {
                    std::vector<arma::vec3> vertex = elem.GetVertex();
                    std::vector<arma::vec3> mid_vertex;

                    arma::vec3 mid1 = (vertex[0] + vertex[1]) / 2;
                    arma::vec3 mid2 = (vertex[1] + vertex[2]) / 2;
                    arma::vec3 mid3 = (vertex[0] + vertex[2]) / 2;
                    arma::vec3 mid = (vertex[0] + vertex[1] + vertex[2]) / 3;

                    for (auto &port : port_list) {
                        std::vector<arma::vec3> port_rect;
                        port_rect = port->GetVecPoints();
                        // if (PolygonInteraction(vertex, port_rect, plane)) {

                        if ((DotInPolygon(port_rect, mid1, plane)) || (DotInPolygon(port_rect, mid2, plane)) ||
                            (DotInPolygon(port_rect, mid3, plane))) {
                            elem.SetPort(static_cast<const Port*>(port));
                            // cout << "port elem index:" << elem.GetIndex() << " at " << port->GetPortName() << endl;
                        }
                    }
                }
            }
        }
        return true;
    }
}
