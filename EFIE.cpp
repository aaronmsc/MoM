#include "EFIE.h"

using namespace std;
namespace mom
{
    namespace internal{
        double g_z_GND = 0;
    }

    void EFIE::UnifyGeoUnit()
    {
        // convert vertexes of elements to normal unit.
        for (size_t i = 0; i < element_list->size(); i++) {
            (*element_list)[i].convertLengthUnit(this->m_p_emoption_->unit, M);
        }


        /*for (size_t i=0; i < this->m_metal_layer_.size(); i++) {
            if (m_metal_layer_[i].IsGND()) {
                internal::g_z_GND =
                    EMOption::convertLengthUnit(m_metal_layer_[i].GetLayerZStart(), this->m_p_emoption_->unit, M);
                this->m_p_emoption_->space_type =
                    this->m_p_emoption_->space_type == FREE_SPACE ? GND_FREESPACE : this->m_p_emoption_->space_type;
            }
            m_metal_layer_[i].SetLayer(
                m_metal_layer_[i].GetLayerIndex(),
                EMOption::convertLengthUnit(m_metal_layer_[i].GetLayerZStart(), this->m_p_emoption_->unit, M),
                EMOption::convertLengthUnit(m_metal_layer_[i].GetLayerThickness(), this->m_p_emoption_->unit, M),
                m_metal_layer_[i].GetLayerEMProp(), m_metal_layer_[i].IsGND(), m_metal_layer_[i].IsMetalLayer());
        }*/

        //std::sort(m_metal_layer_.begin(), m_metal_layer_.end(), LayerZStartASC);
    }

    bool PythonFormat2(std::vector<OutputParamInfor> infor, int type, int port_num, std::string out_file)
    {
        int num_param = 0;


        auto ConvertParamType = [](PortMatType type) {
            switch (type) {
            case PortMatType::S: {
                return "S";
            }
            case PortMatType::Y: {
                return "Y";
            }
            case PortMatType::Z: {
                return "Z";
            }
            case PortMatType::UNKOWN_PORT_TYPE: {
                throw std::runtime_error("Unknow PortMatType!");
            }
            default: {
                throw std::runtime_error("Unknow PortMatType!");
            }
            }
        };

        std::ofstream of(out_file);
        if (!of.is_open()) {
            MSG("%s doesn't exist", out_file.c_str());
            return false;
        }

        of << port_num << endl;
        for (int type_flag = 0x1; type_flag != 0x8; type_flag = type_flag << 1) { // mask
            PortMatType type_tmp = static_cast<PortMatType>(type & type_flag);
            if (type_tmp != UNKOWN_PORT_TYPE) {
                for (int m = 1; m <= port_num; m++) {
                    for (int n = 1; n <= port_num; n++) {
                        char param_port[256];
                        sprintf_s(param_port, "%s(%d,%d)", ConvertParamType(type_tmp), m, n);
                        of << param_port << endl;

                        // x
                        for (size_t i = 0; i < infor.size(); i++) {
                            if (infor[i].param_port == param_port)
                                of << infor[i].freq << "	 ";
                        }
                        of << endl;

                        // of << "S_phase=[ ";
                        for (size_t i = 0; i < infor.size(); i++) {
                            if (infor[i].param_port == param_port)
                                of << infor[i].phase << "	 ";
                        }
                        of << endl;

                        // of << "S_real=[ ";
                        for (size_t i = 0; i < infor.size(); i++) {
                            if (infor[i].param_port == param_port)
                                of << infor[i].real_part << "	 ";
                        }
                        of << endl;

                        // of << "S_image=[ ";
                        for (size_t i = 0; i < infor.size(); i++) {
                            if (infor[i].param_port == param_port)
                                of << infor[i].image_part << "	 ";
                        }
                        of << endl;

                        // of << "S_mag=[ ";
                        for (size_t i = 0; i < infor.size(); i++) {
                            if (infor[i].param_port == param_port)
                                of << infor[i].magnitude << "	 ";
                        }
                        of << endl;

                        // of << "S_mag_DB=[ ";
                        for (size_t i = 0; i < infor.size(); i++) {
                            if (infor[i].param_port == param_port)
                                of << infor[i].magnitude_dB << "	 ";
                        }
                        of << endl;
                        of << endl;
                    }
                }
            }
        }

        of.close();
        return true;
    }

    /*void MinGaussDistance(RWG &src, RWG &field, int idx_figure_s, int idx_figure_f, double &min)
    {
        int i, j;
        double distance;
        const std::vector<arma::vec3> &list_gauss_point_s = (src.GetElement(idx_figure_s))->GetGaussPoint();
        const std::vector<arma::vec3> &list_gauss_point_f = (field.GetElement(idx_figure_f))->GetGaussPoint();
        int num_gauss_nodes_s = (src.GetElement(idx_figure_s))->GetGaussPoint().size();
        int num_gauss_nodes_f = (field.GetElement(idx_figure_f))->GetGaussPoint().size();

        for (i = 0; i < num_gauss_nodes_f; i++) {
            for (j = 0; j < num_gauss_nodes_s; j++) {
                if (list_gauss_point_s[i] == list_gauss_point_f[j]) {
                    continue;
                }
                distance = arma::vec3::Distant(list_gauss_point_s[i], list_gauss_point_f[j]);
                min = distance < min ? distance : min;
            }
        }
    }*/

    void EFIE::Run()
    {
        std::map<double, arma::cx_mat> Y_map;
        // std::map<const Port *, int> port_ix_map;
        // std::set<PortGroup *> port_set;
        // std::map<const Port *, const Port *> pair_port_map;
        // std::string prefix_path =
        // 	this->m_p_emoption_->output_file_name.substr(0, m_p_emoption_->output_file_name.find_last_of("/") + 1);
        UnifyGeoUnit();

        MSG("Options: unit-%d", this->m_p_emoption_->unit);
        MSG("Options: space_type-%d", this->m_p_emoption_->space_type);
        MSG("Options: integral_approach-%d", this->m_p_emoption_->integral_approach);
        MSG("Options: method-%d", this->m_p_emoption_->kernel_method);
        //MSG("g_z_GND-%f", internal::g_z_GND);
        MSG("path-%s", this->m_p_emoption_->output_file_name.c_str());

        MSG("number of basis: %d", basis_list->size());


        if (this->m_p_emoption_->integral_approach == RWG_BASIS) {
            srand((unsigned)time(NULL));
            int num_port = 0;

            std::vector<Element> elem_list(*element_list);

            // This is useful
            /*std::ofstream bf("basis.txt");
            bf << basis_list->size() << endl;
            for (size_t i=0; i < basis_list->size(); i++) {
                bf << (*basis_list)[i][0] << " " << (*basis_list)[i][1] << " 1 2 " << endl;
            }
            bf.close();

            std::ofstream ef("element.txt");
            ef << element_list->size() << endl;
            for (size_t i=0; i < this->element_list->size(); i++) {
                ef << "6 0" << endl;
                for (auto vertex : (*element_list)[i].GetVertex()) {
                    ef << vertex[0] << " " << vertex[1] << " " << vertex[2] << endl;
                }
                ef << "0 0 0" << endl;
                ef << "0 0 0" << endl;
                ef << "0 0 0" << endl;
            }
            ef.close();*/



            // basis和element已知，生成rwg_list
            MSG("start Generate RWG");
            for (size_t i = 0; i < basis_list->size(); i++) {
                internal::PortBasisFlag port_basis_flag =
                    RWG::IsUsefulBasis((*basis_list)[i], *element_list, this->m_p_emoption_->feed_type);
                if (port_basis_flag) {
                    RWG rwg((*basis_list)[i], *element_list);
                    rwg.SetPortBasisFlag(port_basis_flag);
                    rwg_list.push_back(rwg);
                }
            }
            MSG("end Generate RWG");

            // pack RWG
            for (size_t i = 0; i < rwg_list.size(); i++) {
                if (!rwg_list[i].IsHalfbasis())
                    rwg_list[i].CreateCompleteRWG(i);
            }

            for (size_t i = 0; i < rwg_list.size(); i++) {
                if (rwg_list[i].IsHalfbasis())
                    rwg_list[i].CreateHalfRWG(i);
            }

            // #ifdef _DEBUG_ENGINE_ALONE_
            // 		double min_gauss_distance = HUGE_VAL;
            // 		int m, n;
            // 		for (size_t i=0; i < rwg_list.size(); i++) {
            // 			RWG &src = rwg_list[i];
            // 			for (int j = i; j < rwg_list.size(); j++) {
            // 				RWG &field = rwg_list[j];
            // 				bool is_single_basis_f = field.IsHalfbasis();
            // 				bool is_single_basis_s = src.IsHalfbasis();
            // 				if (!is_single_basis_f && !is_single_basis_s) { //源和目标都有两个element
            // 					for (m = 0; m < 2; m++)
            // 						for (n = 0; n < 2; n++) {
            // 							MinGaussDistance(src, field, m, n, min_gauss_distance);
            // 						}
            // 				} else if (is_single_basis_f && is_single_basis_s) { //目标和源都只有一个element
            // 					MinGaussDistance(src, field, 0, 0, min_gauss_distance);
            // 				} else if (is_single_basis_f) { //目标只有一个element
            // 					n = 0;
            // 					for (m = 0; m < 2; m++) {
            // 						MinGaussDistance(src, field, m, n, min_gauss_distance);
            // 					}
            // 				} else { //源只有一个element
            // 					m = 0;
            // 					for (n = 0; n < 2; n++) {
            // 						MinGaussDistance(src, field, m, n, min_gauss_distance);
            // 					}
            // 				}
            // 			}
            // 		}
            // 		Msg::Debug("min_gauss_distance=%lf", min_gauss_distance);
            // #endif
            std::vector<PortGroup> port_group_array;

            //将对应端口与port rwg进行关联
            // TODO :add feed_list judgement
            if (this->m_p_emoption_->feed_type == GAP)
                MSG("Feed Type:GAP");
            else if (this->m_p_emoption_->feed_type == VOLTAGE)
                MSG("Feed Type: VOLTAGE");
            else
                MSG("Unknow Feed Type: %d", this->m_p_emoption_->feed_type);
            for (size_t i = 0; i < this->m_p_emoption_->feed_ports.size(); i++) {
                MSG("Simulation port:");
                MSG("No %d. :{%s}", i, this->m_p_emoption_->feed_ports[i].c_str());
                auto port_str = this->m_p_emoption_->feed_ports[i];
                port_group_array.push_back(PortGroup(port_str));
            }
            MSG("rwg_list.size(): %d", rwg_list.size());
            if (rwg_list.size() == 0)
                MSG("rwg list is empty!");
            if (this->m_p_emoption_->feed_type == GAP || this->m_p_emoption_->feed_type == VOLTAGE) {
                const Port *port;
                for (size_t i = 0; i < rwg_list.size(); i++) {
                    port = NULL;
                    if (rwg_list[i].GetElement(0)->GetPort() != NULL) {
                        port = rwg_list[i].GetElement(0)->GetPort();
                    }
                    else if (!rwg_list[i].IsHalfbasis() && rwg_list[i].GetElement(1)->GetPort() != NULL) {
                        port = rwg_list[i].GetElement(1)->GetPort();
                    }
                    if (port != NULL) {
                        for (size_t j = 0; j < port_group_array.size(); j++) {
                            PortGroup *port_group = &port_group_array[j];
                            //cout << i << ":" << port_group->Find(port->GetPortName()) << endl;
                            switch (port_group->Find(port->GetPortName())) {
                            case PortGroup::NEGATIVE:
                                port_group->port_set_neg.insert(port);
                                port_group->basis_list_neg.push_back(PortBasis(i, port, false));
                                break;
                            case PortGroup::POSITIVE:
                                port_group->port_set_pos.insert(port);
                                port_group->basis_list_pos.push_back(PortBasis(i, port, true));
                                break;
                            case PortGroup::UNKNOWN_DIRECTION:
                                break;
                            case PortGroup::BOTH:
                                break;
                            }
                        }
                    }
                }
            }
            int num_differential_port = 0;
            //计算端口数
            for (auto &iter : port_group_array) {
                iter.idx = num_port;
                for (auto port : iter.port_set_pos) {
                    MSG("port[+]{%s}(%d)", port->GetPortName().c_str(), port->GetIndex());
                }
                for (auto port : iter.port_set_neg) {
                    MSG("port[-]{%s}(%d)", port->GetPortName().c_str(), port->GetIndex());
                }
                num_port++;
                for (auto port_basis : iter.basis_list_pos) {

                    const Port *port = port_basis.port;
                    RWG* rwg = &rwg_list[port_basis.basis_ix];
                    if (port != NULL) {
                        if (!rwg->IsHalfbasis() && rwg->GetElement(1)->GetPort() != NULL) {
                            MSG("port[+]{%s} :   \trwg=%d elemP=%d\t   \telemP=%d", port->GetPortName().c_str(),
                                port_basis.basis_ix, rwg->GetElement(0)->GetIndex(),
                                rwg->GetElement(1)->GetIndex());
                        }
                        else if (!rwg->IsHalfbasis()) {
                            MSG("port[+]{%s} :   \trwg=%d elemP=%d\t   \telem=%d", port->GetPortName().c_str(),
                                port_basis.basis_ix, rwg->GetElement(0)->GetIndex(),
                                rwg->GetElement(1)->GetIndex());
                        }
                        else {
                            MSG("port[+]{%s} :   \trwg=%d elemP=%d\t   ", port->GetPortName().c_str(),
                                port_basis.basis_ix, rwg->GetElement(0)->GetIndex());
                        }
                    }
                    else if (!rwg->IsHalfbasis() && rwg->GetElement(1)->GetPort() != NULL) {
                        port = rwg->GetElement(1)->GetPort();
                        MSG("port[+]{%s} :   \trwg=%d elemP=%d\t   \telem=%d", port->GetPortName().c_str(),
                            port_basis.basis_ix, rwg->GetElement(1)->GetIndex(), rwg->GetElement(0)->GetIndex());
                    }
                }
                MSG("port[%d] total:%d:%d", num_port, static_cast<int>(iter.basis_list_pos.size()),
                    static_cast<int>(iter.basis_list_neg.size()));
            }

            MSG("num_port:%d", num_port);
            if (num_port == 0) {
                MSG("num_port == 0!");
                return;
            }
            // num_port = num_port == 0 ? 1 : num_port;
            this->m_port_num_ = num_port;

            // solver begins
            MoMFactory mom_factory;
            unique_ptr<MoM> moment_method = mom_factory.CreateMoM(&rwg_list, this->m_p_emoption_);
            for (int i = 0; i < num_port; i++) {
                if (this->m_p_emoption_->feed_type == GAP || this->m_p_emoption_->feed_type == VOLTAGE)
                    this->m_F_list_.push_back(arma::cx_vec(rwg_list.size(), arma::fill::zeros));
            }



            //对不同端口加激励
            for (auto &iter : port_group_array) {
                // for (int feed_port = 0; feed_port < num_port; feed_port++) {
                //生成激励矩阵
                if (this->m_p_emoption_->feed_type == GAP)
                    CreateGapFeedMatrixRWG(rwg_list, iter, this->m_F_list_[iter.idx]);
                else if (this->m_p_emoption_->feed_type == VOLTAGE)
                    CreateVFeedMatrixRWG(rwg_list, iter, this->m_F_list_[iter.idx]);
#ifdef _DEBUG_
                this->m_F_list_[iter.idx].save(this->m_p_emoption_->output_file_name + "_F_mat" + std::to_string(iter.idx) +
                    ".txt", arma::raw_ascii);
#endif
            }
            MSG("Create F matrix resolved.");
            for (auto freq : this->m_p_emoption_->frequency) { //对每一个频率进行计算
                MSG(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
                MSG("freq: %f begin:", freq);
                // getchar();
                arma::cx_vec J(num_port);			  //由Ii*basis_fn*feed_direction得到的每个端口的J
                arma::cx_mat Y(num_port, num_port); // Y参数
                if (this->m_p_emoption_->feed_type == GAP || this->m_p_emoption_->feed_type == VOLTAGE)
                    this->m_map_I_list_.insert(
                        make_pair(freq, std::vector<arma::cx_vec>(num_port, arma::cx_vec(rwg_list.size()))));
                std::vector<arma::cx_vec> &I_list = this->m_map_I_list_[freq];
                for (int i = 0; i < num_port; i++) {
                    I_list[i].zeros();
                }

                //生成端口的I矩阵
                // Solver
                if (this->m_p_emoption_->kernel_method == MOM) {
                    if (this->m_p_emoption_->feed_type == GAP || this->m_p_emoption_->feed_type == VOLTAGE)
                        moment_method->Run(freq, this->m_F_list_, I_list);
                }
#ifdef _DEBUG_ENGINE_ALONE_
                for (size_t i = 0; i < I_list.size(); i++) {
                    I_list[i].save(this->m_p_emoption_->output_file_name + "_I_list_" + std::to_string(i) + ".txt");
                }
#endif
                // I_mat => Y_param

                for (int feed_port_num = 0; feed_port_num < num_port; feed_port_num++) {
                    J.zeros();
                    //生成每个端口的J参数
                    if (this->m_p_emoption_->feed_type == GAP)
                        CreateJOfGapFeedPortRWG(rwg_list, port_group_array, I_list[feed_port_num], J);
                    else if (this->m_p_emoption_->feed_type == VOLTAGE)
                        CreateJOfVFeedPortRWG(rwg_list, port_group_array, I_list[feed_port_num], J);

                    J = conj(J);
                    //生成Y参数
                    for (int i = 0; i < num_port; i++) {
                        Y(i, feed_port_num) = J[i];
                    }
                }

                MSG("freq: %f finished!", freq);
                MSG("=========================================\n");
                Y_map.insert(std::make_pair(freq, Y));


                int port_num = GetPortNum();
                arma::cx_vec R_diag = arma::cx_vec(port_num);
                R_diag.fill(this->m_p_emoption_->ref_impedance[0]); // this fill with normal impedance
                arma::cx_mat R_mat = diagmat(R_diag);
                this->port_mat_conv.Initialize(PortMatType::Y, port_num, Y_map, R_mat, this->m_p_emoption_);
                this->port_mat_conv.Save(this->m_p_emoption_->output_file_name, this->m_p_emoption_->port_mat_type,
                    this->m_p_emoption_->data_format);
                // Port display
#ifdef _DEBUG_
                PythonFormat2(this->port_mat_conv.GenerateInfor(this->m_p_emoption_->port_mat_type),
                    this->m_p_emoption_->port_mat_type, port_num, this->m_p_emoption_->output_file_name);
#endif
            }


#ifdef _DEBUG_
            PythonFormat2(this->port_mat_conv.GenerateInfor(this->m_p_emoption_->port_mat_type),
                this->m_p_emoption_->port_mat_type, this->GetPortNum(),
                this->m_p_emoption_->output_file_name + std::string("(interpolated)"));
#endif
        }


        MSG("Save Y_map");
        this->m_map_Y_list_ = Y_map;
        // for (auto Y_ : Y_map) {
        // 	std::cout << "freq:" << Y_.first << "\t\t Y(0,0):" << Y_.second(0, 0) << std::endl;
        // }
        MSG("Run finished!");
    }

    //计算端口gap尺寸
    real EFIE::CalculatePortGap(const PortGroup &port_group, const arma::vec3 &direction_V)
    {
        // PortGroup*port_elem0, *port_elem1;
        // real max_gap = 0;

        // // the max projection on feed direction
        // auto max_projection_on_V = [&direction_V](const arma::vec3 &a, const arma::vec3 &b, const arma::vec3 &c) {
        // 	real max_proj_length = 0;
        // 	arma::vec3 v1 = a - b;
        // 	arma::vec3 v2 = b - c;
        // 	arma::vec3 v3 = a - c;

        // 	max_proj_length = std::abs(v1 ^ direction_V) > std::abs(v2 ^ direction_V) ? std::abs(v1 ^ direction_V)
        // 																			  : std::abs(v2 ^ direction_V);
        // 	max_proj_length = std::abs(v3 ^ direction_V) > max_proj_length ? std::abs(v3 ^ direction_V) : max_proj_length;
        // 	return max_proj_length;
        // };

        // for (auto iter = port_map.lower_bound(port_index); iter != port_map.upper_bound(port_index); iter++) {
        // 	port_elem0 = (iter->second->GetElement(0))->GetPort();
        // 	if (port_elem0 != NULL) {
        // 		const std::vector<arma::vec3> &vertex_tri = iter->second->GetElement(0)->GetVertex();
        // 		max_gap = std::max(max_projection_on_V(vertex_tri[0], vertex_tri[1], vertex_tri[2]), max_gap);
        // 	} else if (!iter->second->IsHalfbasis()) {
        // 		port_elem1 = (iter->second->GetElement(1))->GetPort();
        // 		if (port_elem1 != NULL) {
        // 			const std::vector<arma::vec3> &vertex_tri = iter->second->GetElement(1)->GetVertex();
        // 			max_gap = std::max(max_projection_on_V(vertex_tri[0], vertex_tri[1], vertex_tri[2]), max_gap);
        // 		}
        // 	}
        // }
        // return max_gap;
        return 0;
    }



    void EFIE::CreateJOfVFeedPortRWG(const std::vector<RWG> &rwg_list,
        const std::vector<PortGroup> &port_group_array, const arma::cx_vec &I,
        arma::cx_vec &J)
    {
        arma::vec3 basis_fn;
        arma::vec3 feed_direction;
        int port_index;
        const RWG *rwg = NULL;
        //生成各端口的J参数
        for (auto port_group : port_group_array) {
            port_index = port_group.idx;
            for (auto port_basis : port_group.basis_list_pos) {
                rwg = &rwg_list[port_basis.basis_ix];
                if (rwg->GetPortBasisFlag() == internal::EDGE)
                    J[port_index] += I[port_basis.basis_ix];
                else if (rwg->GetPortBasisFlag() == internal::AREA)
                    J[port_index] += I[port_basis.basis_ix] * 2.0;
            }
        }

    }

    //计算端口的J参数
    // Ji = sum(I_list(i,j)*fj*feed_direction) ,for i  = 1,2,...p; j = 1,2,...num_elemt_on_center_line;
    //根据是否与端口中心点在同一水平线上来计算端口的J值
    void EFIE::CreateJOfGapFeedPortRWG(const std::vector<RWG> &rwg_list,
        const std::vector<PortGroup> &port_group_array, const arma::cx_vec &I,
        arma::cx_vec &J)
    {
        // Element *pElement;
        // int elem_subscript;
        // arma::vec3 basis_fn;
        // arma::vec3 feed_direction;
        // const Port *port0;
        // const Port *port1;
        // int port_index;
        // double length;
        // double error = this->m_p_emoption_->unit == MIL
        // 				   ? 2.54e-10
        // 				   : 1e-5 * std::pow(1e-3, double(3 - (int)this->m_p_emoption_->unit));

        // //生成各端口的J参数
        // for (auto iter = port_map.begin(); iter != port_map.end();) {
        // 	auto end_it = port_map.upper_bound(iter->first);
        // 	if (iter->second->IsHalfbasis())
        // 		port0 = (iter->second->GetElement(0))->GetPort();
        // 	else {
        // 		port0 = (iter->second->GetElement(0))->GetPort();
        // 		port1 = (iter->second->GetElement(1))->GetPort();
        // 	}
        // 	const Port *feed = NULL;
        // 	if (port0) {
        // 		feed = port0;
        // 		feed_direction = port0->GetDirection();
        // 	} else if (port1) {
        // 		feed = port1;
        // 		feed_direction = port1->GetDirection();
        // 	}
        // 	// EMOption::convertGeoUnit(this->m_p_emoption_->unit,M,feed_direction);
        // 	port_index = iter->first;

        // 	arma::vec3 port_center = feed->GetBasPnt();
        // 	EMOption::convertGeoUnit(this->m_p_emoption_->unit, M, port_center);

        // 	auto FeedVertex = [&port_center, &feed_direction, &error](Element *pElement, int origin_point) {
        // 		std::vector<arma::vec3> vertex_tri = pElement->GetVertex();

        // 		arma::vec3 mid = (vertex_tri[(origin_point + 1) % 3] + vertex_tri[(origin_point + 2) % 3]) / 2;

        // 		arma::vec3 center_common(arma::vec3(mid[0], mid[1]));
        // 		if (std::abs((center_common - port_center) ^ feed_direction) < error) {
        // 			std::vector<arma::vec3> vertex_common;
        // 			vertex_common.push_back(vertex_tri[(origin_point + 1) % 3]);
        // 			vertex_common.push_back(vertex_tri[(origin_point + 2) % 3]);
        // 			return vertex_common;
        // 		} else {
        // 			return std::vector<arma::vec3>();
        // 		}
        // 	};

        // 	while (iter != end_it) {
        // 		RWG *rwg = iter->second;
        // 		length = 0;

        // 		elem_subscript = 2;
        // 		if (rwg->GetElement(0)->GetPort() != feed) {
        // 			if (rwg->GetElement(1)->GetPort() != feed) {
        // 				iter++;
        // 				continue;
        // 			}
        // 			elem_subscript = 1;
        // 		} else {
        // 			elem_subscript = 0;
        // 		}
        // 		std::vector<arma::vec3> vertex_common =
        // 			FeedVertex(rwg->GetElement(elem_subscript), rwg->GetOriginalPoint(elem_subscript));
        // 		if (vertex_common.size() == 0) {
        // 			iter++;
        // 			continue;
        // 		}
        // 		pElement = rwg->GetElement(elem_subscript);
        // 		std::vector<arma::vec3> vertex_tri = pElement->GetVertex();

        // 		arma::vec3 rho =
        // 			((vertex_common[0] + vertex_common[1]) / 2 - vertex_tri[rwg->GetOriginalPoint(elem_subscript)]) *
        // 			rwg->GetCurrentDirection(elem_subscript);
        // 		basis_fn = rho * (arma::vec3::Distant(vertex_common[0], vertex_common[1])) / rwg->GetJacobi(elem_subscript);
        // 		J[port_index] += (basis_fn ^ feed_direction) * I[rwg - &rwg_list[0]];
        // 		// cout << rwg - &rwg_list[0] << endl;

        // 		iter++;
        // 	}
        // }
    }



    //生成F矩阵
    void EFIE::CreateVFeedMatrixRWG(const std::vector<RWG> &rwg_list, const PortGroup &feed_port,
        arma::cx_vec &F)
    {
        FComplex _p4ie = pi4 * FComplex(0, 1) / _eta;
        const RWG *rwg = NULL;
        for (auto pos_basis : feed_port.basis_list_pos) {
            rwg = &rwg_list[pos_basis.basis_ix];
            int no_rwg = pos_basis.basis_ix;
            if (rwg->IsHalfbasis()) {
                if (rwg->GetPortBasisFlag() == internal::EDGE)
                    F[no_rwg] = 1.0 * _p4ie;
                else if (rwg->GetPortBasisFlag() == internal::AREA)
                    F[no_rwg] = 2.0 * _p4ie;

            }
            else {
                F[no_rwg] = 2.0 * _p4ie;
            }
        }
    }




    void EFIE::CreateGapFeedMatrixRWG(const std::vector<RWG> &rwg_list, const PortGroup &feed_port_name,
        arma::cx_vec &F)
    {
        // const Port *port_elem0, *port_elem1;
        // std::map<int, real> gap_map;

        // auto Fq_integral = [rwg_list](const int no_rwg, const int no_elem, const arma::vec3 &E_in) {
        // 	arma::vec3 gauss_point, current_gauss_point;
        // 	FComplex integral_result;
        // 	Element *elem_ptr = rwg_list[no_rwg].GetElement(no_elem);
        // 	const std::vector<arma::vec3> &list_current_gauss_point = rwg_list[no_rwg].GetRho(no_elem);
        // 	const std::vector<real> &list_gauss_weight = elem_ptr->GetGaussWeight();
        // 	int num_gauss_point = elem_ptr->GetGaussPoint().size();
        // 	for (size_t i=0; i < num_gauss_point; i++) {
        // 		gauss_point = elem_ptr->GetGaussPoint()[i];
        // 		current_gauss_point = list_current_gauss_point[i];

        // 		integral_result += FComplex(current_gauss_point ^ E_in) * list_gauss_weight[i];
        // 	}
        // 	return integral_result;
        // };
        // // wrong here :FComplex k0 = 2 * PI/ (c_wave_speed / freq);
        // FComplex _p4ie = pi4 * FComplex(0, 1) / _eta;

        // for (auto iter = port_map.begin(); iter != port_map.end();) {
        // 	port_elem0 = iter->second->GetElement(0)->GetPort();
        // 	if (!iter->second->IsHalfbasis())
        // 		port_elem1 = iter->second->GetElement(1)->GetPort();
        // 	else
        // 		port_elem1 = NULL;
        // 	const Port *feed = port_elem0 != NULL ? port_elem0 : port_elem1;
        // 	arma::vec3 port_center = feed->GetBasPnt();
        // 	arma::vec3 port_vec = feed->GetLocalPnt();
        // 	EMOption::convertGeoUnit(this->m_p_emoption_->unit, M, port_center);
        // 	EMOption::convertGeoUnit(this->m_p_emoption_->unit, M, port_vec);
        // 	if (port_elem0 && port_elem0->GetPortName() == feed_port_name) { //只对当前激励端口计算gap
        // 		gap_map.insert(std::make_pair(iter->first, arma::vec3::Distant(port_center, port_vec)));
        // 	} else if (port_elem1 && port_elem1->GetPortName() == feed_port_name)
        // 		gap_map.insert(std::make_pair(iter->first, arma::vec3::Distant(port_center, port_vec)));
        // 	iter = port_map.upper_bound(iter->first);
        // }

        // for (auto iter = gap_map.begin(); iter != gap_map.end(); iter++) {
        // 	Msg::Debug("port gap:%lf", iter->second);
        // }

        // for (auto iter = port_map.begin(); iter != port_map.end();) {
        // 	auto end_it = port_map.upper_bound(iter->first);

        // 	while (iter != end_it) {
        // 		int no_rwg = iter->second - &rwg_list[0];

        // 		port_elem0 = iter->second->GetElement(0)->GetPort();
        // 		if (port_elem0 != NULL && port_elem0->GetIndex() == feed_port) {
        // 			real gap_size = gap_map.find(iter->first)->second;
        // 			F[no_rwg] = Fq_integral(no_rwg, 0, port_elem0->GetDirection() / gap_size) * _p4ie;
        // 		}
        // 		if (!iter->second->IsHalfbasis()) {
        // 			port_elem1 = iter->second->GetElement(1)->GetPort();
        // 			if (port_elem1 != NULL && port_elem1->GetIndex() == feed_port) {
        // 				real gap_size = gap_map.find(iter->first)->second;
        // 				F[no_rwg] += Fq_integral(no_rwg, 1, port_elem1->GetDirection() / gap_size) * _p4ie;
        // 			}
        // 		}
        // 		iter++;
        // 	}
        // }
    }

    //生成F矩阵




}
