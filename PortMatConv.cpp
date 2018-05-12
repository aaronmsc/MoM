#include "PortMatConv.h"

namespace mom
{
    void PortMatConv::Initialize(const PortMatType in_param, const int num_port, const std::map<double, arma::cx_mat> &param_map,
        const arma::cx_mat &Z_in, const EMOption *pOption)
    {
        this->m_in_param_ = in_param;
        this->m_num_port_ = num_port;
        this->m_Z_in_ = Z_in;
        this->m_p_emoption_ = pOption;
        // get G mat
        auto f = [](double &x) { x = 1 / sqrt(x); };
        ;
        arma::cx_mat tmp(arma::diagmat(arma::vec(arma::abs(arma::real(m_Z_in_.diag()))).for_each(f)), arma::mat(arma::size(Z_in)));
        this->m_G_in_ = tmp;

        this->m_conj_Z_in_ = conj(Z_in);
        this->m_inv_G_in_ = inv(this->m_G_in_);

        // add one param_map to m_in_;
        for (auto value : param_map) {
            PortMatSet store_mat;

            store_mat[in_param] = value.second;
            this->m_in_.insert(std::make_pair(value.first, store_mat));
        }
    }



    //生成输出给Python端的数据
    std::vector<OutputParamInfor> PortMatConv::GenerateInfor(PortMatType type)
    {
        std::vector<OutputParamInfor> output_list;
        char temp_str[256];
        for (auto iter : this->m_in_) {
            for (int type_flag = 0x1; type_flag != 0x8; type_flag = type_flag << 1) { // mask
                PortMatType type_tmp = static_cast<PortMatType>(type & type_flag);
                if (type_tmp != UNKOWN_PORT_TYPE) {
                    arma::cx_mat &output_mat = this->at(iter.first, type_tmp);
                    for (int i = 0; i < this->m_num_port_; i++) {
                        for (int j = 0; j < this->m_num_port_; j++) {
                            OutputParamInfor infor;
                            std::complex<double> elem = output_mat(i, j);
                            infor.freq = iter.first;
                            switch (type_tmp) {
                            case S:
                                infor.param_type = "S";
                                sprintf_s(temp_str, "S(%d,%d)", i + 1, j + 1);
                                infor.param_port = temp_str;
                                break;
                            case Y:
                                infor.param_type = "Y";
                                sprintf_s(temp_str, "Y(%d,%d)", i + 1, j + 1);
                                infor.param_port = temp_str;
                                break;
                            case Z:
                                infor.param_type = "Z";
                                sprintf_s(temp_str, "Z(%d,%d)", i + 1, j + 1);
                                infor.param_port = temp_str;
                                break;
                            }
                            infor.real_part = std::real(elem);
                            infor.image_part = std::imag(elem);
                            infor.magnitude = std::abs(elem);
                            infor.phase = std::arg(elem) / PI * 180;
                            infor.magnitude_dB = 20 * log10(std::abs(elem));

                            output_list.push_back(infor);
                        }
                    }
                }
            }
        }

        return output_list;
    }

    std::vector<OutputParamInfor> PortMatConv::GenerateInfor(int type)
    {
        std::vector<OutputParamInfor> result;
        for (int type_flag = 0x1; type_flag != 0x8; type_flag = type_flag << 1) { // mask
            PortMatType type_tmp = static_cast<PortMatType>(type & type_flag);
            if (type_tmp != UNKOWN_PORT_TYPE) {
                std::vector<OutputParamInfor> vec_tmp = this->GenerateInfor(type_tmp);
                std::copy(vec_tmp.begin(), vec_tmp.end(), back_inserter(result));
            }
        }
        return result;
    }

    bool PortMatConv::Save(const std::string file_name, const PortMatType type, const DataFormat format,
        const FrequencyUnit freq_unit, const int n)
    {
        bool retval = true;
        for (int type_flag = 0x1; type_flag != 0x8; type_flag = type_flag << 1) { // mask
            PortMatType type_tmp = static_cast<PortMatType>(type & type_flag);
            if (type_tmp != UNKOWN_PORT_TYPE) {
                retval &= this->SaveImplement(file_name, type_tmp, format, freq_unit, n);
            }
        }
        return retval;
    }

    //按照选项将SYZ参数保存
    bool PortMatConv::SaveImplement(const std::string file_name, const PortMatType out_param, const DataFormat format,
        const FrequencyUnit freq_unit, const int n)
    {
        if (out_param != S && out_param != Y && out_param != Z)
            throw std::runtime_error("illegal PortMatType param");
        std::stringstream ss;
        switch (out_param) {
        case S:
            ss << file_name << ".s" << this->m_num_port_ << "p";
            break;
        case Y:
            ss << file_name << ".y" << this->m_num_port_ << "p";
            break;
        case Z:
            ss << file_name << ".z" << this->m_num_port_ << "p";
            break;
        }

        std::string s = ss.str();
        std::ofstream of(ss.str().c_str());
        std::string str_val1, str_val2, str_param;

        std::ofstream::fmtflags default_fmt = of.flags();

        if (!of.is_open())
            return false;
        time_t ptime;
        time(&ptime);
        std::string current_time(asctime(localtime(&ptime)));
        of << "! Created " << current_time;
        of << "! List of port names : example Port n = <port_name>" << std::endl;

        int ports_num = this->m_p_emoption_->feed_ports.size();
        for (int i = 0; i < ports_num; i++) {
            of << "! Port " << i + 1 << " = " << this->m_p_emoption_->feed_ports[i] << std::endl;
        }
        of << "# ";
        switch (freq_unit) {
        case HZ:
            of << "HZ";
            break;
        case KHZ:
            of << "KHZ";
            break;
        case MHZ:
            of << "MHZ";
            break;
        case GHZ:
            of << "GHZ";
            break;
        default:
            of << "GHZ";
        }

        of << " ";

        switch (out_param) {
        case S:
            of << "S";
            str_param = "S";
            break;
        case Y:
            of << "Y";
            str_param = "Y";
            break;
        case Z:
            of << "Z";
            str_param = "Z";
            break;
        }

        of << " ";

        switch (format) {
        case MA:
            of << "MA";
            str_val1 = "Mag";
            str_val2 = "Ang";
            break;
        case DB:
            of << "DB";
            str_val1 = "DB";
            str_val2 = "Ang";
            break;
        case RI:
            of << "RI";
            str_val1 = "Re";
            str_val2 = "Im";
            break;
        }

        of << " R ";

        of << n;

        of << std::endl;

        of << "! " << this->m_num_port_ << " Port Netword Data";
        of << std::endl;
        of << "! freq";

        //延迟计算
        for (auto pair : this->m_in_) {
            Apply(pair.first, out_param);
        }

        if (this->m_num_port_ <= 2) {
            for (int i = 1; i <= this->m_num_port_; i++)
                for (int j = 1; j <= this->m_num_port_; j++) {
                    of << "\t" << str_val1 << str_param << i << j;
                    of << "\t" << str_val2 << str_param << i << j;
                }
            of << std::endl;

            of.precision(11);
            of.flags(std::ios::left);

            switch (format) {
            case MA: {
                for (auto value : this->m_in_) {
                    of << std::setw(15) << EMOption::convertFreqUnit(value.first, HZ, freq_unit) << " ";
                    for (int i = 0; i < this->m_num_port_; i++) {
                        for (int j = 0; j < this->m_num_port_; j++) {
                            std::complex<double> tmp = (value.second)[out_param](i, j);
                            of << std::setw(15) << std::abs(tmp) << " " << std::setw(15) << std::arg(tmp) / PI * 180 << " "
                                << "\t";
                        }
                    }
                    of << std::endl;
                }
                break;
            }
            case DB: {
                for (auto value : this->m_in_) {
                    of << std::setw(15) << EMOption::convertFreqUnit(value.first, HZ, freq_unit) << " ";
                    for (int i = 0; i < this->m_num_port_; i++) {
                        for (int j = 0; j < this->m_num_port_; j++) {
                            std::complex<double> tmp = (value.second)[out_param](i, j);
                            of << std::setw(15) << 20 * log10(std::abs(tmp)) << " " << std::setw(15)
                                << std::arg(tmp) / PI * 180 << " "
                                << "\t";
                        }
                    }
                    of << std::endl;
                }
                break;
                break;
            }
            case RI: {
                for (auto value : this->m_in_) {
                    of << std::setw(15) << EMOption::convertFreqUnit(value.first, HZ, freq_unit) << " ";
                    for (int i = 0; i < this->m_num_port_; i++) {
                        for (int j = 0; j < this->m_num_port_; j++) {
                            std::complex<double> tmp = (value.second)[out_param](i, j);
                            of << std::setw(15) << std::real(tmp) << " " << std::setw(15) << std::imag(tmp) << " "
                                << "\t";
                        }
                    }
                    of << std::endl;
                }
                break;
            }
            }

        }
        else if (this->m_num_port_ >= 3) {
            int i, j;

            // step means data pair on every row
            int step = std::min(this->m_num_port_, 4);

            for (i = 1; i <= this->m_num_port_; i++) {
                j = 1;
                if (i != 1) {
                    of << "!     ";
                }
                for (; j <= step; j++) {
                    of << "\t" << str_val1 << str_param << i << j;
                    of << "\t" << str_val2 << str_param << i << j;
                }
                of << " ! row " << i;
                of << std::endl;
                for (; j <= this->m_num_port_; j++) {
                    of << "\t" << str_val1 << str_param << i << j;
                    of << "\t" << str_val2 << str_param << i << j;
                    if (j % step == 0)
                        of << std::endl;
                }
                if ((j - 1) % step != 0)
                    of << std::endl;
            }

            of.precision(11);
            of.flags(std::ios::right);

            switch (format) {
            case MA: {
                int i, j;

                int step = std::min(this->m_num_port_, 4);

                for (auto value : this->m_in_) {
                    of << std::setw(15) << EMOption::convertFreqUnit(value.first, HZ, freq_unit) << " ";
                    for (i = 1; i <= this->m_num_port_; i++) {
                        if (i != 1) {
                            of << std::setw(15) << " ";
                        }
                        for (j = 1; j <= this->m_num_port_; j++) {
                            std::complex<double> tmp = (value.second)[out_param].at(i - 1, j - 1);
                            of << std::setw(15) << std::abs(tmp) << " " << std::setw(15) << std::arg(tmp) / PI * 180 << "\t"
                                << " ";
                            if (j % step == 0)
                                of << std::endl;
                        }
                        if ((j - 1) % step != 0)
                            of << std::endl;
                    }
                }

                break;
            }
            case DB: {
                int i, j;

                int step = std::min(this->m_num_port_, 4);

                for (auto value : this->m_in_) {
                    of << std::setw(15) << EMOption::convertFreqUnit(value.first, HZ, freq_unit) << " ";
                    for (i = 1; i <= this->m_num_port_; i++) {
                        if (i != 1) {
                            of << std::setw(15) << " ";
                        }
                        for (j = 1; j <= this->m_num_port_; j++) {
                            std::complex<double> tmp = (value.second)[out_param].at(i - 1, j - 1);
                            of << std::setw(15) << 20 * log10(std::abs(tmp)) << " " << std::setw(15)
                                << std::arg(tmp) / PI * 180 << "\t"
                                << " ";
                            if (j % step == 0)
                                of << std::endl;
                        }
                        if ((j - 1) % step != 0)
                            of << std::endl;
                    }
                }

                break;
            }
            case RI: {
                int i, j;

                int step = std::min(this->m_num_port_, 4);

                for (auto value : this->m_in_) {
                    of << std::setw(15) << EMOption::convertFreqUnit(value.first, HZ, freq_unit) << " ";
                    for (i = 1; i <= this->m_num_port_; i++) {
                        if (i != 1) {
                            of << std::setw(15) << " ";
                        }
                        for (j = 1; j <= this->m_num_port_; j++) {
                            std::complex<double> tmp = (value.second)[out_param].at(i - 1, j - 1);
                            of << std::setw(15) << std::real(tmp) << " " << std::setw(15) << std::imag(tmp) << "\t"
                                << " ";
                            if (j % step == 0)
                                of << std::endl;
                        }
                        if ((j - 1) % step != 0)
                            of << std::endl;
                    }
                }

                break;
            }
            }
        }

        of.flags(default_fmt); //恢复默认格式化输出
        of.close();
        return true;
    }

    //延迟计算
    // here out_param
    bool PortMatConv::Apply(const double freq, const PortMatType out_param)
    {
        auto iter = this->m_in_.find(freq);

        // if freq not in origin frequency or current frequency(which may be interpolated) ,then don't calculate.
        if (iter == this->m_in_.end() ||
            std::find(this->m_p_emoption_->frequency.begin(), this->m_p_emoption_->frequency.end(), freq) ==
            this->m_p_emoption_->frequency.end()) {
            return false;
        }

        for (int type_flag = 0x1; type_flag != 0x8; type_flag = type_flag << 1) { // mask
            PortMatType type_tmp = static_cast<PortMatType>(out_param & type_flag);
            if (type_tmp != UNKOWN_PORT_TYPE) {
                arma::cx_mat &out = iter->second[type_tmp];
                if (out.n_elem == 0) {
                    MatConvInternal(this->m_in_param_, type_tmp, this->m_in_.find(freq)->second[m_in_param_], out);
                }
            }
        }

        return true;
    }

    //返回特定频率的单一参数矩阵
    arma::cx_mat &PortMatConv::at(const double freq, const PortMatType out_param)
    {

        if (out_param != S && out_param != Y && out_param != Z)
            throw std::runtime_error("illegal PortMatType param");
        Apply(freq, out_param);
        return this->m_in_.find(freq)->second[out_param];
    }

    arma::cx_mat &PortMatConv::operator()(const double freq, const PortMatType out_param) { return this->at(freq, out_param); }

    //返回特定频率的所有参数矩阵
    bool PortMatConv::GetPortMat(const double freq, arma::cx_mat &s, arma::cx_mat &y, arma::cx_mat &z)
    {
        auto iter = this->m_in_.find(freq);
        if (iter == this->m_in_.end())
            return false;

        PortMatSet &port_set = iter->second;

        //延迟计算
        if (port_set[S].n_elem == 0) {
            MatConvInternal(this->m_in_param_, S, port_set[m_in_param_], port_set[S]);
        }
        if (port_set[Y].n_elem == 0) {
            MatConvInternal(this->m_in_param_, Y, port_set[m_in_param_], port_set[Y]);
        }
        if (port_set[Z].n_elem == 0) {
            MatConvInternal(this->m_in_param_, Z, port_set[m_in_param_], port_set[Z]);
        }

        s = port_set[S];
        y = port_set[Y];
        z = port_set[Z];
        return true;
    }

    //内部调用的方法
    void PortMatConv::MatConvInternal(const PortMatType in_param, const PortMatType out_param, const arma::cx_mat &in,
        arma::cx_mat &out)
    {
        arma::cx_mat e(this->m_num_port_, this->m_num_port_);
        e.eye();
        switch (in_param) {
        case S: {
            if (out_param == PortMatType::Y) { // S => Y
                arma::cx_mat s2y(this->m_num_port_, this->m_num_port_);
                s2y = m_inv_G_in_ * inv(in * m_Z_in_ + m_conj_Z_in_) * (e - in) * m_G_in_;
                out = std::move(s2y);
            }
            else if (out_param == PortMatType::Z) {
                arma::cx_mat s2z(this->m_num_port_, this->m_num_port_);
                s2z = m_inv_G_in_ * inv(e - in) * (in * m_Z_in_ + m_conj_Z_in_) * m_G_in_;
                out = std::move(s2z);
            }
            else {
                ;
            }
            break;
        }
        case Y: {
            if (out_param == PortMatType::S) {
                arma::cx_mat y2s(this->m_num_port_, this->m_num_port_);
                y2s = m_G_in_ * (e - m_conj_Z_in_ * in) * inv(e + m_Z_in_ * in) * m_inv_G_in_;
                out = std::move(y2s);
            }
            else if (out_param == PortMatType::Z) {
                arma::cx_mat y2z(this->m_num_port_, this->m_num_port_);
                y2z = inv(in);
                out = std::move(y2z);
            }
            else {
                ;
            }
            break;
        }
        case Z: {
            if (out_param == PortMatType::S) {
                arma::cx_mat z2s(this->m_num_port_, this->m_num_port_);
                z2s = m_G_in_ * (in - m_conj_Z_in_) * inv(in + m_Z_in_) * m_inv_G_in_;
                out = std::move(z2s);
            }
            else if (out_param == PortMatType::Y) {
                arma::cx_mat z2y(this->m_num_port_, this->m_num_port_);
                z2y = inv(in);
                out = std::move(z2y);
            }
            else {
                ;
            }
            break;
        }
        }
    }
}
