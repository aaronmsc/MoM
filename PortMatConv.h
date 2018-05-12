#ifndef _PORT_MAT_CONV_H
#define _PORT_MAT_CONV_H


#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>
#include "Constant.h"
#include "EMOption.h"
#include "PortMatStd.h"

namespace mom
{

    typedef struct OutputParamInfor OutputParamInfor;

    class PortMatConv
    {
    public:
        PortMatConv() {}
        void Clear() { m_in_.clear(); }

        void Initialize(const PortMatType in_param, const int num_port, const std::map<double, arma::cx_mat> &param_map,
            const arma::cx_mat &Z_in, const EMOption *pOption);
        bool Save(const std::string file_name, const PortMatType type = S, const DataFormat format = MA,
            const FrequencyUnit freq_unit = GHZ, const int n = 50);
        bool SaveImplement(const std::string file_name, const PortMatType out_param, const DataFormat format,
            const FrequencyUnit freq_unit, const int n);
        std::vector<OutputParamInfor> GenerateInfor(PortMatType type);
        std::vector<OutputParamInfor> GenerateInfor(int type);


        arma::cx_mat &at(const double freq, const PortMatType out_param);
        arma::cx_mat &operator()(const double freq, const PortMatType out_param);
        bool GetPortMat(const double freq, arma::cx_mat &s, arma::cx_mat &y, arma::cx_mat &z);

    private:
        void MatConvInternal(const PortMatType in_param, const PortMatType out_param, const arma::cx_mat &in, arma::cx_mat &out);
        bool Apply(const double freq, const PortMatType out_param);

        const EMOption *m_p_emoption_;
        std::map<double, PortMatSet> m_in_;
        int m_num_port_;		 //端口数
        PortMatType m_in_param_; //输入类型
        arma::cx_mat m_Z_in_;			 //输入阻抗矩阵 Z矩阵
        arma::cx_mat m_G_in_;			 //有阻抗矩阵转化得到的 G矩阵
        arma::cx_mat m_conj_Z_in_;
        arma::cx_mat m_inv_G_in_;
    };

    typedef struct OutputParamInfor {
        double freq;
        std::string param_type;
        std::string param_port;
        double magnitude;
        double magnitude_dB;
        double phase;
        double real_part;
        double image_part;
        OutputParamInfor(const OutputParamInfor &rhs)
        {
            this->freq = rhs.freq;
            this->param_type = rhs.param_type;
            this->param_port = rhs.param_port;
            this->magnitude = rhs.magnitude;
            this->magnitude_dB = rhs.magnitude_dB;
            this->phase = rhs.phase;
            this->real_part = rhs.real_part;
            this->image_part = rhs.image_part;
        }
        OutputParamInfor() {}
        int operator==(const OutputParamInfor &rhs) const
        {

            if (this->freq == rhs.freq && this->param_type == rhs.param_type && this->magnitude == rhs.magnitude &&
                this->magnitude_dB == rhs.magnitude_dB && this->phase == rhs.phase && this->real_part == rhs.real_part &&
                this->image_part == rhs.image_part &&
                this->param_port == rhs.param_port) {
                return 1;
            }
            return 0;
        }
    } OutputParamInfor;
}
#endif
