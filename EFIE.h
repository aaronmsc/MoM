#ifndef _EFIE_H_
#define _EFIE_H_
#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <armadillo>
#include "Typedef.h"
#include "MoM.h"
#include "Port.h"
#include "RWG.h"
#include "Element.h"
#include "Basis.h"
#include "EMOption.h"
#include "PortMatConv.h"


namespace mom
{

    class EFIE
    {
    private:
        EMOption *m_p_emoption_;
        std::vector<Element> *element_list;
        std::vector<Basis> *basis_list;
        std::map<real, std::vector<arma::cx_vec>> m_map_I_list_; // map of freq to I_list
        std::map<real, arma::cx_mat> m_map_Y_list_;
        std::vector<arma::cx_vec> m_F_list_;
        std::vector<RWG> rwg_list;
        // F_list
        int m_port_num_;

        void CreateVFeedMatrixRWG(const std::vector<RWG> &rwg_list, const struct PortGroup &feed_port, arma::cx_vec &F);
        void CreateGapFeedMatrixRWG(const std::vector<RWG> &rwg_list, const PortGroup &feed_port, arma::cx_vec &F);
        void CreateJOfVFeedPortRWG(const std::vector<RWG> &rwg_list, const std::vector<PortGroup> &port_group_array,
            const arma::cx_vec &I, arma::cx_vec &J);
        void CreateJOfGapFeedPortRWG(const std::vector<RWG> &rwg_list, const std::vector<PortGroup> &port_group_array,
            const arma::cx_vec &I, arma::cx_vec &J);
        real CalculatePortGap(const PortGroup &port_group, const arma::vec3 &direction_V);

    public:
        PortMatConv port_mat_conv;
        EFIE() {}
        inline void SetEMOption(EMOption *sp) { this->m_p_emoption_ = sp; }
        void Clear()
        {
            m_map_I_list_.clear();
            m_map_Y_list_.clear();
            m_F_list_.clear();
            rwg_list.clear();
            port_mat_conv.Clear();
        }



        inline void SetElemAndBasisList(std::vector<Element> *element_list, std::vector<Basis> *basis_list)
        {
            this->element_list = element_list;
            this->basis_list = basis_list;
        }

        void Run();
        int GetPortNum() const { return this->m_port_num_; }
        void UnifyGeoUnit();
    };
}
#endif
