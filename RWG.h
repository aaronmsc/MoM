#ifndef _RWG_H_
#define _RWG_H_
#include <cassert>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>


#include "Basis.h"
#include "EMOption.h"
#include "Element.h"
#include  "Port.h"


namespace mom
{

    namespace internal
    {
        enum PortBasisFlag { ABORT = 0, NORMAL, EDGE, AREA };
    }

    struct RWGGroup;

    class RWG
    {
    private:
        Element *element[2];
        bool is_halfbasis;
        real common_length;
        int index;
        // the index of the vertexes in element which does not belong to the common edge
        int original_point[2];
        int current_direction[2];
        internal::PortBasisFlag port_basis_flag;
        std::vector<arma::vec3> rho_list[2];
        std::vector<arma::vec3> diff_rho_list[2];

        // the original index of the artificial basis

        bool CheckSubscript(int subscript) const;
        std::vector<arma::vec3> CalculateRho(int subscript);
        std::vector<arma::vec3> CalculateDiffRho(int subscript);

    public:
        RWG();
        RWG(const Basis &basis, std::vector<Element> &element_list);
        ~RWG() {}

        inline void SetPortBasisFlag(internal::PortBasisFlag port_basis_flag) { this->port_basis_flag = port_basis_flag; }
        inline internal::PortBasisFlag GetPortBasisFlag() const { return this->port_basis_flag; }


        void CreateCompleteRWG(int original_index);
        void CreateHalfRWG(int original_index);

        bool IsHalfbasis() const;
        Element *GetElement(int subscript) const;
        real GetCommonLength() const;
        int GetOriginalPoint(int subscript) const;
        arma::vec3 GetElementCentre(int subscript) const;

        //
        const std::vector<arma::vec3> &GetRho(int subscript) const
        {
            // assert(CheckSubscript(subscript));
            return this->rho_list[subscript];
        }
        const std::vector<arma::vec3> &GetDiffRho(int subscript) const
        {
            // assert(CheckSubscript(subscript));
            return this->diff_rho_list[subscript];
        }
        arma::vec3 GetGaussianIntegral(int subscript, const std::vector<arma::vec3> &rho) const;
        arma::vec3 GetDiffGaussianIntegral(int subscript, const std::vector<arma::vec3> &diff_rho) const;

        real GetDivergence(int subscript) const;

        inline static internal::PortBasisFlag IsUsefulBasis(const Basis &basis, std::vector<Element> &element_list,
            FeedType feed_type)
        {
            if (feed_type == GAP)
                return internal::NORMAL;
            else if (feed_type == VOLTAGE) {
                if (basis[0] < 0 && element_list[basis[1]].GetPort()) {
                    // MSG("catch Edge RWG [%d,%d]", basis[0], basis[1]);
                    return internal::EDGE;
                }
                else if (basis[1] < 0 && element_list[basis[0]].GetPort()) {
                    // MSG("catch Edge RWG [%d,%d]", basis[0], basis[1]);
                    return internal::EDGE;
                }
                else if (basis[0] >= 0 && basis[1] >= 0) {
                    if (element_list[basis[0]].GetPort() || element_list[basis[1]].GetPort()) {
                        // MSG("catch Area RWG [%d,%d]", basis[0], basis[1]);
                        return internal::AREA;
                    }
                    else
                        return internal::NORMAL;
                }
                else {
                    MSG("abord [ %d,%d ]", basis[0], basis[1]);
                    return internal::ABORT;
                }
            }
            return internal::ABORT;
        }




        real GetCurrentDirection(int subscript) const;

        arma::vec3 GetPointByParam(int subscript, real u0, real u1, real u2) const;

        real GetJacobi(int subscript) const;
        real GetRadius(int subscript) const;
    };
}
#endif // RWGBASIS_H_
