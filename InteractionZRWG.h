#ifndef _INTERACTIONZ_RWG_HPP_
#define _INTERACTIONZ_RWG_HPP_
#include <cmath>
#include <unordered_map>
#include "EMOption.h"
#include "Constant.h"
#include "Typedef.h"
#include "RWG.h"

namespace mom
{





    struct CoordinateCompare {
        bool operator()(double a, double b) const
        {
            if (std::abs(a - b) < 1e-12) {
                return false;
            }
            return a < b;
        }
    };
    namespace internal{
        extern double g_z_GND;
        struct LocalConstant {
            real lambda;
            real k0;
            FComplex ik;
            real frequency;
            LocalConstant(real freq)
            {
                this->lambda = c_wave_speed / freq;
                this->k0 = 2 * PI / lambda;
                this->ik = FComplex(0, 1) * k0;
                this->frequency = freq;
            }
        };


        enum Axis { X_AXIS, Y_AXIS, Z_AXIS };
        inline real ReflectDistance(const arma::vec3 &field, const arma::vec3 &src, internal::Axis axis)
        {
            arma::vec3 reflect_point(field);
            switch (axis) {
            case X_AXIS:
                reflect_point[0] = -reflect_point[0];
                break;
            case Y_AXIS:
                reflect_point[1] = -reflect_point[1];
                break;
            case Z_AXIS:
                reflect_point[2] = 2 * internal::g_z_GND - reflect_point[2];
                break;
            default:
                break;
            }
            return arma::norm(reflect_point - src, 2);
        }
        inline void GAInterpolation(const real distance, const internal::LocalConstant &local_constant, arma::cx_vec &GA)
        {
#ifdef _DEBUG_ENGINE_ALONE_
            if (distance == 0)
                Msg::Warning("distance == 0!");
#endif
            GA[0] = GA[1] = GA[2] = std::exp(local_constant.ik * distance) / distance;
        }

        inline void GReflectInterpolation(const real distance, const internal::LocalConstant &local_constant, Axis axis,
            arma::cx_vec &G_reflect)
        {
            switch (axis) {
            case X_AXIS:
                G_reflect[0] = std::exp(local_constant.ik * distance) / distance;
                G_reflect[1] = G_reflect[2] = -G_reflect[0];
                break;
            case Y_AXIS:
                G_reflect[1] = std::exp(local_constant.ik * distance) / distance;
                G_reflect[0] = G_reflect[2] = -G_reflect[1];
                break;
            case Z_AXIS:
                G_reflect[2] = std::exp(local_constant.ik * distance) / distance;
                G_reflect[1] = G_reflect[0] = -G_reflect[2];
                break;
            default:
                break;
            }
        }
    }

    // struct hash_pair {
    // 	size_t operator()(const pair<int, int> &Pair) const { return ((Pair.first) << 15) ^ (Pair.second); }
    // };

    class InteractionZRWG
    {
    protected:
        SpaceType m_space_;
        // std::unordered_map<std::pair<int, int>,,hash_pair >

        inline FComplex IntegralZMomRWG(const RWG &field, const RWG &src,
            const internal::LocalConstant &local_constant);
        inline FComplex IntegralZMomRWGDiffKernel(const RWG &field, const RWG &src, const int idx_figure_f,
            const int idx_figure_s, const internal::LocalConstant &local_constant);
        inline FComplex IntegralOverROfTri(const std::vector<arma::vec3> &vertexes, const arma::vec3 &gauss_point);

        inline FComplex IntegralSingularFigureRWGFreespace(const RWG &field, const RWG &src, int idx_figure_f,
            int idx_figure_s, const internal::LocalConstant &local_constant);
        inline FComplex IntegralNonsingularFigureRWGFreespace(const RWG &field, const RWG &src, int idx_figure_f,
            int idx_figure_s,
            const internal::LocalConstant &local_constant);
        inline FComplex IntegralSingularFigureRWGGNDLayer(const RWG &field, const RWG &src, int idx_figure_f,
            int idx_figure_s, const internal::LocalConstant &local_constant);

        inline FComplex IntegralNonsingularFigureRWGGNDLayer(const RWG &field, const RWG &src, int idx_figure_f,
            int idx_figure_s,
            const internal::LocalConstant &local_constant);

    public:
        InteractionZRWG() {}
        InteractionZRWG(SpaceType space)
            :m_space_(space) {}
        void Init(SpaceType space)
        {
            this->m_space_ = space;
        }

        FComplex Calculate(const RWG &field, const RWG &src, const real freq);
    };


}
#endif
