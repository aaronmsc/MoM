#ifndef _PORT_MAT_STD_H
#define _PORT_MAT_STD_H
#include "EMOption.h"
namespace mom
{
    class PortMatSet
    {
    public:
        arma::cx_mat port_mat[3];
        PortMatSet() {}
        PortMatSet(const PortMatSet &rhs)
        {
            port_mat[0] = rhs.port_mat[0];
            port_mat[1] = rhs.port_mat[1];
            port_mat[2] = rhs.port_mat[2];
        }
        PortMatSet(PortMatSet &&prvalue)
        {
            port_mat[0] = std::move(prvalue.port_mat[0]);
            port_mat[1] = std::move(prvalue.port_mat[1]);
            port_mat[2] = std::move(prvalue.port_mat[2]);
        }
        PortMatSet(int dim)
        {
            for (size_t i = 0; i < 3; i++)
                port_mat[i] = arma::cx_mat(dim, dim);
        }
        inline void Resize(int dim)
        {
            for (size_t i = 0; i < 3; i++)
                port_mat[i].resize(dim, dim);
        }
        inline arma::cx_mat &operator[](PortMatType port_type)
        {
            if (port_type & 0x7)
                return this->port_mat[port_type >> 1];
            else
                throw std::runtime_error("Unknow port type found");
        }
        inline arma::cx_mat operator[](PortMatType port_type) const
        {
            if (port_type & 0x7)
                return this->port_mat[port_type >> 1];
            else
                throw std::runtime_error("Unknow port type found");
        }
        inline arma::cx_mat &operator[](int i)
        {
            if (i & 0x7)
                return this->port_mat[i >> 1];
            else
                throw std::runtime_error("Unknow port type found");
        }
        inline arma::cx_mat operator[](int i) const
        {
            if (i & 0x7)
                return this->port_mat[i >> 1];
            else
                throw std::runtime_error("Unknow port type found");
        }
    };
}
#endif // PORT_MAT_STD_H
