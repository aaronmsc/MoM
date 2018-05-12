#ifndef _PORT_H_               
#define _PORT_H_
#include <iostream>
#include <armadillo>
#include <algorithm>
#include <vector>
#include <set>
#include <list>
#include "EMOption.h"

namespace mom {

    class Port
    {

    private:
        //arma::vec3 mBas, mVec;
        int m_index_;
        std::vector<arma::vec3> m_rect_port_;
        std::string m_port_name_;
    private:
        //void GetPoints(std::vector<arma::vec3> &VecPnt);
        //void CalcRectPoints(const arma::vec3 &P1i, const arma::vec3 &Pi, arma::vec3 &xPnt, arma::vec3 &yPnt);
    public:


        Port();
        Port(const std::vector<arma::vec3>& rect, const uint32_t& index, const std::string portname = "");
        //Port(const arma::vec3& bas, const arma::vec3& vec, const uint32_t& index, const std::string portname, const bool polarity);
        Port(const Port &pot);

        virtual ~Port();

        Port& operator=(const Port &rhs);



        bool	operator==(const Port &pt)const;

        //void InitPort_ByVec(const arma::vec3& bas, const arma::vec3& vec, const int index, const std::string& portname, const bool polarity);
        void InitPortByRect(const std::vector<arma::vec3>& vPnt, const int index, const std::string& portname);
        void SetPortName(const std::string& portname);
        std::string GetPortName()const;

        std::vector<arma::vec3> GetVecPoints()const;

        int GetIndex()const;


        void SetIndex(const int index);

        //arma::vec3& GetBasPnt();

        //const arma::vec3 GetBasPnt()const;

        //arma::vec3 GetDirection()const;

        std::vector<arma::vec3> GetPortPnts();
        void  SetRectPort(const  std::vector<arma::vec3>& rectangular);

    };
    typedef struct PortBasis {
        bool voltage_direction;
        const Port *port;
        int basis_ix;
        PortBasis(int basis_ix, const Port *port, bool voltage_direction)
            : basis_ix(basis_ix), port(port), voltage_direction(voltage_direction)
        {
        }
    } PortBasis;
    typedef struct PortGroup {
        std::set<const Port *> port_set_pos;
        std::set<const Port *> port_set_neg;
        std::list<PortBasis> basis_list_pos;
        std::list<PortBasis> basis_list_neg;
        std::vector<std::string> string_array_pos;
        std::vector<std::string> string_array_neg;
        FeedType feed_type;
        int idx;
        enum PortDirection { UNKNOWN_DIRECTION = 0, POSITIVE = 1, NEGATIVE = 2, BOTH = 3 };
        PortGroup(const std::string &s)
        {
            size_t semicolon_pos = s.find_last_of(";");
            size_t colon_pos = s.find_last_of(":");

            bool is_differential = semicolon_pos != std::string::npos;

            std::vector<std::string> v1, v2;
            std::string s1 = s.substr(0, is_differential ? semicolon_pos : colon_pos);
            std::string s2 = s.substr(semicolon_pos + 1, colon_pos - semicolon_pos - 1);
            std::string s3 = s.substr(colon_pos);
            std::string si;

            this->feed_type = VOLTAGE;
            std::stringstream ss(s1);

            while (getline(ss, si, ',')) {
                this->string_array_pos.push_back(si + s3);
            }
            if (is_differential) {
                std::stringstream ss(s2);
                while (getline(ss, si, ','))
                    this->string_array_neg.push_back(si + s3);
            }
        }
        PortDirection Find(const std::string &s, PortDirection direction = BOTH)
        {
            int direction_flag = static_cast<int>(direction);
            if ((direction_flag & 0x1) &&
                find(string_array_pos.begin(), string_array_pos.end(), s) != string_array_pos.end())
                return POSITIVE;
            if ((direction_flag & 0x2) &&
                find(string_array_neg.begin(), string_array_neg.end(), s) != string_array_neg.end())
                return NEGATIVE;
            return UNKNOWN_DIRECTION;
        }
    } PortGroup;
}
#endif //_PORT_H_
