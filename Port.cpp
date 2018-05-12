#include"Port.h"

namespace mom {
    Port::Port()
    {
    }
    Port::Port(const std::vector<arma::vec3>& rect, const uint32_t & index, const std::string portname)
    {
        InitPortByRect(rect, index, portname);
    }
    Port::Port(const Port & rhs)
    {
        InitPortByRect(rhs.m_rect_port_, rhs.m_index_, rhs.m_port_name_);
    }
    Port::~Port()
    {
    }
    Port&  Port::operator=(const Port & rhs)
    {
        InitPortByRect(rhs.m_rect_port_, rhs.m_index_, rhs.m_port_name_);
        return *this;
    }
    bool Port::operator==(const Port & rhs) const
    {
        if (this->m_port_name_ == rhs.m_port_name_ && this->m_index_ == rhs.m_index_)
            return true;
        return false;
    }
    void Port::InitPortByRect(const std::vector<arma::vec3>& rectangular, const int index, const std::string & portname)
    {
        for (auto point : rectangular) {
            this->m_rect_port_.push_back(point);
        }
        this->m_port_name_ = portname;
        this->m_index_ = index;
    }
    void Port::SetPortName(const std::string & portname)
    {
        this->m_port_name_ = portname;
    }
    std::string Port::GetPortName() const
    {
        return this->m_port_name_;
    }
    std::vector<arma::vec3> Port::GetVecPoints()const
    {
        return this->m_rect_port_;
    }
    int Port::GetIndex() const
    {
        return this->m_index_;
    }
    void Port::SetIndex(const int index)
    {
        this->m_index_ = index;
    }
    std::vector<arma::vec3> Port::GetPortPnts()
    {
        return std::vector<arma::vec3>();
    }
    void Port::SetRectPort(const std::vector<arma::vec3>& rectangular)
    {
        std::copy(rectangular.begin(), rectangular.end(), this->m_rect_port_.begin());
    }
}
