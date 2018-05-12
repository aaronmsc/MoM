#include "Basis.h"

namespace mom {

    Basis::Basis() {}

    Basis::Basis(const int element_index_a, const int element_index_b)
    {
        element_index[0] = element_index_a;
        element_index[1] = element_index_b;
    }

    Basis::~Basis() {}

    int Basis::operator[](const int subscript) const
    {
        if (subscript != 0 && subscript != 1)
            return -1;

        return this->element_index[subscript];
    }
}
