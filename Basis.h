#ifndef _BASIS_H_
#define _BASIS_H_

namespace mom
{
    class Basis
    {
    protected:
        int index;
        int element_index[2];

    public:
        Basis();
        Basis(const int element_index_a, const int element_index_b);
        ~Basis();

        // Get the index by subscript
        int operator[](const int subscript) const;
    };
}
#endif // _BASIS_H_
