#include<iostream>
#include<vector>
#include<armadillo>
#include<set>
#include<map>
#include "Constant.h"
#include "Element.h"

using namespace std;
using namespace MoM;
using namespace arma;

void TestSquareCapacity();
bool CutSquare(const double error, Element test_sqr, std::vector<Element> &differential_square_vector);


int main()
{
    TestSquareCapacity();

    return 0;
}

struct SqrtComp {
    bool operator()(const Element& lhs, const Element& rhs)const {
        auto a = lhs.GetCenter();
        auto b = rhs.GetCenter();
        if (a[0] < b[0])
            return true;
        else if (a[0] == b[0] && a[1] < b[1])
            return true;
        else
            return false;
    }
};

void TestSquareCapacity()
{
    double length = EMOption::convertLengthUnit(1.0, MM, MM);
    std::vector<vec3> vec_square;
    vec_square.push_back({ -length, -length, 0 });
    vec_square.push_back({ length, -length, 0 });
    vec_square.push_back({ length, length, 0 });
    vec_square.push_back({ -length, length, 0 });

    //获得离散数据
    Element square(0, vec_square, NULL);
    std::vector<Element> differential_square_vector;
    CutSquare(length / 8, square, differential_square_vector);
    if (differential_square_vector.size() == 0) {
        differential_square_vector.push_back(square);
    }

    std::set<double, std::greater<>> s;
    for (auto sqr : differential_square_vector) {
        auto center = sqr.GetCenter();
        s.insert(center[1]);
    }

    map<double, vector<Element>, std::greater<>> m;
    for (auto y : s) {
        m.insert(make_pair(y, vector<Element>()));
    }
    for (auto sqr : differential_square_vector) {
        auto center = sqr.GetCenter();
        auto it = m.find(center[1]);
        (*it).second.push_back(sqr);
    }

    for (auto& it : m) {
        sort(it.second.begin(), it.second.end(), SqrtComp());
    }

    vector<Element> sorted_vector;
    for (auto it : m) {
        cout << it.first << endl;
        for (auto t : it.second)
            sorted_vector.push_back(t);
    }


    cout << "cut size:" << sorted_vector.size() << endl;

    //计算[Lmn]矩阵
    auto &basis_list = sorted_vector;
    int N = basis_list.size();
    int num = sqrt(N);
    double a = length;
    double b = a / sqrt(N);
    mat L(N, N);
    mat G(N, 1, fill::ones), X(N, 1);


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                L(i, j) = 2 * b*log(1 + sqrt(2.0)) / (PI*epsilon0);
            }
            else {
                auto center_f = basis_list[i].GetCenter();
                auto center_s = basis_list[j].GetCenter();
                L(i, j) = b*b / (PI*epsilon0*norm(center_f - center_s, 2));
            }
        }
    }
    L.save("L_mat.txt", raw_ascii);

    solve(X, L, G);

    double sum = 0;
    for (int i = 0; i < X.size(); i++) {
        sum += X[i] * 2 * b * 2 * b;
    }
    cout << sum / (2 * a) *1e12 << "(uuF)" << endl;
    X.reshape(num, num);
    X.save("alpha.mat", raw_ascii);

}


bool CutSquare(const double error, Element test_sqr, std::vector<Element> &differential_square_vector) {
    vector<vec3> vertexes = test_sqr.GetVertex();
    if (norm(vertexes[0] - vertexes[1], 2) <= error) {
        return true;
    }
    auto mid1 = (vertexes[0] + vertexes[1]) / 2.0;
    auto mid2 = (vertexes[1] + vertexes[2]) / 2.0;
    auto mid3 = (vertexes[2] + vertexes[3]) / 2.0;
    auto mid4 = (vertexes[0] + vertexes[3]) / 2.0;
    auto center = (mid1 + mid3) / 2;
    Element square1(0, { vertexes[0], mid1, center, mid4 }, NULL);
    Element square2(0, { vertexes[1], mid2, center, mid1 }, NULL);
    Element square3(0, { vertexes[2], mid3, center, mid2 }, NULL);
    Element square4(0, { vertexes[3], mid4, center, mid3 }, NULL);

    bool flag = false;
    flag = flag || CutSquare(error, square1, differential_square_vector);
    flag = flag || CutSquare(error, square2, differential_square_vector);
    flag = flag || CutSquare(error, square3, differential_square_vector);
    flag = flag || CutSquare(error, square4, differential_square_vector);
    if (flag == true) {
        differential_square_vector.push_back(square1);
        differential_square_vector.push_back(square2);
        differential_square_vector.push_back(square3);
        differential_square_vector.push_back(square4);
    }
    return false;
}
