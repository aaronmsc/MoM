#include<iostream>
#include<vector>
#include<armadillo>
#include<set>
#include<map>
#include "Constant.h"
#include "Element.h"
#include "Basis.h"
#include "EFIE.h"
#include "DataInput.h"

using namespace std;
using namespace mom;
//using namespace arma;

void TestSquareCapacity();
bool CutSquare(const double error, Element test_sqr, std::vector<Element> &differential_square_vector);
void MoMTest();

int main()
{
    //TestSquareCapacity();
    MoMTest();
    system("pause");
    return 0;
}



void MoMTest()
{



    std::fstream files("file_list.txt");
    if (!files.is_open()) {
        MSG("file_list.txt don't exist. please create files");
        system("pause");
        exit(1);
    }
    int output_num, input_num, input_select, output_select;
    std::string output_prefix[128], input_prefix[128], in_select_prefix, out_select_prefix, out_file;
    files >> input_num;
    for (int i = 0; i < input_num; i++) {
        files >> input_prefix[i];
    }
    files >> output_num;
    for (int i = 0; i < output_num; i++) {
        files >> output_prefix[i];
    }
    files.close();

    MSG("Please choose input prefix:(Choose 0 for default input_file=1,output_file=./matlab.txt)");
    for (int i = 0; i < input_num; i++) {
        MSG("[%d]%s", i + 1, input_prefix[i].c_str());
    }
    std::cin >> input_select;
    if (input_select > 0) {
        MSG("Please choose output prefix:");
        for (int i = 0; i < output_num; i++) {
            MSG("[%d]%s", i + 1, output_prefix[i].c_str());
        }
        std::cin >> output_select;
        if (input_select > input_num || output_select > input_num || output_select <= 0) {
            MSG("Invalid input , please check your input");
            system("pause");
            exit(1);
        }
        MSG("Please input output file name:");
        std::cin >> out_file;
        input_select--;
        output_select--;
        out_file = output_prefix[output_select] + out_file;
    }
    else if (input_select == 0) {
        input_select = output_select = 0;
        out_file = output_prefix[output_select] + "matlab.txt";
    }
    else {
        MSG("Invalid input , please check your input");
        system("pause");
        exit(1);
    }
    in_select_prefix = input_prefix[input_select];

    std::fstream in(in_select_prefix + "rwg_option.txt");
    if (!in.is_open()) {
        MSG("rwg_option file doesn't exist");
        system("pause");
        exit(1);
    }

    mom::EFIE *p_EFIE;
    p_EFIE = new mom::EFIE();

    SpaceType space_type;
    KernelMethod kernel_method;
    AccurLevel accur_level;
    IntegralApproach integral_approach;
    Unit unit;
    FeedType feed_type;
    std::vector<Basis> basis_list;
    std::vector<Element> elem_list;
    std::vector<Port *> port_list;
    std::vector<std::string> feed_ports;
    std::vector<double> freq_def;
    double freq_tmp;
    std::vector<std::string> files_name;
    std::string name_tmp;
    std::string line;
    std::stringstream ss;
    std::string port_pair;

    // input option
    SkipComment(in, ss);
    space_type = static_cast<SpaceType>(GetInt(ss));
    SkipComment(in, ss);
    kernel_method = static_cast<KernelMethod>(GetInt(ss));
    SkipComment(in, ss);
    integral_approach = static_cast<IntegralApproach>(GetInt(ss));
    SkipComment(in, ss);
    unit = static_cast<Unit>(GetInt(ss));
    SkipComment(in, ss);
    accur_level = static_cast<AccurLevel>(GetInt(ss));
    SkipComment(in, ss);
    feed_type = static_cast<FeedType>(GetInt(ss));
    SkipComment(in, ss);
    while (std::getline(ss, port_pair, ' ')) {
        if (!port_pair.empty())
            feed_ports.push_back(port_pair);
    }


    SkipComment(in, ss);
    for (size_t i = 0; i < 3; i++) {

        ss >> freq_tmp;
        freq_def.push_back(freq_tmp);
    }

    // files
    while (in >> name_tmp, !in.eof()) {
        SkipComment(in, ss);
        files_name.push_back(name_tmp);
    }

    // Options Instantiate
    // Linear distribution
    EMOption* p_option = new EMOption(freq_def, feed_ports, { FComplex(50, 0) }, unit, space_type, kernel_method, feed_type, LINEAR_STEP, TRIANGLE, { "Y_Parameter", "S_Parameter" }, out_file);
    // Logarithm distribution
    // EMOption* p_option(new EMOption(freq_def, { "net1.P1" }, { FComplex(50, 0) }, unit, space_type,
    // 												IC_STANDARD, kernel_method, VOLTAGE, LOGARITHM, TRIANGLE,
    // 												RECTANGULAR, { "Y_Parameter" }, out_file));
    p_EFIE->SetEMOption(p_option);


    DataInput data_input;

    MSG("Start ReadData");
    data_input.ReadData(in_select_prefix, files_name, basis_list, elem_list, port_list);
    MSG("End ReadData");

    // run
    p_EFIE->SetElemAndBasisList(&elem_list, &basis_list);

    time_t begin, end;

    begin = time(NULL);
    p_EFIE->Run();
    end = time(NULL);
    cout << "time cost: { " << (end - begin) / 60 / 60 << ":" << ((end - begin) / 60) % 60 << ":" << (end - begin) % 60
        << "}" << endl;
    int port_num = p_EFIE->GetPortNum();
    // MatlabFormat(port_mat.GenerateInfor(S+Y),S+Y,port_num);
    // PythonFormat(p_EFIE->port_mat_conv.GenerateInfor(S), S, port_num, out_file);
    // uniq_pOutput->Save(out_file);
    if (p_EFIE) {
        delete p_EFIE;
        p_EFIE = 0;
    }
    if (p_option) {
        delete p_option;
        p_option = 0;
    }

    for (size_t i = 0; i < port_list.size(); i++) {
        delete port_list[i];
    }
    in.close();
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
    std::vector<arma::vec3> vec_square;
    vec_square.push_back({ -length, -length, 0 });
    vec_square.push_back({ length, -length, 0 });
    vec_square.push_back({ length, length, 0 });
    vec_square.push_back({ -length, length, 0 });

    //获得离散数据
    Element square(0, vec_square, NULL);
    std::vector<Element> differential_square_vector;
    CutSquare(length / 2, square, differential_square_vector);
    if (differential_square_vector.size() == 0) {
        differential_square_vector.push_back(square);
    }

    //对square按在正方形中的位置从左上到右下进行排序
    std::set<double, std::greater<>> s_y;
    std::set<double, std::less<>> s_x;
    for (auto sqr : differential_square_vector) {
        auto center = sqr.GetCenter();
        s_y.insert(center[1]);
        s_x.insert(center[0]);
    }

    map<double, vector<Element>, std::greater<>> m;
    for (auto y : s_y) {
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
    int num = static_cast<int>(std::sqrt(N));
    double a = length;
    double b = a / std::sqrt(N);
    arma::mat L(N, N);
    arma::mat G(N, 1, arma::fill::ones), X(N, 1);


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                L(i, j) = 2 * b*std::log(1 + std::sqrt(2.0)) / (PI*epsilon0);
            }
            else {
                auto center_f = basis_list[i].GetCenter();
                auto center_s = basis_list[j].GetCenter();
                L(i, j) = b*b / (PI*epsilon0*norm(center_f - center_s, 2));
            }
        }
    }
    //L.save("output/L_mat.txt", raw_ascii);

    //求解线性方程组
    solve(X, L, G);

    //计算正方板的电容
    double sum = 0;
    for (size_t i = 0; i < X.size(); i++) {
        sum += X[i] * 2 * b * 2 * b;
    }
    cout << sum / (2 * a) *1e12 << "(uuF)" << endl;



    //存储格点信息
    ofstream file_X("output/figure_info_X.txt");
    if (!file_X.is_open()) {
        cout << "cann't open figure_info.txt" << endl;
        exit(1);
    }
    for (auto value : s_x) {
        file_X << value << endl;
    }
    file_X.close();

    ofstream file_Y("output/figure_info_Y.txt");
    if (!file_Y.is_open()) {
        cout << "cann't open figure_info.txt" << endl;
        exit(1);
    }
    for (auto value : s_y) {
        file_Y << value << endl;
    }
    file_Y.close();

    X.reshape(num, num);
    X.save("output/figure_info_Z.txt", arma::raw_ascii);


}

bool CutSquare(const double error, Element test_sqr, std::vector<Element> &differential_square_vector) {
    vector<arma::vec3> vertexes = test_sqr.GetVertex();
    if (arma::norm(vertexes[0] - vertexes[1], 2) <= error) {
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
