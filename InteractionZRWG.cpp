#include "InteractionZRWG.h"
namespace mom {
    FComplex InteractionZRWG::Calculate(const RWG &field, const RWG &src, const real freq)
    {
        internal::LocalConstant local_constant(freq);
        return this->IntegralZMomRWG(field, src, local_constant) * local_constant.k0;
    }

    //求Z矩阵中的元素
    //对场点、源点的basis function的两边的图形的高斯点进行积分
    FComplex InteractionZRWG::IntegralZMomRWG(const RWG &field, const RWG &src,
        const internal::LocalConstant &local_constant)
    {
        FComplex integral_result = 0.0;
        int i, j;
        bool is_single_basis_f = field.IsHalfbasis();
        bool is_single_basis_s = src.IsHalfbasis();
        if (!is_single_basis_f && !is_single_basis_s) { //源和目标都有两个element
            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++) {
                    integral_result += IntegralZMomRWGDiffKernel(field, src, j, i, local_constant);
                }
        }
        else if (is_single_basis_f && is_single_basis_s) { //目标和源都只有一个element
            integral_result += IntegralZMomRWGDiffKernel(field, src, 0, 0, local_constant);
        }
        else if (is_single_basis_f) { //目标只有一个element
            j = 0;
            for (i = 0; i < 2; i++) {
                integral_result += IntegralZMomRWGDiffKernel(field, src, j, i, local_constant);
            }
        }
        else { //源只有一个element
            i = 0;
            for (j = 0; j < 2; j++) {
                integral_result += IntegralZMomRWGDiffKernel(field, src, j, i, local_constant);
            }
        }
        return integral_result;
    }

    //对于freespace和分层情况调用不同的方法
    FComplex InteractionZRWG::IntegralZMomRWGDiffKernel(const RWG &field, const RWG &src,
        const int idx_figure_f, const int idx_figure_s,
        const internal::LocalConstant &local_constant)
    {
        FComplex integral_result;

        int idx_elem_s = (src.GetElement(idx_figure_s))->GetIndex();
        int idx_elem_f = (field.GetElement(idx_figure_f))->GetIndex();

        if (this->m_space_ == FREE_SPACE) { //处理自由空间问题
            if (idx_elem_s == idx_elem_f) { //奇异情况，即两个basis function的同一图形时
                integral_result =
                    IntegralSingularFigureRWGFreespace(field, src, idx_figure_f, idx_figure_s, local_constant);
            }

            else { //非奇异情况

                integral_result =
                    IntegralNonsingularFigureRWGFreespace(field, src, idx_figure_f, idx_figure_s, local_constant);
            }
        }
        else if (this->m_space_ == GND_FREESPACE) {
            if (idx_elem_s == idx_elem_f) { //奇异情况，即两个basis function的同一图形时
                integral_result = IntegralSingularFigureRWGGNDLayer(field, src, idx_figure_f, idx_figure_s, local_constant);
                //std::cout << "Singular GND:" << integral_result << std::endl;
            }

            else { //非奇异情况
                integral_result =
                    IntegralNonsingularFigureRWGGNDLayer(field, src, idx_figure_f, idx_figure_s, local_constant);
                //std::cout << "NonSingular GND:" << integral_result << std::endl;
            }
        }
        return integral_result;
    }

    //用于三角形奇异情况的处理，即计算1/R当R为0的情况
    FComplex InteractionZRWG::IntegralOverROfTri(const std::vector<arma::vec3> &vertexes,
        const arma::vec3 &gauss_point)
    {
        FComplex integral_result = 0;
        auto singular_integral = [&integral_result](arma::vec3 p1, arma::vec3 p2, arma::vec3 p3) {
            arma::vec3 v_1_3, v_3_2;
            real jacob = 0, a = 0, b = 0, c = 0;
            v_1_3 = p1 - p3;
            v_3_2 = p3 - p2;
            //?
            jacob = norm(arma::cross(v_1_3, v_3_2), 2);
            a = dot(v_1_3, v_1_3);
            b = 2 * arma::dot(v_3_2, v_1_3);
            c = arma::dot(v_3_2, v_3_2);

            integral_result += jacob * (log(2 * sqrt(c * (a + b + c)) + 2 * c + b) - log(2 * sqrt(c * a) + b)) / sqrt(c);
        };
        singular_integral(gauss_point, vertexes[0], vertexes[1]);
        singular_integral(gauss_point, vertexes[1], vertexes[2]);
        singular_integral(gauss_point, vertexes[2], vertexes[0]);
        return integral_result;
    }

    //自由空间奇异情况处理
    //奇异情况将Green函数分成两部分进行计算：
    // Green = (exp(ik*R)-1)/R + 1/R
    //然后分别对两部分进行积分
    FComplex
        InteractionZRWG::IntegralSingularFigureRWGFreespace(const RWG &field, const RWG &src,
            int idx_figure_f, int idx_figure_s,
            const internal::LocalConstant &local_constant)
    {
        int i, j;
        FComplex integral_result, integral_tmp1, one_Rth, green, green0, z0_for_loss;
        real div_basis_fn_s = 0, div_basis_fn_f = 0;
        real jacob_f = 0, distance = 0, divj_dot_divj = 0, j_dot_j = 0;
        std::vector<arma::vec3> list_basis_fn_s, list_basis_fn_f;
        std::vector<arma::vec3> list_gauss_point_s, list_gauss_point_f;
        std::vector<real> list_gauss_weight_s, list_gauss_weight_f;
        int num_gauss_nodes_s = (src.GetElement(idx_figure_s))->GetGaussPoint().size();
        int num_gauss_nodes_f = (field.GetElement(idx_figure_f))->GetGaussPoint().size();
        real omega = 2 * pi * local_constant.frequency;

        //提取rwg内部需要的信息：basis fn List、div basis fn List、图形的高斯点、权重 List
        list_basis_fn_s = src.GetRho(idx_figure_s);
        list_basis_fn_f = field.GetRho(idx_figure_f);

        div_basis_fn_s = src.GetDivergence(idx_figure_s) / local_constant.k0;
        div_basis_fn_f = field.GetDivergence(idx_figure_f) / local_constant.k0;

        list_gauss_point_s = (src.GetElement(idx_figure_s))->GetGaussPoint();
        list_gauss_point_f = (field.GetElement(idx_figure_f))->GetGaussPoint();

        list_gauss_weight_s = (src.GetElement(idx_figure_s))->GetGaussWeight();
        list_gauss_weight_f = (field.GetElement(idx_figure_f))->GetGaussWeight();

        jacob_f = field.GetJacobi(idx_figure_f);
        divj_dot_divj = div_basis_fn_s * div_basis_fn_f;
        for (i = 0; i < num_gauss_nodes_f; i++) {
            j_dot_j = arma::dot(list_basis_fn_s[i], list_basis_fn_f[i]);
            for (j = 0; j < num_gauss_nodes_s; j++) {
                distance = arma::norm(list_gauss_point_f[i] - list_gauss_point_s[j], 2);

                if (distance < 1e-12) { //如果两个高斯点近似相同  //?精度要求
                    integral_tmp1 = (dot(list_basis_fn_f[i], list_basis_fn_s[j]) - divj_dot_divj) * local_constant.ik;
                }

                else { //两个高斯点不同
                       //(exp(ik*R)-1)/R部分
                    green = std::exp(local_constant.ik * distance) / distance;
                    green0 = FComplex(1.0 / distance);

                    integral_tmp1 = dot(list_basis_fn_f[i], list_basis_fn_s[j]) * green - j_dot_j * green0;
                    integral_tmp1 -= divj_dot_divj * (green - green0);
                }
                // integral_result += (src.GetElement())->GetGaussWeight()[i] *
                // (field.GetElement())->GetGaussWeight()[j] * integral_tmp1 * src.GetJacobi()*
                // field.GetJacobi();
                // with complete basis function
                integral_result += list_gauss_weight_f[i] * list_gauss_weight_s[j] * integral_tmp1;
            }

            //奇异部分，即计算1/R的积分
            one_Rth = IntegralOverROfTri((src.GetElement(idx_figure_s))->GetVertex(), list_gauss_point_f[i]);
            integral_result += list_gauss_weight_f[i] *
                (((j_dot_j - divj_dot_divj) * one_Rth + j_dot_j * z0_for_loss) / jacob_f); //+/jacob_f;
        }
        return integral_result;
    }

    //自由空间非奇异情况精确处理
    // A* = integral(exp(local_constant.ik*R)/R * Gauss_Weight)
    // Ad = integral(r*G*r' * Gauss_Weight)
    // A = A*-Ad
    FComplex
        InteractionZRWG::IntegralNonsingularFigureRWGFreespace(const RWG &field, const RWG &src,
            int idx_figure_f, int idx_figure_s,
            const internal::LocalConstant &local_constant)
    {
        int i, j;
        FComplex integral_result, integral_tmp1, green;
        real distance = 0;
        real div_basis_fn_s = 0, div_basis_fn_f = 0;
        const std::vector<real> &list_gauss_weight_s = (src.GetElement(idx_figure_s))->GetGaussWeight();
        const std::vector<real> &list_gauss_weight_f = (field.GetElement(idx_figure_f))->GetGaussWeight();
        const std::vector<arma::vec3> &list_gauss_point_s = (src.GetElement(idx_figure_s))->GetGaussPoint();
        const std::vector<arma::vec3> &list_gauss_point_f = (field.GetElement(idx_figure_f))->GetGaussPoint();
        int num_gauss_nodes_s = list_gauss_point_s.size();
        int num_gauss_nodes_f = list_gauss_point_f.size();

        //提取rwg内部需要的信息：basis fn List、div basis fn List、图形的高斯点、高斯权重List
        std::vector<arma::vec3> list_basis_fn_s = src.GetRho(idx_figure_s);
        std::vector<arma::vec3> list_basis_fn_f = field.GetRho(idx_figure_f);

        //三角形高斯点的div_basis_fn都相同
        div_basis_fn_s = src.GetDivergence(idx_figure_s) / local_constant.k0;
        div_basis_fn_f = field.GetDivergence(idx_figure_f) / local_constant.k0;

        //对所有高斯点求积分
        for (i = 0; i < num_gauss_nodes_s; i++)
            for (j = 0; j < num_gauss_nodes_f; j++) {
                distance = arma::norm(list_gauss_point_s[i] - list_gauss_point_f[j], 2);

                // Green's Function
                green = std::exp(local_constant.ik * distance) / distance;

                // A*
                integral_tmp1 = dot(list_basis_fn_s[i], list_basis_fn_f[j]) * green;

                // Ad
                integral_tmp1 -= div_basis_fn_s * green * div_basis_fn_f;

                // integral_result += list_gauss_weight_s[i] * list_gauss_weight_f[j] * integral_tmp1 *
                // (src.GetElement())->GetJacobi()*(field.GetElement())->GetJacobi();	//complete basis function
                integral_result += list_gauss_weight_s[i] * list_gauss_weight_f[j] * integral_tmp1;
            }
        integral_result = integral_result;
        return integral_result;
    }


    // GND奇异情况处理
    //奇异情况将Green函数分成两部分进行计算：
    // Green = (exp(ik*R)-1)/R + 1/R
    //然后分别对两部分进行积分
    FComplex
        InteractionZRWG::IntegralSingularFigureRWGGNDLayer(const RWG &field, const RWG &src,
            int idx_figure_f, int idx_figure_s,
            const internal::LocalConstant &local_constant)
    {
        int i, j;
        FComplex integral_result, integral_tmp1, one_Rth, green0, z0_for_loss, GA_divj, G_reflect_divj, j_tmp2, j_tmp3;
        arma::cx_vec GA_j(3), G_reflect_j(3);
        real div_basis_fn_s = 0, div_basis_fn_f = 0;
        real jacob_f = 0, distance = 0, divj_dot_divj = 0, j_dot_j = 0, reflect_distance = 0;
        const std::vector<real> &list_gauss_weight_s = (src.GetElement(idx_figure_s))->GetGaussWeight();
        const std::vector<real> &list_gauss_weight_f = (field.GetElement(idx_figure_f))->GetGaussWeight();
        const std::vector<arma::vec3> &list_gauss_point_s = (src.GetElement(idx_figure_s))->GetGaussPoint();
        const std::vector<arma::vec3> &list_gauss_point_f = (field.GetElement(idx_figure_f))->GetGaussPoint();
        int num_gauss_nodes_s = list_gauss_point_s.size();
        int num_gauss_nodes_f = list_gauss_point_f.size();
        real omega = 2 * pi * local_constant.frequency;

        FComplex _i_we0 = -FComplex(0, 1) / omega / epsilon0;
        FComplex _iwu0 = -FComplex(0, 1) * omega * miu0;

        //提取rwg内部需要的信息：basis fn List、div basis fn List、图形的高斯点、高斯权重List
        const std::vector<arma::vec3> &list_basis_fn_s = src.GetRho(idx_figure_s);
        const std::vector<arma::vec3> &list_basis_fn_f = field.GetRho(idx_figure_f);

        //三角形高斯点的div_basis_fn都相同
        div_basis_fn_s = src.GetDivergence(idx_figure_s) / local_constant.k0;
        div_basis_fn_f = field.GetDivergence(idx_figure_f) / local_constant.k0;
        jacob_f = field.GetJacobi(idx_figure_f);
        divj_dot_divj = div_basis_fn_s * div_basis_fn_f;
        for (i = 0; i < num_gauss_nodes_f; i++) {
            j_dot_j = dot(list_basis_fn_s[i], list_basis_fn_f[i]);
            for (j = 0; j < num_gauss_nodes_s; j++) {
                // R and reflect_R
                distance = arma::norm(list_gauss_point_f[i] - list_gauss_point_s[j], 2);
                reflect_distance =
                    ReflectDistance(list_gauss_point_f[i], list_gauss_point_s[j], internal::Axis::Z_AXIS);

                // G_reflect_j and G_reflect_divj
                GReflectInterpolation(reflect_distance, local_constant, internal::Axis::Z_AXIS, G_reflect_j);
                G_reflect_divj = -std::exp(local_constant.ik * reflect_distance) / reflect_distance;
                if (distance < 1e-12) { //如果两个高斯点近似相同  //?精度要求
                    integral_tmp1 = (arma::dot(list_basis_fn_f[i], list_basis_fn_s[j]) - divj_dot_divj) * local_constant.ik;
                }

                else { //两个高斯点不同
                       //(exp(ik*R)-1)/R部分
                       // GA_j and GA_divj and green0
                    GAInterpolation(distance, local_constant, GA_j);
                    green0 = FComplex(1.0 / distance);
                    GA_divj = std::exp(local_constant.ik * distance) / distance;

                    integral_tmp1 = list_basis_fn_f[i][0] * GA_j[0] * list_basis_fn_s[j][0] +
                        list_basis_fn_f[i][1] * GA_j[1] * list_basis_fn_s[j][1] +
                        list_basis_fn_f[i][2] * GA_j[2] * list_basis_fn_s[j][2] - j_dot_j * green0;

                    integral_tmp1 -= divj_dot_divj * (GA_divj - green0);
                }
                // integral_result += (src.GetElement())->GetGaussWeight()[i] *
                // (field.GetElement())->GetGaussWeight()[j] * integral_tmp1 * src.GetJacobi()*
                // field.GetJacobi();
                // with complete basis function
                integral_tmp1 += list_basis_fn_f[i][0] * G_reflect_j[0] * list_basis_fn_s[j][0] +
                    list_basis_fn_f[i][1] * G_reflect_j[1] * list_basis_fn_s[j][1] +
                    list_basis_fn_f[i][2] * G_reflect_j[2] * list_basis_fn_s[j][2];
                integral_tmp1 -= divj_dot_divj * (G_reflect_divj);
                integral_result += list_gauss_weight_f[i] * list_gauss_weight_s[j] * integral_tmp1;
            }

            //奇异部分，即计算1/R的积分
            one_Rth = IntegralOverROfTri((src.GetElement(idx_figure_s))->GetVertex(), list_gauss_point_f[i]);
            integral_result += list_gauss_weight_f[i] *
                (((j_dot_j - divj_dot_divj) * one_Rth + j_dot_j * z0_for_loss) / jacob_f); //+/jacob_f;
        }
        return integral_result;
    }

    // GND非奇异情况精确处理
    // A* = integral(exp(local_constant.ik*R)/R * Gauss_Weight)
    // Ad = integral(r*G*r' * Gauss_Weight)
    // A = A*-Ad
    FComplex
        InteractionZRWG::IntegralNonsingularFigureRWGGNDLayer(const RWG &field, const RWG &src,
            int idx_figure_f, int idx_figure_s,
            const internal::LocalConstant &local_constant)
    {
        int i, j;
        FComplex integral_result, integral_tmp1, GA_divj, G_reflect_divj, Gphi;
        arma::cx_vec GA_j(3), G_reflect_j(3);
        real distance = 0, reflect_distance = 0;
        real div_basis_fn_s = 0, div_basis_fn_f = 0;
        const std::vector<real> &list_gauss_weight_s = (src.GetElement(idx_figure_s))->GetGaussWeight();
        const std::vector<real> &list_gauss_weight_f = (field.GetElement(idx_figure_f))->GetGaussWeight();
        const std::vector<arma::vec3> &list_gauss_point_s = (src.GetElement(idx_figure_s))->GetGaussPoint();
        const std::vector<arma::vec3> &list_gauss_point_f = (field.GetElement(idx_figure_f))->GetGaussPoint();
        int num_gauss_nodes_s = list_gauss_point_s.size();
        int num_gauss_nodes_f = list_gauss_point_f.size();

        //提取rwg内部需要的信息：basis fn List、div basis fn List、图形的高斯点、高斯权重List
        const std::vector<arma::vec3> &list_basis_fn_s = src.GetRho(idx_figure_s);
        const std::vector<arma::vec3> &list_basis_fn_f = field.GetRho(idx_figure_f);

        //三角形高斯点的div_basis_fn都相同
        div_basis_fn_s = src.GetDivergence(idx_figure_s) / local_constant.k0;
        div_basis_fn_f = field.GetDivergence(idx_figure_f) / local_constant.k0;

        //对所有高斯点求积分
        for (i = 0; i < num_gauss_nodes_s; i++)
            for (j = 0; j < num_gauss_nodes_f; j++) {
                // R and reflect R
                distance = arma::norm(list_gauss_point_s[i] - list_gauss_point_f[j], 2);
                reflect_distance =
                    ReflectDistance(list_gauss_point_f[j], list_gauss_point_s[i], internal::Axis::Z_AXIS);

                // Green's Function
                // GA_j and G_reflect_j
                GAInterpolation(distance, local_constant, GA_j);
                GReflectInterpolation(reflect_distance, local_constant, internal::Axis::Z_AXIS, G_reflect_j);

                // A*
                integral_tmp1 = list_basis_fn_f[j][0] * (GA_j[0] + G_reflect_j[0]) * list_basis_fn_s[i][0] +
                    list_basis_fn_f[j][1] * (GA_j[1] + G_reflect_j[1]) * list_basis_fn_s[i][1] +
                    list_basis_fn_f[j][2] * (GA_j[2] + G_reflect_j[2]) * list_basis_fn_s[i][2];


                GA_divj = GA_j[0];
                G_reflect_divj = G_reflect_j[0];

                Gphi = GA_divj + G_reflect_divj;
                integral_tmp1 -= div_basis_fn_s * Gphi * div_basis_fn_f;

                integral_result += list_gauss_weight_s[i] * list_gauss_weight_f[j] * integral_tmp1;
            }
        integral_result = integral_result;
        return integral_result;
    }
}
