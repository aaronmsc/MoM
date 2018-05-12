#include "RWG.h"


namespace mom
{

    using namespace std;
    const double DINF = numeric_limits<double>::max();
    const double NDINF = -DINF;

    bool RWG::CheckSubscript(int subscript) const
    {
        if (IsHalfbasis()) {
            return subscript == 0;
        }
        else {
            return subscript <= 1 && subscript >= 0;
        }
    }

    RWG::RWG() {}

    RWG::RWG(const Basis &basis, vector<Element> &element_list)
    {
        if (basis[0] < 0) {
            is_halfbasis = true;
            element[0] = &element_list[basis[1]];
            element[1] = nullptr;
        }
        else if (basis[1] < 0) {
            is_halfbasis = true;
            element[0] = &element_list[basis[0]];
            element[1] = nullptr;
        }
        else {
            is_halfbasis = false;
            element[0] = &element_list[basis[0]];
            element[1] = &element_list[basis[1]];
        }
        // parent = nullptr;
    }



    void RWG::CreateCompleteRWG(int original_index)
    {
        if (is_halfbasis) {
            MSG("this method is only for complete basis");
            return;
        }

        // if (element[0]->GetPort() != NULL || (!is_halfbasis && element[1]->GetPort() != NULL))
        // 	cout << "Position1:" << endl;
        // if (element[0]->GetPort() != NULL) {
        // 	Msg::Info("[%d] elem1 port=%d", element[0]->GetIndex(), element[0]->GetPort()->GetIndex());
        // }
        // if (!is_halfbasis && element[1]->GetPort() != NULL) {
        // 	Msg::Info("[%d] elem2 port=%d", element[1]->GetIndex(), element[1]->GetPort()->GetIndex());
        // }
        // get the original points
        bool is_distinct[2][3];
        memset(is_distinct, true, sizeof(is_distinct));

        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                if (arma::approx_equal(element[0]->GetVertex()[i], element[1]->GetVertex()[j], "absdiff", 1e-13)) {
                    is_distinct[0][i] = false;
                    is_distinct[1][j] = false;
                    break;
                }
            }
        }
        for (size_t i = 0; i < 2; i++) {

            bool original_flag = false;
            for (size_t j = 0; j < 3; j++) {
                if (is_distinct[i][j]) {
                    original_point[i] = j;
                    // save the edge which is visited
                    element[i]->visited[j] = true;
                    original_flag = true;
                    break;
                }
            }
            if (original_flag == false) {
                MSG("original point not exist!");
                for (arma::vec3 &p : element[0]->GetVertex()) {
                    MSG(" element 0 :%f %f %f", p[0], p[1], p[2]);
                }
                for (arma::vec3 &p : element[1]->GetVertex()) {
                    MSG(" element 1 :%f %f %f", p[0], p[1], p[2]);
                }
                exit(0);
            }
        }

        // get the length of the common edge

        common_length = arma::norm(element[0]->GetVertex()[(original_point[0] + 1) % 3] -
            element[0]->GetVertex()[(original_point[0] + 2) % 3], 2);
        this->index = original_index;
        this->current_direction[0] = 1;
        this->current_direction[1] = port_basis_flag == internal::AREA ? 1 : -1;
        rho_list[0] = this->CalculateRho(0);
        rho_list[1] = this->CalculateRho(1);
        diff_rho_list[0] = this->CalculateDiffRho(0);
        diff_rho_list[1] = this->CalculateDiffRho(1);
    }

    void RWG::CreateHalfRWG(int original_index)
    {
        if (!is_halfbasis) {
            MSG("this method is only for half basis");
            return;
        }

        // find the original point which is not visited
        for (size_t i = 0; i < 3; i++) {
            if (!element[0]->visited[i]) {
                element[0]->visited[i] = true;
                original_point[0] = i;
                break;
            }
        }

        // get the length of the common edge
        common_length = arma::norm(element[0]->GetVertex()[(original_point[0] + 1) % 3] -
            element[0]->GetVertex()[(original_point[0] + 2) % 3], 2);
        this->index = original_index;
        this->index = original_index;
        this->current_direction[0] = 1;
        this->current_direction[1] = -1;
        rho_list[0] = this->CalculateRho(0);
        diff_rho_list[0] = this->CalculateDiffRho(0);
    }

    bool RWG::IsHalfbasis() const { return this->is_halfbasis; }

    Element *RWG::GetElement(int subscript) const
    {
        assert(CheckSubscript(subscript));
        return element[subscript];
    }

    real RWG::GetCurrentDirection(int subscript) const { return this->current_direction[subscript]; }

    real RWG::GetCommonLength() const { return this->common_length; }

    int RWG::GetOriginalPoint(int subscript) const
    {
        assert(CheckSubscript(subscript));
        return this->original_point[subscript];
    }
    arma::vec3 RWG::GetElementCentre(int subscript) const
    {
        assert(CheckSubscript(subscript));
        return GetPointByParam(subscript, 0.3333333333, 0.3333333333, 0.3333333334);
    }

    real RWG::GetDivergence(int subscript) const
    {
        assert(CheckSubscript(subscript));

        return GetCurrentDirection(subscript) * 2;
    }

    vector<arma::vec3> RWG::CalculateRho(int subscript)
    {
        assert(CheckSubscript(subscript));
        vector<arma::vec3> rho;
        const vector<arma::vec3> &gauss_points = element[subscript]->GetGaussPoint();
        // Msg::Info("start cal rho %d", subscript);
        for (auto it = gauss_points.begin(); it != gauss_points.end(); it++) {
            arma::vec3 temp =
                (*it - element[subscript]->GetVertex()[original_point[subscript]]) * GetCurrentDirection(subscript);
            rho.push_back(temp);
        }
        // Msg::Info("start cal rho %d", subscript);

        return rho;
    }
    vector<arma::vec3> RWG::CalculateDiffRho(int subscript)
    {
        assert(CheckSubscript(subscript));
        vector<arma::vec3> rho;
        const vector<arma::vec3> &gauss_points = element[subscript]->GetDiffGaussPoint();
        // Msg::Info("start cal rho %d", subscript);
        for (auto it = gauss_points.begin(); it != gauss_points.end(); it++) {
            arma::vec3 temp =
                (*it - element[subscript]->GetVertex()[original_point[subscript]]) * GetCurrentDirection(subscript);
            rho.push_back(temp);
        }
        // Msg::Info("start cal rho %d", subscript);

        return rho;
    }
    arma::vec3 RWG::GetGaussianIntegral(int subscript, const vector<arma::vec3> &rho) const
    {
        assert(CheckSubscript(subscript));

        arma::vec3 ret({ 0, 0, 0 });
        const vector<real> &gauss_weight = element[subscript]->GetGaussWeight();
        for (auto i = 0; i != rho.size(); i++) {
            ret = ret + (rho[i] * gauss_weight[i]);
        }
        return ret;
    }
    arma::vec3 RWG::GetDiffGaussianIntegral(int subscript, const vector<arma::vec3> &rho) const
    {
        assert(CheckSubscript(subscript));

        arma::vec3 ret({ 0, 0, 0 });
        const vector<real> &gauss_weight = element[subscript]->GetDiffGaussWeight();
        for (auto i = 0; i != rho.size(); i++) {
            ret = ret + (rho[i] * gauss_weight[i]);
        }
        return ret;
    }


    arma::vec3 RWG::GetPointByParam(int subscript, real u0, real u1, real u2) const
    {
        assert(CheckSubscript(subscript));
        auto vertex = element[subscript]->GetVertex();
        return vertex[0] * u0 + vertex[1] * u1 + vertex[2] * u2;
    }



    real RWG::GetJacobi(int subscript) const
    {
        assert(CheckSubscript(subscript));
        vector<arma::vec3> vertex = element[subscript]->GetVertex();
        return arma::norm(arma::cross((vertex[0] - vertex[1]), (vertex[0] - vertex[2])), 2);
    }
    real RWG::GetRadius(int subscript) const
    {
        assert(CheckSubscript(subscript));
        return arma::norm(GetElementCentre(subscript) - element[subscript]->GetVertex().front(), 2);
    }
}
