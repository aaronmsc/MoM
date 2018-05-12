#ifndef _MOM_H_
#define _MOM_H_
#ifndef MOM_HPP_
#define MOM_HPP_

#include "Constant.h"
#include "RWG.h"
#include <time.h>
#include <memory>

#include "EMOption.h"
#include "InteractionZRWG.h"

//#define _DEBUG_ENGINE_ALONE_

//可能出现的编译问题：
// 1.前置函数声明


#endif // MOM_HPP_
namespace mom {
    class MoM {
    public:
        virtual void Run(const real freq, const std::vector<arma::cx_vec> &F_list, std::vector<arma::cx_vec> &I_list) = 0;
        virtual ~MoM() {};
    };

    class MoMRWG :public MoM {
        InteractionZRWG m_inter_Z_;
    private:
        const std::vector<RWG> *m_rwg_list_;
        const EMOption *m_emoption_;

    public:
        MoMRWG(const std::vector<RWG>* rwg_list, const EMOption* emoption)
        {
            this->m_rwg_list_ = rwg_list;
            this->m_emoption_ = emoption;
            this->m_inter_Z_.Init(emoption->space_type);
        }
        void Run(const real freq, const std::vector<arma::cx_vec> &F_list, std::vector<arma::cx_vec> &I_list) {
            int num_Bas = m_rwg_list_->size();
            arma::cx_mat Z(num_Bas, num_Bas, arma::fill::zeros);


            FComplex _p4ie = pi4 * FComplex(0, 1) / _eta;





            //填充Z矩阵
            if (this->m_emoption_->space_type == FREE_SPACE || this->m_emoption_->space_type == GND_FREESPACE) {
                // #pragma omp parallel for
                for (int i = 0; i < num_Bas; i++) {
                    for (int j = i; j < num_Bas; j++) {
                        Z(i, j) = this->m_inter_Z_.Calculate((*m_rwg_list_)[j], (*m_rwg_list_)[i], freq);
                        if (j < 5 && i < 5)
                            MSG("No(%d,%d): (%.6f,%.6f)\n", i, j, (Z(i, j) / _p4ie).real(), (Z(i, j) / _p4ie).imag());
                    }
                }

                for (int i = 1; i < num_Bas; i++)
                    for (int j = 0; j < i; j++)
                        Z(i, j) = Z(j, i);
                // Z.block(0,0,7,7).print("Z(0-8):");
            }


#ifdef _DEBUG_EFIE_ALONE_
            FComplex sum, abs_sum;
            for (size_t i = 0; i < num_Bas; i++)
                for (size_t j = 0; j < num_Bas; j++) {
                    sum += std::abs(Z(i, j));
                }
            sum /= _p4ie;
            MSG("sumZ=(%lf,%lf)\n", sum.real(), sum.imag());

            Z.save(this->m_emoption_->output_file_name + "_" + "Z_mat");
            // Z.save(prefix_path + "Z_mat.txt" + std::to_string(freq));

#endif

            for (size_t i = 0; i < F_list.size(); i++) {
                arma::solve(I_list[i], Z, F_list[i]);
            }
            I_list[0].save("I_list.txt", arma::raw_ascii);
            std::complex<double> sumI;
            for (size_t i = 0; i < I_list[0].size(); i++) {
                sumI += I_list[0][i];
            }
            MSG("sumI:(%.6f,%.6f)\n", sumI.real(), sumI.imag());
        }
        ~MoMRWG() {}
    };
    class MoMFactory {
    public:
        std::unique_ptr<MoM> CreateMoM(const std::vector<RWG>* rwg_list, const EMOption* emoption) {
            return std::make_unique<MoMRWG>(rwg_list, emoption);
        }
    };
}
#endif
