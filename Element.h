#ifndef _ELEMENT_H_
#define _ELEMENT_H_
#include <vector>
#include <armadillo>
#include "Typedef.h"
#include "EMOption.h"
#include "uniform_func.h"


namespace MoM
{

	class FPort;
	class FPoint;

	class Element
	{
	private:
		// index of the element
		int index;

		// Coordinate of vertex of the element
		std::vector<arma::vec3> element_vertex;
		arma::vec3 center_point;

		// if on port, pointer of the port
		const FPort *port;

		std::vector<arma::vec3> gauss_points;
		std::vector<real> gauss_weight;

		std::vector<arma::vec3> diff_gauss_points;
		std::vector<real> diff_gauss_weight;

		static bool init_flag;
		static double sym_gauss_weight_tet[9][100];
		static double sym_gauss_coord_tet[9][4][100];
		static int num_sym_gauss_nodes_tet[9];
		static void GenerateWeightsAndCoordinatesForSymGaussianQuadrature();

	public:
		//Element();

		Element(int index, const std::vector<arma::vec3> &element_vertex, const FPort *port = nullptr);

		~Element();

		std::vector<bool> visited;

		// Get the values of the member variables
		int GetIndex() const;
		void SetIndex(int SetIndex) { this->index = SetIndex; }
		std::vector<arma::vec3> &GetVertex();
		arma::vec3 GetCenter() { return this->center_point; }
		void Show() 
		{
			for (int i = 0; i < static_cast<int>(this->element_vertex.size());i++) {
				std::cout << "Vertex[" << i << "]: (" << element_vertex[i][0] << " , " << element_vertex[i][1] << " , " << element_vertex[i][2] << ")" << std::endl;
			}
			std::cout << std::endl;
		}
		const FPort *GetPort() const;
		void SetPort(const FPort *port);
		const std::vector<arma::vec3> &GetGaussPoint() const;
		const std::vector<real> &GetGaussWeight() const;
		const std::vector<arma::vec3> &GetDiffGaussPoint() const;
		const std::vector<real> &GetDiffGaussWeight() const;
		void convertLengthUnit(Unit from, Unit to);
		arma::vec3 GetPointByParam(real u0, real u1, real u2) const;
		arma::vec3 GetPointByParam(real u0, real u1, real u2, real u3) const;
	};
}

#endif // _ELEMENT_H_INCLUDED_