#include "Element.h"

namespace MoM
{

	double Element::sym_gauss_weight_tet[9][100];
	double Element::sym_gauss_coord_tet[9][4][100];
	int Element::num_sym_gauss_nodes_tet[9];
	bool Element::init_flag = (GenerateWeightsAndCoordinatesForSymGaussianQuadrature(), true);

	//Element::Element() {}

	Element::Element(const int index, const std::vector<arma::vec3> &element_vertex, const FPort *port)
	{
		this->index = index;
		for (auto i = element_vertex.begin(); i != element_vertex.end(); i++) {
			this->element_vertex.push_back(*i);
		}
		this->port = port;

		//
		int basis_count = 0;
		if (element_vertex.size() == 4)
			basis_count = 4;
		else if (element_vertex.size() == 3)
			basis_count = 3;

		for (auto vertex: element_vertex) {
			this->center_point += vertex;
		}
		this->center_point /= element_vertex.size();

		for (int i = 0; i < basis_count; i++)
			visited.push_back(false);

		int gauss_accuracy = 3;

		if (element_vertex.size() == 4) {
			int gauss_num = num_sym_gauss_nodes_tet[gauss_accuracy];
			for (int i = 0; i < gauss_num; i++) {
				gauss_points.push_back(this->GetPointByParam(
					sym_gauss_coord_tet[gauss_accuracy][0][i], sym_gauss_coord_tet[gauss_accuracy][1][i],
					sym_gauss_coord_tet[gauss_accuracy][2][i], sym_gauss_coord_tet[gauss_accuracy][3][i]));
				gauss_weight.push_back(sym_gauss_weight_tet[gauss_accuracy][i]);
			}

		}
		else if (element_vertex.size() == 3) {

			gauss_points.push_back(this->GetPointByParam(0.666666666666667, 0.166666666666667, 0.166666666666667));
			gauss_points.push_back(this->GetPointByParam(0.166666666666667, 0.666666666666667, 0.166666666666667));
			gauss_points.push_back(this->GetPointByParam(0.166666666666667, 0.166666666666667, 0.666666666666667));
			for (int i = 0; i < 3; i++)
				gauss_weight.push_back(1.0 / 6);

			// gauss_points.push_back(this->GetPointByParam(0.33333333, 0.33333333, 0.33333334));
			// gauss_points.push_back(this->GetPointByParam(3.0 / 5, 1.0 / 5, 1.0 / 5));
			// gauss_points.push_back(this->GetPointByParam(1.0 / 5, 3.0 / 5, 1.0 / 5));
			// gauss_points.push_back(this->GetPointByParam(1.0 / 5, 1.0 / 5, 3.0 / 5));
			// gauss_weight.push_back(-9.0 / 32);
			// gauss_weight.push_back(25.0 / 96);
			// gauss_weight.push_back(25.0 / 96);
			// gauss_weight.push_back(25.0 / 96);
			
			// diff_gauss_points.push_back(this->GetPointByParam(0.333333333, 0.333333333, 0.333333334));
			// diff_gauss_points.push_back(this->GetPointByParam(0.059715871789770, 0.470142064105115,
			// 0.470142064105115)); diff_gauss_points.push_back(this->GetPointByParam(0.797426985353088,
			// 0.101286507323456, 0.101286507323456)); diff_gauss_weight.push_back(0.225);
			// diff_gauss_weight.push_back(0.132394152788506);
			// diff_gauss_weight.push_back(0.125939180544827);
		}
		return;
	}

	Element::~Element() {}

	int Element::GetIndex() const { return this->index; }

	std::vector<arma::vec3> &Element::GetVertex() { return this->element_vertex; }
	const FPort *Element::GetPort() const { return this->port; }
	void Element::SetPort(const FPort *port) { this->port = port; }
	const std::vector<arma::vec3> &Element::GetGaussPoint() const { return gauss_points; }
	const std::vector<real> &Element::GetGaussWeight() const { return gauss_weight; }
	const std::vector<arma::vec3> &Element::GetDiffGaussPoint() const { return diff_gauss_points; }
	const std::vector<real> &Element::GetDiffGaussWeight() const { return diff_gauss_weight; }
	inline arma::vec3 Element::GetPointByParam(real u0, real u1, real u2) const
	{
		return element_vertex[0] * u0 + element_vertex[1] * u1 + element_vertex[2] * u2;
	}

	inline arma::vec3 Element::GetPointByParam(real u0, real u1, real u2, real u3) const
	{
		return element_vertex[0] * u0 + element_vertex[1] * u1 + element_vertex[2] * u2 + element_vertex[3] * u3;
	}

	void Element::convertLengthUnit(Unit from, Unit to)
	{
		for (auto &point : this->element_vertex) {
			EMOption::convertGeoUnit(from, to, point);
		}

		for (auto &point : this->gauss_points) {
			EMOption::convertGeoUnit(from, to, point);
		}
		for (auto &point : this->diff_gauss_points) {
			EMOption::convertGeoUnit(from, to, point);
		}
	}

	

	void Element::GenerateWeightsAndCoordinatesForSymGaussianQuadrature()
	{
		double sw, w, c[4];
		int i, j, k, l, m;
		// Generate the weights and coordinates of degree N=2.
		w = 0.25;
		c[0] = 0.5854101966249685;
		c[1] = 0.1381966011250105;
		c[2] = 0.1381966011250105;
		c[3] = 0.1381966011250105;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[2][i] = w;
			sym_gauss_coord_tet[2][0][i] = c[i % 4];
			sym_gauss_coord_tet[2][1][i] = c[(i + 1) % 4];
			sym_gauss_coord_tet[2][2][i] = c[(i + 2) % 4];
			sym_gauss_coord_tet[2][3][i] = c[(i + 3) % 4];
		}
		num_sym_gauss_nodes_tet[2] = 4;

		// Generate the weights and coordinates of degree N=3.
		w = -0.8;
		c[0] = 0.25;
		c[1] = 0.25;
		c[2] = 0.25;
		c[3] = 0.25;
		sym_gauss_weight_tet[3][0] = w;
		sym_gauss_coord_tet[3][0][0] = c[0 % 4];
		sym_gauss_coord_tet[3][1][0] = c[(0 + 1) % 4];
		sym_gauss_coord_tet[3][2][0] = c[(0 + 2) % 4];
		sym_gauss_coord_tet[3][3][0] = c[(0 + 3) % 4];

		w = 0.45;
		c[0] = 0.5;
		c[1] = 0.1666666666666667;
		c[2] = 0.1666666666666667;
		c[3] = 0.1666666666666667;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[3][i + 1] = w;
			sym_gauss_coord_tet[3][0][i + 1] = c[i % 4];
			sym_gauss_coord_tet[3][1][i + 1] = c[(i + 1) % 4];
			sym_gauss_coord_tet[3][2][i + 1] = c[(i + 2) % 4];
			sym_gauss_coord_tet[3][3][i + 1] = c[(i + 3) % 4];
		}
		num_sym_gauss_nodes_tet[3] = 5;

		// Generate the weights and coordinates of degree N=4.
		w = 0.5037379410012282e-1;
		c[0] = 0.7716429020672371;
		c[1] = 0.7611903264425430e-1;
		c[2] = 0.7611903264425430e-1;
		c[3] = 0.7611903264425430e-1;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[4][i] = w;
			sym_gauss_coord_tet[4][0][i] = c[i % 4];
			sym_gauss_coord_tet[4][1][i] = c[(i + 1) % 4];
			sym_gauss_coord_tet[4][2][i] = c[(i + 2) % 4];
			sym_gauss_coord_tet[4][3][i] = c[(i + 3) % 4];
		}

		w = 0.6654206863329239e-1;
		c[0] = 0.1197005277978019;
		c[1] = 0.7183164526766925e-1;
		c[2] = 0.4042339134672644;
		c[3] = 0.4042339134672644;
		i = 4;
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = 0; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[4][i] = w;
								sym_gauss_coord_tet[4][j][i] = c[0];
								sym_gauss_coord_tet[4][k][i] = c[1];
								sym_gauss_coord_tet[4][l][i] = c[2];
								sym_gauss_coord_tet[4][m][i] = c[3];
								break;
							}
						}
					i++;
				}
			}
		num_sym_gauss_nodes_tet[4] = 16;

		// Generate the weights and coordinates of degree N=5.
		w = 0.1884185567365411;
		c[0] = 0.25;
		c[1] = 0.25;
		c[2] = 0.25;
		c[3] = 0.25;
		sym_gauss_weight_tet[5][0] = w;
		sym_gauss_coord_tet[5][0][0] = c[0 % 4];
		sym_gauss_coord_tet[5][1][0] = c[(0 + 1) % 4];
		sym_gauss_coord_tet[5][2][0] = c[(0 + 2) % 4];
		sym_gauss_coord_tet[5][3][0] = c[(0 + 3) % 4];

		w = 0.6703858372604275e-1;
		c[0] = 0.7316369079576180;
		c[1] = 0.8945436401412733e-1;
		c[2] = 0.8945436401412733e-1;
		c[3] = 0.8945436401412733e-1;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[5][i + 1] = w;
			sym_gauss_coord_tet[5][0][i + 1] = c[i % 4];
			sym_gauss_coord_tet[5][1][i + 1] = c[(i + 1) % 4];
			sym_gauss_coord_tet[5][2][i + 1] = c[(i + 2) % 4];
			sym_gauss_coord_tet[5][3][i + 1] = c[(i + 3) % 4];
		}

		w = 0.4528559236327399e-1;
		c[0] = 0.1325810999384657;
		c[1] = 0.2454003792903000e-1;
		c[2] = 0.4214394310662522;
		c[3] = 0.4214394310662522;
		i = 5;
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = 0; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[5][i] = w;
								sym_gauss_coord_tet[5][j][i] = c[0];
								sym_gauss_coord_tet[5][k][i] = c[1];
								sym_gauss_coord_tet[5][l][i] = c[2];
								sym_gauss_coord_tet[5][m][i] = c[3];
								break;
							}
						}
					i++;
				}
			}
		num_sym_gauss_nodes_tet[5] = 17;

		// Generate the weights and coordinates of degree N=6.
		w = 0.9040129046014750e-1;
		c[0] = 0.25;
		c[1] = 0.25;
		c[2] = 0.25;
		c[3] = 0.25;
		sym_gauss_weight_tet[6][0] = w;
		sym_gauss_coord_tet[6][0][0] = c[0 % 4];
		sym_gauss_coord_tet[6][1][0] = c[(0 + 1) % 4];
		sym_gauss_coord_tet[6][2][0] = c[(0 + 2) % 4];
		sym_gauss_coord_tet[6][3][0] = c[(0 + 3) % 4];

		w = 0.1911983427899124e-1;
		c[0] = 0.8277192480479295;
		c[1] = 0.5742691731735683e-1;
		c[2] = 0.5742691731735683e-1;
		c[3] = 0.5742691731735683e-1;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[6][i + 1] = w;
			sym_gauss_coord_tet[6][0][i + 1] = c[i % 4];
			sym_gauss_coord_tet[6][1][i + 1] = c[(i + 1) % 4];
			sym_gauss_coord_tet[6][2][i + 1] = c[(i + 2) % 4];
			sym_gauss_coord_tet[6][3][i + 1] = c[(i + 3) % 4];
		}

		w = 0.4361493840666568e-1;
		c[0] = 0.5135188412556341e-1;
		c[1] = 0.4860510285706072;
		c[2] = 0.2312985436519147;
		c[3] = 0.2312985436519147;
		i = 5;
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = 0; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[6][i] = w;
								sym_gauss_coord_tet[6][j][i] = c[0];
								sym_gauss_coord_tet[6][k][i] = c[1];
								sym_gauss_coord_tet[6][l][i] = c[2];
								sym_gauss_coord_tet[6][m][i] = c[3];
								break;
							}
						}
					i++;
				}
			}

		w = 0.2581167596199161e-1;
		c[0] = 0.2967538129690260;
		c[1] = 0.6081079894015281;
		c[2] = 0.4756909881472290e-1;
		c[3] = 0.4756909881472290e-1;
		i = 17;
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = 0; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[6][i] = w;
								sym_gauss_coord_tet[6][j][i] = c[0];
								sym_gauss_coord_tet[6][k][i] = c[1];
								sym_gauss_coord_tet[6][l][i] = c[2];
								sym_gauss_coord_tet[6][m][i] = c[3];
								break;
							}
						}
					i++;
				}
			}
		num_sym_gauss_nodes_tet[6] = 29;

		// N=5,6,7,8. These data are digested from P. Keast, "Moderate-degree tetrahedral quadrature formulas," Computer
		// Methods
		// in Applied Mechanics and Engineering, vo. 55, pp339-348, 1986.

		// Generate the weights and coordinates of degree N=4.
		sw = -0.131555555555555550e-1 + 4 * 0.762222222222222222e-2 + 6 * 0.248888888888888880e-1;
		w = -0.131555555555555550e-1;
		c[0] = 0.250000000000000000e+0;
		c[1] = 0.250000000000000000e+0;
		c[2] = 0.250000000000000000e+0;
		c[3] = 0.250000000000000000e+0;
		sym_gauss_weight_tet[4][0] = w / sw;
		sym_gauss_coord_tet[4][0][0] = c[0 % 4];
		sym_gauss_coord_tet[4][1][0] = c[(0 + 1) % 4];
		sym_gauss_coord_tet[4][2][0] = c[(0 + 2) % 4];
		sym_gauss_coord_tet[4][3][0] = c[(0 + 3) % 4];

		w = 0.762222222222222222e-2;
		c[0] = 0.714285714285714285e-1;
		c[1] = 0.714285714285714285e-1;
		c[2] = 0.714285714285714285e-1;
		c[3] = 0.785714285714285714e+0;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[4][i + 1] = w / sw;
			sym_gauss_coord_tet[4][0][i + 1] = c[i % 4];
			sym_gauss_coord_tet[4][1][i + 1] = c[(i + 1) % 4];
			sym_gauss_coord_tet[4][2][i + 1] = c[(i + 2) % 4];
			sym_gauss_coord_tet[4][3][i + 1] = c[(i + 3) % 4];
		}

		w = 0.248888888888888880e-1;
		c[0] = 0.399403576166799219e+0;
		c[1] = 0.399403576166799219e+0;
		c[2] = 0.100596423833200785e+0;
		c[3] = 0.100596423833200785e+0;
		i = 5;
		for (j = 0; j < 3; j++)
			for (k = j; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = l; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[4][i] = w / sw;
								sym_gauss_coord_tet[4][j][i] = c[0];
								sym_gauss_coord_tet[4][k][i] = c[1];
								sym_gauss_coord_tet[4][l][i] = c[2];
								sym_gauss_coord_tet[4][m][i] = c[3];
								i++;
								break;
							}
						}
				}
			}
		num_sym_gauss_nodes_tet[4] = 11;

		// Generate the weights and coordinates of degree N=5.
		sw = 1 * 0.302836780970891856e-1 + 4 * 0.602678571428571597e-2 + 4 * 0.116452490860289742e-1 +
			6 * 0.109491415613864534e-1;

		w = 0.302836780970891856e-1;
		c[0] = 0.250000000000000000e+0;
		c[1] = 0.250000000000000000e+0;
		c[2] = 0.250000000000000000e+0;
		c[3] = 0.250000000000000000e+0;
		sym_gauss_weight_tet[5][0] = w / sw;
		sym_gauss_coord_tet[5][0][0] = c[0 % 4];
		sym_gauss_coord_tet[5][1][0] = c[(0 + 1) % 4];
		sym_gauss_coord_tet[5][2][0] = c[(0 + 2) % 4];
		sym_gauss_coord_tet[5][3][0] = c[(0 + 3) % 4];

		w = 0.602678571428571597e-2;
		c[0] = 0.333333333333333333e+0;
		c[1] = 0.333333333333333333e+0;
		c[2] = 0.333333333333333333e+0;
		c[3] = 0.000000000000000000e+0;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[5][i + 1] = w / sw;
			sym_gauss_coord_tet[5][0][i + 1] = c[i % 4];
			sym_gauss_coord_tet[5][1][i + 1] = c[(i + 1) % 4];
			sym_gauss_coord_tet[5][2][i + 1] = c[(i + 2) % 4];
			sym_gauss_coord_tet[5][3][i + 1] = c[(i + 3) % 4];
		}

		w = 0.116452490860289742e-1;
		c[0] = 0.909090909090909091e-1;
		c[1] = 0.909090909090909091e-1;
		c[2] = 0.909090909090909091e-1;
		c[3] = 0.727272727272727273e+0;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[5][i + 5] = w / sw;
			sym_gauss_coord_tet[5][0][i + 5] = c[i % 4];
			sym_gauss_coord_tet[5][1][i + 5] = c[(i + 1) % 4];
			sym_gauss_coord_tet[5][2][i + 5] = c[(i + 2) % 4];
			sym_gauss_coord_tet[5][3][i + 5] = c[(i + 3) % 4];
		}

		i = 9;
		w = 0.109491415613864534e-1;
		c[0] = 0.665501535736642813e-1;
		c[1] = 0.665501535736642813e-1;
		c[2] = 0.433449846426335728e+0;
		c[3] = 0.433449846426335728e+0;
		for (j = 0; j < 3; j++)
			for (k = j; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = l; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[5][i] = w / sw;
								sym_gauss_coord_tet[5][j][i] = c[0];
								sym_gauss_coord_tet[5][k][i] = c[1];
								sym_gauss_coord_tet[5][l][i] = c[2];
								sym_gauss_coord_tet[5][m][i] = c[3];
								i++;
								break;
							}
						}
				}
			}
		num_sym_gauss_nodes_tet[5] = 15;

		// Generate the weights and coordinates of degree N=6.
		sw = 4 * 0.665379170969464506e-2 + 4 * 0.167953517588677620e-2 + 4 * 0.922619692394239843e-2 +
			12 * 0.803571428571428248e-2;

		w = 0.665379170969464506e-2;
		c[0] = 0.214602871259151684e+0;
		c[1] = 0.214602871259151684e+0;
		c[2] = 0.214602871259151684e+0;
		c[3] = 0.356191386222544953e+0;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[6][i] = w / sw;
			sym_gauss_coord_tet[6][0][i] = c[i % 4];
			sym_gauss_coord_tet[6][1][i] = c[(i + 1) % 4];
			sym_gauss_coord_tet[6][2][i] = c[(i + 2) % 4];
			sym_gauss_coord_tet[6][3][i] = c[(i + 3) % 4];
		}

		w = 0.167953517588677620e-2;
		c[0] = 0.406739585346113397e-1;
		c[1] = 0.406739585346113397e-1;
		c[2] = 0.406739585346113397e-1;
		c[3] = 0.877978124396165982e+0;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[6][i + 4] = w / sw;
			sym_gauss_coord_tet[6][0][i + 4] = c[i % 4];
			sym_gauss_coord_tet[6][1][i + 4] = c[(i + 1) % 4];
			sym_gauss_coord_tet[6][2][i + 4] = c[(i + 2) % 4];
			sym_gauss_coord_tet[6][3][i + 4] = c[(i + 3) % 4];
		}

		w = 0.922619692394239843e-2;
		c[0] = 0.322337890142275646e+0;
		c[1] = 0.322337890142275646e+0;
		c[2] = 0.322337890142275646e+0;
		c[3] = 0.329863295731730594e-1;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[6][i + 8] = w / sw;
			sym_gauss_coord_tet[6][0][i + 8] = c[i % 4];
			sym_gauss_coord_tet[6][1][i + 8] = c[(i + 1) % 4];
			sym_gauss_coord_tet[6][2][i + 8] = c[(i + 2) % 4];
			sym_gauss_coord_tet[6][3][i + 8] = c[(i + 3) % 4];
		}

		w = 0.803571428571428248e-2;
		c[0] = 0.636610018750175299e-1;
		c[1] = 0.636610018750175299e-1;
		c[2] = 0.269672331458315867e+0;
		c[3] = 0.603005664791649076e+0;
		i = 12;
		for (j = 0; j < 4; j++)
			for (k = j; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = 0; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[6][i] = w / sw;
								sym_gauss_coord_tet[6][j][i] = c[0];
								sym_gauss_coord_tet[6][k][i] = c[1];
								sym_gauss_coord_tet[6][l][i] = c[2];
								sym_gauss_coord_tet[6][m][i] = c[3];
								i++;
							}
						}
				}
			}
		num_sym_gauss_nodes_tet[6] = 24;

		// Generate the weights and coordinates of degree N=7.
		sw = 6 * 0.970017636684296702e-3 + 1 * 0.182642234661087939e-1 + 4 * 0.105999415244141609e-1 -
			4 * 0.625177401143299494e-1 + 4 * 0.489142526307353653e-2 + 12 * 0.275573192239850917e-1;

		w = 0.105999415244141609e-1;
		c[0] = 0.782131923303186549e-1;
		c[1] = 0.782131923303186549e-1;
		c[2] = 0.782131923303186549e-1;
		c[3] = 0.765360423009044044e+0;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[7][i] = w / sw;
			sym_gauss_coord_tet[7][0][i] = c[i % 4];
			sym_gauss_coord_tet[7][1][i] = c[(i + 1) % 4];
			sym_gauss_coord_tet[7][2][i] = c[(i + 2) % 4];
			sym_gauss_coord_tet[7][3][i] = c[(i + 3) % 4];
		}

		w = -0.625177401143299494e-1;
		c[0] = 0.121843216663904411e+0;
		c[1] = 0.121843216663904411e+0;
		c[2] = 0.121843216663904411e+0;
		c[3] = 0.634470350008286765e+0;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[7][i + 4] = w / sw;
			sym_gauss_coord_tet[7][0][i + 4] = c[i % 4];
			sym_gauss_coord_tet[7][1][i + 4] = c[(i + 1) % 4];
			sym_gauss_coord_tet[7][2][i + 4] = c[(i + 2) % 4];
			sym_gauss_coord_tet[7][3][i + 4] = c[(i + 3) % 4];
		}

		w = 0.489142526307353653e-2;
		c[0] = 0.332539164446420554e+0;
		c[1] = 0.332539164446420554e+0;
		c[2] = 0.332539164446420554e+0;
		c[3] = 0.238250666073834549e-2;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[7][i + 8] = w / sw;
			sym_gauss_coord_tet[7][0][i + 8] = c[i % 4];
			sym_gauss_coord_tet[7][1][i + 8] = c[(i + 1) % 4];
			sym_gauss_coord_tet[7][2][i + 8] = c[(i + 2) % 4];
			sym_gauss_coord_tet[7][3][i + 8] = c[(i + 3) % 4];
		}

		w = 0.275573192239850917e-1;
		c[0] = 0.100000000000000000e+0;
		c[1] = 0.100000000000000000e+0;
		c[2] = 0.200000000000000000e+0;
		c[3] = 0.600000000000000000e+0;
		i = 12;
		for (j = 0; j < 4; j++)
			for (k = j; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = 0; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[7][i] = w / sw;
								sym_gauss_coord_tet[7][j][i] = c[0];
								sym_gauss_coord_tet[7][k][i] = c[1];
								sym_gauss_coord_tet[7][l][i] = c[2];
								sym_gauss_coord_tet[7][m][i] = c[3];
								i++;
							}
						}
				}
			}

		w = 0.182642234661087939e-1;
		c[0] = 0.250000000000000000e+0;
		c[1] = 0.250000000000000000e+0;
		c[2] = 0.250000000000000000e+0;
		c[3] = 0.250000000000000000e+0;
		sym_gauss_weight_tet[7][i] = w / sw;
		sym_gauss_coord_tet[7][0][i] = c[0 % 4];
		sym_gauss_coord_tet[7][1][i] = c[(0 + 1) % 4];
		sym_gauss_coord_tet[7][2][i] = c[(0 + 2) % 4];
		sym_gauss_coord_tet[7][3][i] = c[(0 + 3) % 4];

		i = 25;
		w = 0.970017636684296702e-3;
		c[0] = 0.500000000000000000e+0;
		c[1] = 0.500000000000000000e+0;
		c[2] = 0.000000000000000000e+0;
		c[3] = 0.000000000000000000e+0;
		for (j = 0; j < 3; j++)
			for (k = j; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = l; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[7][i] = w / sw;
								sym_gauss_coord_tet[7][j][i] = c[0];
								sym_gauss_coord_tet[7][k][i] = c[1];
								sym_gauss_coord_tet[7][l][i] = c[2];
								sym_gauss_coord_tet[7][m][i] = c[3];
								i++;
								break;
							}
						}
				}
			}

		num_sym_gauss_nodes_tet[7] = 31;

		// Generate the weights and coordinates of degree N=8.
		sw = -1 * 0.393270066412926145e-1 + 4 * 0.408131605934270525e-2 + 4 * 0.658086773304341943e-3 +
			6 * 0.438425882512284693e-2 + 6 * 0.138300638425098166e-1 + 12 * 0.424043742468372453e-2 +
			12 * 0.223873973961420164e-2;

		w = -0.393270066412926145e-1;
		c[0] = 0.250000000000000000e+0;
		c[1] = 0.250000000000000000e+0;
		c[2] = 0.250000000000000000e+0;
		c[3] = 0.250000000000000000e+0;
		i = 0;
		sym_gauss_weight_tet[8][i] = w / sw;
		sym_gauss_coord_tet[8][0][i] = c[0 % 4];
		sym_gauss_coord_tet[8][1][i] = c[(0 + 1) % 4];
		sym_gauss_coord_tet[8][2][i] = c[(0 + 2) % 4];
		sym_gauss_coord_tet[8][3][i] = c[(0 + 3) % 4];

		w = 0.408131605934270525e-2;
		c[0] = 0.127470936566639015e+0;
		c[1] = 0.127470936566639015e+0;
		c[2] = 0.127470936566639015e+0;
		c[3] = 0.617587190300082967e+0;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[8][i + 1] = w / sw;
			sym_gauss_coord_tet[8][0][i + 1] = c[i % 4];
			sym_gauss_coord_tet[8][1][i + 1] = c[(i + 1) % 4];
			sym_gauss_coord_tet[8][2][i + 1] = c[(i + 2) % 4];
			sym_gauss_coord_tet[8][3][i + 1] = c[(i + 3) % 4];
		}

		w = 0.658086773304341943e-3;
		c[0] = 0.320788303926322960e-1;
		c[1] = 0.320788303926322960e-1;
		c[2] = 0.320788303926322960e-1;
		c[3] = 0.903763508822103123e+0;
		for (i = 0; i < 4; i++) {
			sym_gauss_weight_tet[8][i + 5] = w / sw;
			sym_gauss_coord_tet[8][0][i + 5] = c[i % 4];
			sym_gauss_coord_tet[8][1][i + 5] = c[(i + 1) % 4];
			sym_gauss_coord_tet[8][2][i + 5] = c[(i + 2) % 4];
			sym_gauss_coord_tet[8][3][i + 5] = c[(i + 3) % 4];
		}

		i = 9;
		w = 0.438425882512284693e-2;
		c[0] = 0.497770956432810185e-1;
		c[1] = 0.497770956432810185e-1;
		c[2] = 0.450222904356718978e+0;
		c[3] = 0.450222904356718978e+0;
		for (j = 0; j < 3; j++)
			for (k = j; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = l; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[8][i] = w / sw;
								sym_gauss_coord_tet[8][j][i] = c[0];
								sym_gauss_coord_tet[8][k][i] = c[1];
								sym_gauss_coord_tet[8][l][i] = c[2];
								sym_gauss_coord_tet[8][m][i] = c[3];
								i++;
								break;
							}
						}
				}
			}

		i = 15;
		w = 0.138300638425098166e-1;
		c[0] = 0.183730447398549945e+0;
		c[1] = 0.183730447398549945e+0;
		c[2] = 0.316269552601450060e+0;
		c[3] = 0.316269552601450060e+0;
		for (j = 0; j < 3; j++)
			for (k = j; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = l; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[8][i] = w / sw;
								sym_gauss_coord_tet[8][j][i] = c[0];
								sym_gauss_coord_tet[8][k][i] = c[1];
								sym_gauss_coord_tet[8][l][i] = c[2];
								sym_gauss_coord_tet[8][m][i] = c[3];
								i++;
								break;
							}
						}
				}
			}

		w = 0.424043742468372453e-2;
		c[0] = 0.231901089397150906e+0;
		c[1] = 0.231901089397150906e+0;
		c[2] = 0.229177878448171174e-1;
		c[3] = 0.513280033360881072e+0;
		i = 21;
		for (j = 0; j < 4; j++)
			for (k = j; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = 0; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[8][i] = w / sw;
								sym_gauss_coord_tet[8][j][i] = c[0];
								sym_gauss_coord_tet[8][k][i] = c[1];
								sym_gauss_coord_tet[8][l][i] = c[2];
								sym_gauss_coord_tet[8][m][i] = c[3];
								i++;
							}
						}
				}
			}

		w = 0.223873973961420164e-2;
		c[0] = 0.379700484718286102e-1;
		c[1] = 0.379700484718286102e-1;
		c[2] = 0.730313427807538396e+0;
		c[3] = 0.193746475248804382e+0;
		i = 33;
		for (j = 0; j < 4; j++)
			for (k = j; k < 4; k++) {
				if (j != k) {
					for (l = 0; l < 4; l++)
						for (m = 0; m < 4; m++) {
							if (j != l && j != m && k != l && k != m && l != m) {
								sym_gauss_weight_tet[8][i] = w / sw;
								sym_gauss_coord_tet[8][j][i] = c[0];
								sym_gauss_coord_tet[8][k][i] = c[1];
								sym_gauss_coord_tet[8][l][i] = c[2];
								sym_gauss_coord_tet[8][m][i] = c[3];
								i++;
							}
						}
				}
			}
		num_sym_gauss_nodes_tet[8] = 45;
	}
}