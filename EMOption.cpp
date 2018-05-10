#include "EMOption.h"
//#include "MeshRule.h"
namespace MoM
{

	EMOption::EMOption(const std::vector<double> &frequency, const std::vector<std::string> &feed_ports,
		const std::vector<FComplex> &ref_impedance, const Unit unit, const SpaceType space_type,
		 const KernelMethod kernel_method, const FeedType feed_type,
		const FreqType freq_type, const MeshType mesh_type,
		const std::vector<std::string> &output_content, const std::string &output_file_name)
		: frequency(frequency), feed_ports(feed_ports), ref_impedance(ref_impedance), unit(unit), space_type(space_type),
		kernel_method(kernel_method), feed_type(feed_type), freq_type(freq_type), mesh_type(mesh_type),
		 output_content(output_content),
		output_file_name(output_file_name)
	{
		// internal arguments initial

		this->frequency = convertFreq2(this->frequency, this->freq_type);
		this->data_format = RI;
		this->freq_unit_out = HZ;
		switch (mesh_type) {
		case TRIANGLE: {
			this->integral_approach = RWG;
			break;
		}
		}


		if (output_content.size() == 0)
			printf("No output type is chose");
		this->port_mat_type = UNKOWN_PORT_TYPE;
		for (auto content : this->output_content) {
			if (content == "S_Parameter") {
				this->port_mat_type = static_cast<PortMatType>(static_cast<int>(this->port_mat_type) + 1);
			}
			else if (content == "Y_Parameter") {
				this->port_mat_type = static_cast<PortMatType>(static_cast<int>(this->port_mat_type) + 2);
			}
			else if (content == "Z_Parameter") {
				this->port_mat_type = static_cast<PortMatType>(static_cast<int>(this->port_mat_type) + 4);
			}
		}
		if (static_cast<int>(this->port_mat_type) > 7) {
			printf("same parameters were input!");
			this->port_mat_type = SYZ;
		}
	}
}