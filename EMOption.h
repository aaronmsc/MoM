#ifndef _EMOPTION_H_
#define _EMOPTION_H_


#include <vector>
#include <algorithm>
#include <armadillo>
#include "Typedef.h"
#include "uniform_func.h"


namespace mom
{

    // Input Options
    enum CoordinatePlane { XY_PLANE = 0, YZ_PLANE, XZ_PLANE };
    enum Unit { NM = 0, UM, MM, M, MIL };
    enum SpaceType { FREE_SPACE = 0, GND_FREESPACE };
    enum AccurLevel { PCB_LOW = 0, PCB_STANDARD, PCB_HIGH, IC_LOW, IC_STANDARD, IC_HIGH };
    enum KernelMethod { MOM = 0 };
    enum IntegralApproach { RWG_BASIS = 0 };
    enum FeedType { VOLTAGE = 0, GAP };
    enum FreqType { LINEAR_STEP = 0, LOGARITHM, DISCRETE };
    enum MeshType { TRIANGLE = 0 };

    // Output Options
    enum FrequencyUnit { HZ = 0, KHZ, MHZ, GHZ };
    enum DataFormat { DB = 0, RI, MA };
    enum PortMatType { UNKOWN_PORT_TYPE = 0, S = 1, Y = 2, SY = 3, Z = 4, SZ = 5, YZ = 6, SYZ = 7 };

    struct EMOption {

        // Input
        Unit unit;
        SpaceType space_type;
        KernelMethod kernel_method;
        FeedType feed_type;
        FreqType freq_type;
        MeshType mesh_type;
        std::vector<real> frequency;
        std::vector<std::string> feed_ports;

        // Output
        std::vector<FComplex> ref_impedance;
        std::string output_file_name;
        std::vector<std::string> output_content;

        // internal
        IntegralApproach integral_approach;
        FrequencyUnit freq_unit_out;
        DataFormat data_format;
        int plot_point_num;
        PortMatType port_mat_type;


        EMOption() {}
        EMOption(const std::vector<double> &frequency, const std::vector<std::string> &feed_ports = { "net1.P1" },
            const std::vector<FComplex> &ref_impedance = { FComplex(50, 0) }, const Unit unit = MM,
            const SpaceType space_type = GND_FREESPACE,
            const KernelMethod kernel_method = MOM, const FeedType feed_type = VOLTAGE,
            const FreqType freq_type = LINEAR_STEP, const MeshType mesh_type = TRIANGLE,
            const std::vector<std::string> &output_content = { "S_Parameter" },
            const std::string &output_file_name = "simulation");

        // convert arbitrary freq to uniform freq
        // freq_def format:	LINEAR_STEP:[begin,end,step]
        //					DISCRETE   :[num1,num2,...]
        //					LOGARITHM  :[begin,end,num_of_freq] Formula:freq_list = [begin*power(end/begin,i/(n-1)) for i in
        // range(0,freq_def.size())]
        static std::vector<real> convertFreq2(std::vector<real> freq_def, const FreqType freq_type = LINEAR_STEP)
        {
            MSG("start convertFreq2");
            std::vector<real> freq_list;
            for (auto freq : freq_def) {
                if (freq == 0) {
                    MSG("frequency should not have been 0!");
                    exit(1);
                }
            }
            if (freq_type == LINEAR_STEP) { //传递上下界和步长的情况
                real tmp_freq;
                if (freq_def.size() != 3)
                    MSG("frequency define unknow");
                if (freq_def[0] > freq_def[1])
                    std::swap(freq_def[0], freq_def[1]);

                if (freq_def[2] < 0) {
                    tmp_freq = freq_def[1];
                    MSG("frequency step < 0! ");
                    do {
                        freq_list.push_back(tmp_freq);
                        tmp_freq += freq_def[2];
                    } while (tmp_freq > freq_def[0]);
                    if (freq_def[0] != freq_def[1] && tmp_freq == freq_def[0])
                        freq_list.push_back(tmp_freq);
                }
                else if (freq_def[2] > 0) {
                    tmp_freq = freq_def[0];
                    do {
                        freq_list.push_back(tmp_freq);
                        tmp_freq += freq_def[2];
                    } while (tmp_freq < freq_def[1]);
                    if (freq_def[0] != freq_def[1] && tmp_freq == freq_def[1])
                        freq_list.push_back(tmp_freq);
                }
                else {
                    freq_list.push_back(freq_def[0]);
                    freq_list.push_back(freq_def[1]);
                }
            }
            else if (freq_type == DISCRETE) { //传递频率点的情况
                for (auto freq : freq_def) {
                    freq_list.push_back(freq);
                }
            }
            else if (freq_type == LOGARITHM) { // logarithm
                if (freq_def.size() != 3 || freq_def[2] <= 0)
                    MSG("frequency defined unknown");
                if (freq_def[0] > freq_def[1])
                    std::swap(freq_def[0], freq_def[1]);
                for (size_t i = 0; i < freq_def[2]; i++) {
                    freq_list.push_back(freq_def[0] *
                        pow(freq_def[1] / freq_def[0], static_cast<double>(i) / (freq_def[2] - 1)));
                }
            }
            std::sort(freq_list.begin(), freq_list.end());
            return freq_list;
        }

        // freq Unit convert
        static double convertFreqUnit(double value, FrequencyUnit from = HZ, FrequencyUnit to = HZ)
        {
            if (from != to) {
                switch (from) {
                case HZ: {
                    switch (to) {
                    case KHZ:
                        value *= 1e-3;
                        break;
                    case MHZ:
                        value *= 1e-6;
                        break;
                    case GHZ:
                        value *= 1e-9;
                        break;
                    default:
                        break;
                    }
                    break;
                }

                case KHZ: {
                    switch (to) {
                    case HZ:
                        value *= 1e3;
                        break;
                    case MHZ:
                        value *= 1e-3;
                        break;
                    case GHZ:
                        value *= 1e-6;
                        break;
                    default:
                        break;
                    }
                    break;
                }

                case MHZ: {
                    switch (to) {
                    case HZ:
                        value *= 1e6;
                        break;
                    case KHZ:
                        value *= 1e3;
                        break;
                    case GHZ:
                        value *= 1e-3;
                        break;
                    default:
                        break;
                    }
                    break;
                }

                case GHZ: {
                    switch (to) {
                    case HZ:
                        value *= 1e9;
                        break;
                    case KHZ:
                        value *= 1e6;
                        break;
                    case MHZ:
                        value *= 1e3;
                        break;
                    default:
                        break;
                    }
                    break;
                }
                default:
                    break;
                }
            }
            return value;
        }

        // arma::vec3 convert
        static void convertGeoUnit(Unit from, Unit to, arma::vec3 &point)
        {
            point[0] = convertLengthUnit(point[0], from, to);
            point[1] = convertLengthUnit(point[1], from, to);
            point[2] = convertLengthUnit(point[2], from, to);
        }

        // length Unit convert
        static double convertLengthUnit(double value, Unit from = M, Unit to = M)
        {
            if (from != to) {
                switch (from) {
                case M: { // m To others unit
                    switch (to) {
                    case MIL:
                        value *= (1.0 / 2.54) * 1e5;
                        break;
                    case MM:
                        value *= 1e3;
                        break;
                    case UM:
                        value *= 1e6;
                        break;
                    case NM:
                        value *= 1e9;
                        break;
                    default:
                        break;
                    }
                    break;
                }

                case MIL: { // mil to other unit
                    switch (to) {
                    case M:
                        value *= 2.54e-5;
                        break;
                    case MM:
                        value *= 2.54e-2;
                        break;
                    case UM:
                        value *= 2.54e1;
                        break;
                    case NM:
                        value *= 2.54e4;
                        break;
                    default:
                        break;
                    }
                    break;
                }

                case MM: { // mm To others unit
                    switch (to) {
                    case M:
                        value *= 1e-3;
                        break;
                    case MIL:
                        value *= (1.0 / 2.54) * 1e2;
                        break;
                    case UM:
                        value *= 1e3;
                        break;
                    case NM:
                        value *= 1e6;
                        break;
                    default:
                        break;
                    }
                    break;
                }

                case UM: { // um To others unit
                    switch (to) {
                    case M:
                        value *= 1e-6;
                        break;
                    case MIL:
                        value *= (1.0 / 2.54) * 1e-1;
                        break;
                    case MM:
                        value *= 1e-3;
                        break;
                    case NM:
                        value *= 1e3;
                        break;
                    default:
                        break;
                    }
                    break;
                }

                case NM: { // nm To others unit
                    switch (to) {
                    case M:
                        value *= 1e-9;
                        break;
                    case MIL:
                        value *= (1.0 / 2.54) * 1e-4;
                        break;
                    case MM:
                        value *= 1e-6;
                        break;
                    case UM:
                        value *= 1e-3;
                        break;
                    default:
                        break;
                    }
                    break;
                }
                default:
                    break;
                }
            }
            return value;
        }
    };
}

#endif // EMOPTION_H_
