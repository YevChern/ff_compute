#include "pteros/analysis/task_plugin.h"
#include "pteros/analysis/trajectory_reader.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

#include "cmath"
#include <fstream>
#include <iomanip>
#include "FFMaps.h"
#include "main.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

const double FOUR_PI_SQUARED = 157.91367041742973790135185599802L;

// Auxilary functions for FF calculations

double getXrayStrength(std::string name, float q, std::map<char, vector<double> > map, float charge) {
    std::map<char,vector<double> >::const_iterator constants_itr = map.find(char(name[0]));
    float fraction = 0.0;
    double fi = 0.0;
    vector<double> constants = constants_itr->second;
    for (int i = 0; i < 4; ++i){
        fi += constants.at(i)*exp(-constants.at(i+4)*q*q/FOUR_PI_SQUARED);
        // Add up total number of electrons of the atom
        fraction += constants.at(i);
    }
    fi += constants.at(8);
    fraction += constants.at(8);
    // Determine scaling factor due to partial charge
    fraction /= fraction + charge;
    fi *= fraction;
    return fi;
}

double getIonXrayStrength(std::string name, float q, std::map<std::string, vector<double> > map, float charge) {
    std::map<string,vector<double> >::const_iterator constants_itr = map.find(name);
    float fraction = 0.0;
    double fi = 0.0;
    vector<double> constants = constants_itr->second;
    for (int i = 0; i < 4; ++i){
        fi += constants.at(i)*exp(-constants.at(i+4)*q*q/FOUR_PI_SQUARED);
        // Add up total number of electrons of the atom
        fraction += constants.at(i);
    }
    fi += constants.at(8);
    fraction += constants.at(8);
    // Determine scaling factor due to partial charge
    fraction /= fraction + charge;
    fi *= fraction;
    return fi;
}

// Charge delocalization-corrected x-ray scattering strength for water
double getWaterXrayStrength(std::string name, float q, std::map<char, vector<double> > map){
    float delta = 0.22;                     // nm^-1
    if (char(name[0]) == 'H') {
        float alpha_H = -0.48;
        double fi = getXrayStrength("H", q, map) * (1.0 + alpha_H*exp(-q*q/(2*delta*delta)));
        return fi;
    } else if (char(name[0]) == 'O') {
        float alpha_O = 0.12;
        double fi = getXrayStrength("O", q, map) * (1.0 + alpha_O*exp(-q*q/(2*delta*delta)));
        return fi;
    } else {
        throw pteros::Pteros_error("getWaterXrayStrength got the 'name' != 'O' or 'H'");
    }
}

double getNeutronStrength(std::string name, std::map<char, double> map){
    std::map<char, double>::const_iterator neutronSL_itr = map.find(char(name[0]));
    if (neutronSL_itr==map.end())
        throw pteros::Pteros_error("Warning: Atom type " + name + " not found in neutron atomic FF map!");

    double fi = neutronSL_itr->second;
    return fi;
}

// Same as getNeutronStrength, but for ions (as we can't rely on a single char to define a type)
double getIonNeutronStrength(std::string name, std::map<std::string, double> map){
    std::map<string, double>::const_iterator neutronSL_itr = map.find(name);
    if (neutronSL_itr==map.end())
        throw pteros::Pteros_error("Warning: Atom type " + name + " not found in neutron atomic FF map!");

    double fi = neutronSL_itr->second;
    return fi;
}

double getWaterNeutronStrength(std::string name, float deut, std::map<char, double > map){
    if (char(name[0]) == 'H') {
        double fi = getNeutronStrength("H", map)*(1 - deut) + getNeutronStrength("D", map)*deut;
        return fi;
    } else if (char(name[0]) == 'O') {
        double fi = getNeutronStrength("O", map);
        return fi;
    } else {
        throw pteros::Pteros_error("getWaterNeutronStrength got the 'name' != 'O' or 'H'");
    }
}


TASK_PARALLEL(FF_compute)
public:
    string help() override {
        return  "Purpose:\n"
                "\tPut purpose of your plugin here\n"
                "Output:\n"
                "\tDescription of its output\n"
                "Options:\n"
                "\tAny options";
    }
protected:
    void before_spawn() override {
        // Initialize constants
        cutoff = options("cutoff").as_float();
        w_dens = options("w_dens").as_float();
        w_dens_sqr = options("w_dens_sqr").as_float();
        d_parts = options("d_parts").as_floats();

        // Get particle indices for water, origin, ions
        std::fstream file_water(options("water_ind").as_string(), std::ios_base::in);
        int buf = 0;
        while (file_water >> buf){
            particles_water.push_back(buf);
        }
        file_water.close();

        std::fstream file_origin(options("orig_ind").as_string(), std::ios_base::in);
        while (file_origin >> buf){
            particles_for_origin.push_back(buf);
        }
        file_origin.close();

        string file_ions_str = options("ions_ind","").as_string();
        if (file_ions_str.size()>0){
            std::fstream file_ions(file_ions_str, std::ios_base::in);
            while (file_ions >> buf){
                particles_ion.push_back(buf);
            }
            file_ions.close();
        }

        // Get charges
        string file_charges_str = options("charges","").as_string();
        if (file_charges_str.size()>0){
            std::fstream file_charges(file_charges_str, std::ios_base::in);
            float charge_buf = 0.0;
            while (file_charges >> charge_buf){
                charge.push_back(charge_buf);
            }
            file_charges.close();
        } else {
            for (int i=0; i<system.select_all().size(); ++i){
                charge.push_back(0.0);
            }
        }

        // Get exchangeable hydrogens
        string file_exch_h_str = options("exch_h","").as_string();
        if (file_exch_h_str.size()>0){
            std::fstream file_exch_h(file_exch_h_str, std::ios_base::in);
            while (file_exch_h >> buf){
                particles_exch_h.push_back(buf);
            }
            file_exch_h.close();
        }

        // Initialize q values (either from file or defined by range and step)
        // From file:
        float fact = options("-q_fact","1.0").as_float();
        string file_q_xray_str = options("q_xray_file","").as_string();
        if (file_q_xray_str.size()>0){
            std::fstream file_q_xray(file_q_xray_str, std::ios_base::in);
            float q_xray_buf = 0.0;
            while (file_q_xray >> q_xray_buf){
                qs_xray.push_back(q_xray_buf*fact);
            }
            file_q_xray.close();
        }
        string file_q_neutron_str = options("q_neutron_file","").as_string();
        if (file_q_neutron_str.size()>0){
            std::fstream file_q_neutron(file_q_neutron_str, std::ios_base::in);
            float q_neutron_buf = 0.0;
            while (file_q_neutron >> q_neutron_buf){
                qs_neutron.push_back(q_neutron_buf*fact);
            }
            file_q_neutron.close();
        }

        // Setup maps for the helpers getXray/NeutronStrength()
        // Xray constants map
        aff_constants = create_aff_map();
        // Xray constants map (ions)
        aff_constants_ions = create_aff_ion_map();
        // Neutron constants map
        neutronFFmap = create_nsld_map();
        // Neutron constants map (ions)
        neutronFFmap_ions = create_nsld_ion_map();

        // Calculate scattering strength for individual atoms
        // First for all atoms in the same way, next, taking water/ion into account
        for (int i=0; i<system.select_all().size(); ++i){
            // X-ray
            vector<float> xrayStrength_tmp(qs_xray.size(), 0.0);
            for (int j=0; j<qs_xray.size(); ++j){
                xrayStrength_tmp[j] = getXrayStrength(system.atom(i).name, qs_xray[j], aff_constants, charge[i]);
            }
            xrayAtomStrength.push_back(xrayStrength_tmp);
            // Neutron
            vector<float> neutronStrength_tmp(d_parts.size(), 0.0);
            for (int j=0; j<d_parts.size(); ++j){
                neutronStrength_tmp[j] = getNeutronStrength(system.atom(i).name, neutronFFmap);
            }
            neutronAtomStrength.push_back(neutronStrength_tmp);
        }
        // Water:
        for (int i=0; i<particles_water.size(); ++i){
            // X-ray
            vector<float> xrayStrength_tmp(qs_xray.size(), 0.0);
            for (int j=0; j<qs_xray.size(); ++j){
                xrayStrength_tmp[j] = getWaterXrayStrength(system.atom(particles_water[i]).name, qs_xray[j], aff_constants);
            }
            xrayAtomStrength[particles_water[i]] = xrayStrength_tmp;
            // Neutron
            vector<float> neutronStrength_tmp(d_parts.size(), 0.0);
            for (int j=0; j<d_parts.size(); ++j){
                neutronStrength_tmp[j] = getWaterNeutronStrength(system.atom(particles_water[i]).name, d_parts[j], neutronFFmap);
            }
            neutronAtomStrength[particles_water[i]] = neutronStrength_tmp;
        }
        // Exchangeable hydrogens:
        for (int i=0; i<particles_exch_h.size(); ++i){
            // Only neutron strength is affected
            vector<float> neutronStrength_tmp(d_parts.size(), 0.0);
            for (int j=0; j<d_parts.size(); ++j){
                neutronStrength_tmp[j] = getWaterNeutronStrength("H", d_parts[j], neutronFFmap);
            }
            neutronAtomStrength[particles_exch_h[i]] = neutronStrength_tmp;
        }
        // Ions:
        for (int i=0; i<particles_ion.size(); ++i){
            // X-ray
            vector<float> xrayStrength_tmp(qs_xray.size(), 0.0);
            for (int j=0; j<qs_xray.size(); ++j){
                xrayStrength_tmp[j] = getIonXrayStrength(system.atom(particles_ion[i]).name, qs_xray[j], aff_constants_ions, charge[particles_ion[i]]);
            }
            xrayAtomStrength[particles_ion[i]] = xrayStrength_tmp;
            // Neutron
            vector<float> neutronStrength_tmp(d_parts.size(), 0.0);
            for (int j=0; j<d_parts.size(); ++j){
                neutronStrength_tmp[j] = getIonNeutronStrength(system.atom(particles_ion[i]).name, neutronFFmap_ions);
            }
            neutronAtomStrength[particles_ion[i]] = neutronStrength_tmp;
        }
    }

    void pre_process() override {
        // Xray initialize FF components
        A_real_xray.reserve(qs_xray.size());
        A_complex_xray.reserve(qs_xray.size());
        B_real_xray.reserve(qs_xray.size());
        A_sqr_xray.reserve(qs_xray.size());
        B_sqr_xray.reserve(qs_xray.size());
        for (int i=0; i<qs_xray.size(); ++i){
            A_real_xray.push_back(0.0);
            A_complex_xray.push_back(0.0);
            B_real_xray.push_back(0.0);
            A_sqr_xray.push_back(0.0);
            B_sqr_xray.push_back(0.0);
        }
        // Neutron
        A_real_neutron.reserve(qs_neutron.size()*d_parts.size());
        A_complex_neutron.reserve(qs_neutron.size()*d_parts.size());
        B_real_neutron.reserve(qs_neutron.size()*d_parts.size());
        A_sqr_neutron.reserve(qs_neutron.size()*d_parts.size());
        B_sqr_neutron.reserve(qs_neutron.size()*d_parts.size());
        for (int i=0; i<qs_neutron.size()*d_parts.size(); ++i){
            A_real_neutron.push_back(0.0);
            A_complex_neutron.push_back(0.0);
            B_real_neutron.push_back(0.0);
            A_sqr_neutron.push_back(0.0);
            B_sqr_neutron.push_back(0.0);
        }
    }

    void process_frame(const Frame_info &info) override {
        // Get B
        // X-ray
        for (int i=0; i<qs_xray.size(); ++i){
            const float wXrayStrength = 2.0 * getXrayStrength("H", qs_xray[i], aff_constants) * (1.0 + (-0.48)*exp(-qs_xray[i]*qs_xray[i]/(2.0*0.22*0.22))) +
                                              getXrayStrength("O", qs_xray[i], aff_constants) * (1.0 + 0.12*exp(-qs_xray[i]*qs_xray[i]/(2.0*0.22*0.22)));
            B_real_xray[i] += 2.0 * w_dens * system.box(0).extent(0) * system.box(0).extent(1) * wXrayStrength * sin(qs_xray[i]*cutoff) / qs_xray[i];
            B_sqr_xray[i] += w_dens_sqr * (2.0*system.box(0).extent(0) * system.box(0).extent(1) * wXrayStrength * sin(qs_xray[i]*cutoff) / qs_xray[i])
                                        * (2.0*system.box(0).extent(0) * system.box(0).extent(1) * wXrayStrength * sin(qs_xray[i]*cutoff) / qs_xray[i]);
        }
        // Neutron
        for (int i=0; i<qs_neutron.size(); ++i){
            // A_neutron and B_neutron components are of a size qs_neutron.size()*d_parts.size() as we have to calculate neutron FF for each deuteration value
            for (int d=0; d<d_parts.size(); ++d){
                const float w_neutr_scatt_streng = getNeutronStrength("O", neutronFFmap) + 2.0*(d_parts[d] * getNeutronStrength("D", neutronFFmap) +
                                                                                                (1.0 - d_parts[d]) * getNeutronStrength("H", neutronFFmap));
                B_real_neutron[d*qs_neutron.size() + i] += 2.0 * w_dens * system.box(0).extent(0) * system.box(0).extent(1) * w_neutr_scatt_streng * sin(qs_neutron[i]*cutoff) / qs_neutron[i];
                B_sqr_neutron[d*qs_neutron.size() + i] += w_dens_sqr * (2.0*system.box(0).extent(0) * system.box(0).extent(1) * w_neutr_scatt_streng * sin(qs_neutron[i]*cutoff) / qs_neutron[i])
                                                                     * (2.0*system.box(0).extent(0) * system.box(0).extent(1) * w_neutr_scatt_streng * sin(qs_neutron[i]*cutoff) / qs_neutron[i]);
            }
        }

        // Get bilayer COM coordinates
        float orig = 0.0;
        float total_m = 0.0;
        for (int i=0; i<particles_for_origin.size(); ++i){
            total_m += system.atom(particles_for_origin[i]).mass;
            orig += system.atom(particles_for_origin[i]).mass * system.xyz(particles_for_origin[i])(2);
        }
        orig /= total_m;

        // Calculate FF components
        vector<double> A_real_xray_current(A_real_xray.size(), 0.0), A_complex_xray_current(A_complex_xray.size(), 0.0),
                       A_real_neutron_current(A_real_neutron.size(), 0.0), A_complex_neutron_current(A_complex_neutron.size(), 0.0);
        for (int k=0; k<system.select_all().size(); ++k){
            // Wrap the coordinates of the atom around bilayer COM
            float z = system.xyz(k)(2) - orig;
            if (fabs(z) > system.box(0).extent(2)/2.0) {
                if (z<0.0) {
                    z = z - system.box(0).extent(2)*floor((z - 0.5*(system.box(0).extent(2)))/system.box(0).extent(2));
                }
                if (z>=0.0) {
                    z = z - system.box(0).extent(2)*floor((z + 0.5*(system.box(0).extent(2)))/system.box(0).extent(2));
                }
            }
            if (fabs(z)<=cutoff){
                // X-ray
                for (int i=0; i<qs_xray.size(); ++i){
                    A_real_xray_current[i] += xrayAtomStrength[k][i] * cos(qs_xray[i]*z);
                    A_complex_xray_current[i] += xrayAtomStrength[k][i] * sin(qs_xray[i]*z);
                }
                // Neutron
                for (int d=0; d<d_parts.size(); ++d){
                    for (int i=0; i<qs_neutron.size(); ++i){
                        A_real_neutron_current[d*qs_neutron.size() + i] += neutronAtomStrength[k][d] * cos(qs_neutron[i]*z);
                        A_complex_neutron_current[d*qs_neutron.size() + i] += neutronAtomStrength[k][d] * sin(qs_neutron[i]*z);
                    }
                }
            }
        }
        for (int i=0; i<qs_xray.size(); ++i){
            A_real_xray[i] += A_real_xray_current[i];
            A_complex_xray[i] += A_real_xray_current[i];
        }
        for (int d=0; d<d_parts.size(); ++d){
            for (int i=0; i<qs_neutron.size(); ++i){
                A_real_neutron[d*qs_neutron.size() + i] += A_real_neutron_current[d*qs_neutron.size() + i];
                A_complex_neutron[d*qs_neutron.size() + i] += A_complex_neutron_current[d*qs_neutron.size() + i];
            }
        }

        for (int i=0; i<qs_xray.size(); ++i){
            A_sqr_xray[i] += A_real_xray_current[i]*A_real_xray_current[i] + A_complex_xray_current[i]*A_complex_xray_current[i];
        }
        for (int d=0; d<d_parts.size(); ++d){
            for (int i=0; i<qs_neutron.size(); ++i){
                A_sqr_neutron[d*qs_neutron.size() + i] += A_real_neutron_current[d*qs_neutron.size() + i]*A_real_neutron_current[d*qs_neutron.size() + i]
                                                        + A_complex_neutron_current[d*qs_neutron.size() + i]*A_complex_neutron_current[d*qs_neutron.size() + i];
            }
        }

        log->info("Frame " + to_string(info.valid_frame));
    }

    void post_process(const Frame_info &info) override {

    }

    void collect_data(const std::vector<std::shared_ptr<Task_base>>& tasks, int n_frames) override {
        // Collect data and normalize
        vector<double> A_real_xray_out(qs_xray.size(), 0.0), A_complex_xray_out(qs_xray.size(), 0.0), A_sqr_xray_out(qs_xray.size(), 0.0),
                       B_real_xray_out(qs_xray.size(), 0.0), B_sqr_xray_out(qs_xray.size(), 0.0),
                       A_real_neutron_out(qs_neutron.size()*d_parts.size(), 0.0), A_complex_neutron_out(qs_neutron.size()*d_parts.size(), 0.0),
                       A_sqr_neutron_out(qs_neutron.size()*d_parts.size(), 0.0),
                       B_real_neutron_out(qs_neutron.size()*d_parts.size(), 0.0), B_sqr_neutron_out(qs_neutron.size()*d_parts.size(), 0.0);
        for (int i=0;i<tasks.size();++i){
            for (int j=0; j<qs_xray.size(); ++j){
                A_real_xray_out[j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->A_real_xray[j];
                A_complex_xray_out[j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->A_complex_xray[j];
                A_sqr_xray_out[j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->A_sqr_xray[j];
                B_real_xray_out[j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->B_real_xray[j];
                B_sqr_xray_out[j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->B_sqr_xray[j];
            }
            for (int d=0; d<d_parts.size(); ++d){
                for (int j=0; j<qs_neutron.size(); ++j){
                    A_real_neutron_out[d*qs_neutron.size() + j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->A_real_neutron[d*qs_neutron.size() + j];
                    A_complex_neutron_out[d*qs_neutron.size() + j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->A_complex_neutron[d*qs_neutron.size() + j];
                    A_sqr_neutron_out[d*qs_neutron.size() + j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->A_sqr_neutron[d*qs_neutron.size() + j];
                    B_real_neutron_out[d*qs_neutron.size() + j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->B_real_neutron[d*qs_neutron.size() + j];
                    B_sqr_neutron_out[d*qs_neutron.size() + j] += (1.0/n_frames) * std::static_pointer_cast<FF_compute>(tasks[i])->B_sqr_neutron[d*qs_neutron.size() + j];
                }
            }
        }
        // Output
        string out_pref_xray = options("out_pref", "").as_string();
        if (out_pref_xray.size() == 0) out_pref_xray = "xray";
        ofstream out_file_xray;
        out_file_xray.open(string(out_pref_xray+".xff"));
        for (int i=0; i<qs_xray.size(); ++i){
            const double F_total = sqrt(abs(A_real_xray_out[i]*A_real_xray_out[i] + A_complex_xray_out[i]*A_complex_xray_out[i]
                                        + B_real_xray_out[i]*B_real_xray_out[i] - 2.0*A_real_xray_out[i]*B_real_xray_out[i]
                                        + A_sqr_xray_out[i] - A_real_xray_out[i]*A_real_xray_out[i] - A_complex_xray_out[i]*A_complex_xray_out[i]
                                        - B_sqr_xray_out[i] + B_real_xray_out[i]*B_real_xray_out[i]));
            out_file_xray << qs_xray[i] << "    " << F_total << endl;
        }
        out_file_xray.close();
        for (int d=0; d<d_parts.size(); ++d){
            string out_pref_neutron = options("out_pref", "").as_string();
            if (out_pref_neutron.size() == 0) out_pref_neutron = "neutron";
            ofstream out_file_neutron;
            std::stringstream stream;
            stream << std::fixed << std::setprecision(2) << d_parts[d];
            string out_name_neutron = out_pref_neutron + stream.str() +"D.nff";
            out_file_neutron.open(out_name_neutron);
            for (int i=0; i<qs_neutron.size(); ++i){
                const double F_total = sqrt(abs(A_real_neutron_out[d*qs_neutron.size() + i]*A_real_neutron_out[d*qs_neutron.size() + i]
                                          + A_complex_neutron_out[d*qs_neutron.size() + i]*A_complex_neutron_out[d*qs_neutron.size() + i]
                                          + B_real_neutron_out[d*qs_neutron.size() + i]*B_real_neutron_out[d*qs_neutron.size() + i]
                                          - 2.0*A_real_neutron_out[d*qs_neutron.size() + i]*B_real_neutron_out[d*qs_neutron.size() + i]
                                          + A_sqr_neutron_out[d*qs_neutron.size() + i]
                                          - A_real_neutron_out[d*qs_neutron.size() + i]*A_real_neutron_out[d*qs_neutron.size() + i]
                                          - A_complex_neutron_out[d*qs_neutron.size() + i]*A_complex_neutron_out[d*qs_neutron.size() + i]
                                          - B_sqr_neutron_out[d*qs_neutron.size() + i]
                                          + B_real_neutron_out[d*qs_neutron.size() + i]*B_real_neutron_out[d*qs_neutron.size() + i]));
                out_file_neutron << qs_neutron[i] << "    " << F_total << endl;
            }
            out_file_neutron.close();
        }

        log->info("Finished");
    }

private:
    vector<double> A_real_xray, A_real_neutron, B_real_xray, B_real_neutron, A_complex_xray, A_complex_neutron, A_sqr_xray, A_sqr_neutron, B_sqr_xray, B_sqr_neutron;
    vector<int> particles_for_origin;
    vector<int> particles_water;
    vector<int> particles_ion;
    vector<int> particles_exch_h;
    vector<float> charge;

    vector<float> qs_xray, qs_neutron;
    float cutoff;
    float w_dens, w_dens_sqr;
    vector<float> d_parts;

    map<char, vector<double> > aff_constants;
    map<std::string, vector<double> > aff_constants_ions;
    map<char, double> neutronFFmap;
    map<std::string, double> neutronFFmap_ions;

    vector<vector<float> > xrayAtomStrength, neutronAtomStrength;

};

int main(int argc, char** argv){
    try {
        Options options;
        parse_command_line(argc,argv,options);
        Trajectory_reader engine(options);
        auto task = new FF_compute(options);
        engine.add_task(task);
        cout << "-------------------------------------------------------------" << endl;
        cout << "  This is stand-alone Pteros analysis plugin 'FF_compute'" << endl;
        cout << "-------------------------------------------------------------" << endl;
        if(!options.has("f") && !options.has("help")){
            cout << "Usage:" << endl;
            cout << "\tpteros_FF_compute -f <files> <task options>" << endl;
            cout << "\n\tFor specific task options use '-help task'" << endl;
            cout << "\tFor trajectory processing options use '-help traj'" << endl;
            cout << "\tFor all available options use '-help all' or just '-help'" << endl;
            return 1;
        }
        if(options.has("help")){
            string help = options("help","").as_string();
            if(help=="traj"){
                cout << engine.help() << endl;
            } else if(help=="task"){
                cout << task->help() << endl;
            } else {
                cout << task->help() << endl << endl;
                cout << engine.help() << endl;
            }
            return 1;
        }
        engine.run();
    } catch (const std::exception& e) {
        LOG()->error(e.what());
    } catch(...) {
        LOG()->error("Unknown error");
    }
}

