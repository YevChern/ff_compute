#include "FFMaps.h"

// Setup atomic scattering constants maps for the helpers getXray/NeutronStrength()

// Xray constants map
std::map<char, vector<double> > create_aff_map() {
    std::map<char, vector<double> > aff_constants;
    //build the vector to insert, then put in the map, do for each element
    //taken from: E.N. Maslen, A.G. Fox, M.A. O’Keefe Interpretation of diffracted intensities (6.1) A.J.C. Wilson, E. Prince (Eds.), International Table for Crystallography, Kluwer Academic Publishers, Dordrecht, The Netherlands (1999), pp. 547-584
    vector<double> hydrogen;
    hydrogen.push_back(0.493); //
    hydrogen.push_back(0.323);
    hydrogen.push_back(0.14);
    hydrogen.push_back(0.041);
    hydrogen.push_back(10.511);
    hydrogen.push_back(26.126);
    hydrogen.push_back(3.142);
    hydrogen.push_back(57.8);
    hydrogen.push_back(0.003);
    aff_constants.insert(make_pair('H', hydrogen));

    vector<double> deuterium;
    deuterium.push_back(0.493); //a1
    deuterium.push_back(0.323); //a2
    deuterium.push_back(0.14); //a3
    deuterium.push_back(0.041); //a4
    deuterium.push_back(10.511); //b1
    deuterium.push_back(26.126); //b2
    deuterium.push_back(3.142); //b3
    deuterium.push_back(57.8); //b4
    deuterium.push_back(0.003); //c
    aff_constants.insert(make_pair('D', deuterium));

    vector<double> carbon;
    carbon.push_back(2.31); //a1
    carbon.push_back(1.02); //a2
    carbon.push_back(1.589); //a3
    carbon.push_back(0.865); //a4
    carbon.push_back(20.844); //b1
    carbon.push_back(10.208); //b2
    carbon.push_back(0.569); //b3
    carbon.push_back(51.651); //b4
    carbon.push_back(0.216); //c
    aff_constants.insert(make_pair('C', carbon));

    vector<double> nitrogen;
    nitrogen.push_back(12.213); //a1
    nitrogen.push_back(3.132); //a2
    nitrogen.push_back(2.013); //a3
    nitrogen.push_back(1.166); //a4
    nitrogen.push_back(0.006); //b1
    nitrogen.push_back(9.893); //b2
    nitrogen.push_back(28.997); //b3
    nitrogen.push_back(0.583); //b4
    nitrogen.push_back(-11.524); //c
    aff_constants.insert(make_pair('N', nitrogen));

    vector<double> oxygen;
    oxygen.push_back(3.049); //a1
    oxygen.push_back(2.287); //a2
    oxygen.push_back(1.546); //a3
    oxygen.push_back(0.867); //a4
    oxygen.push_back(13.277); //b1
    oxygen.push_back(5.701); //b2
    oxygen.push_back(0.324); //b3
    oxygen.push_back(32.909); //b4
    oxygen.push_back(0.251); //c
    aff_constants.insert(make_pair('O', oxygen));

    vector<double> phosphorus;
    phosphorus.push_back(6.435); //a1
    phosphorus.push_back(4.179); //a2
    phosphorus.push_back(1.78); //a3
    phosphorus.push_back(1.491); //a4
    phosphorus.push_back(1.907); //b1
    phosphorus.push_back(27.157); //b2
    phosphorus.push_back(0.526); //b3
    phosphorus.push_back(68.164); //b4
    phosphorus.push_back(1.115); //c
    aff_constants.insert(make_pair('P', phosphorus));
    return aff_constants;
}

// Xray ions map
std::map<std::string, vector<double> > create_aff_ion_map() {
    std::map<std::string, vector<double> > aff_constants_ions;
    vector<double> potassium_plus;
    potassium_plus.push_back(7.9578); //a1
    potassium_plus.push_back(7.4917); //a2
    potassium_plus.push_back(6.359); //a3
    potassium_plus.push_back(1.1915); //a4
    potassium_plus.push_back(12.6331); //b1
    potassium_plus.push_back(0.7674); //b2
    potassium_plus.push_back(-0.002); //b3
    potassium_plus.push_back(31.9128); //b4
    potassium_plus.push_back(-4.9978); //c
    aff_constants_ions.insert(make_pair("K+", potassium_plus));

    vector<double> chloride_minus;
    chloride_minus.push_back(18.2915); //a1
    chloride_minus.push_back(7.2084); //a2
    chloride_minus.push_back(6.5337); //a3
    chloride_minus.push_back(2.3386); //a4
    chloride_minus.push_back(0.0066); //b1
    chloride_minus.push_back(1.1717); //b2
    chloride_minus.push_back(19.5424); //b3
    chloride_minus.push_back(60.4486); //b4
    chloride_minus.push_back(-16.378); //c
    aff_constants_ions.insert(make_pair("Cl-", chloride_minus));

    vector<double> sodium_plus;
    sodium_plus.push_back(3.2565); //a1
    sodium_plus.push_back(3.9362); //a2
    sodium_plus.push_back(1.3998); //a3
    sodium_plus.push_back(1.0032); //a4
    sodium_plus.push_back(2.6671); //b1
    sodium_plus.push_back(6.1153); //b2
    sodium_plus.push_back(0.2001); //b3
    sodium_plus.push_back(14.039); //b4
    sodium_plus.push_back(0.404); //c
    aff_constants_ions.insert(make_pair("Na+", sodium_plus));
    return aff_constants_ions;
}

// Neutron map
std::map<char, double> create_nsld_map() {
    std::map<char, double> neutronFFmap;
    neutronFFmap.insert(pair<char,double>('H',-3.7409E-5)); //taken from: "Koester, L., Nistier, W.: Z. Phys. A 272 (1975) 189."
    neutronFFmap.insert(pair<char,double>('D',6.67E-5));
    neutronFFmap.insert(pair<char,double>('C',6.6484E-5));  //taken from: "Koester, L., Nistier, W.: Z. Phys. A 272 (1975) 189."
    neutronFFmap.insert(pair<char,double>('N',9.36E-5));    //taken from: "Koester, L., Knopf, K., Waschkowski, W.: Z. Phys. A 277 (1976) 77."
    neutronFFmap.insert(pair<char,double>('O',5.805E-5));   //taken from: "Koester, L., Knopf, K., Waschkowski, W.: Z. Phys. A 292 (1979) 95."
    neutronFFmap.insert(pair<char,double>('P',5.17E-5));
    return neutronFFmap;
}

// Neutron map (ions)
std::map<std::string, double> create_nsld_ion_map() {
    std::map<std::string, double> neutronFFmap_ions;
    neutronFFmap_ions.insert(pair<std::string,double>("K+",3.67E-5));
    neutronFFmap_ions.insert(pair<std::string,double>("Cl-",9.5792E-5));
    neutronFFmap_ions.insert(pair<std::string,double>("Na+",3.63E-5));
    return neutronFFmap_ions;
}
