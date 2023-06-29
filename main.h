#include <string>

double getXrayStrength(std::string name, float q, std::map<char, std::vector<double> > map, float charge=0.0);
double getIonXrayStrength(std::string name, float q, std::map<std::string, std::vector<double> > map, float charge);
double getWaterXrayStrength(std::string name, float q, std::map<char, std::vector<double> > map);
double getNeutronStrength(std::string name, std::map<char, double> map);
double getIonNeutronStrength(std::string name, std::map<std::string, double> map);
double getWaterNeutronStrength(std::string name, float deut, std::map<char, double> map);