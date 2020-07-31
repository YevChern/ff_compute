#include <map>
#include <vector>

using namespace std;

std::map<char, vector<double> > create_aff_map();
std::map<std::string, vector<double> > create_aff_ion_map();
std::map<char, double> create_nsld_map();
std::map<std::string, double> create_nsld_ion_map();
