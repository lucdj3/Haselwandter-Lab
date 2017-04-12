#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>

void polar_angle_tilt(const char* inFName, const char* outFName, const char* energyFName, const double polarAngleChange);

//MAIN DEFINITIONS FOR CONTROLLING INCLUSION
double incl_cen;
double incl_rad;
double incl_ang;
double cont_ang;
double incl_hght;
// Defining values for the 'sphere' that a portion of will be used to represent the inclusion
double sph_rad;
double sphcen_to_inclcen;
double sph_hght;
double sph_dist;
double incl_thk;
// xplane*(x-prodist) + zplane*(z-prohght) = 0
double x_plane;
double z_plane;

int main(int argc, char** argv){
	if(argc != 5){
    	std::cout << "**Not enough arguments. Remember to put input file, mesh output file,"
	    	<<" energy output file, and polar angle change (in degrees)**" << std::endl;
    	return 1;
  	}

  	polar_angle_tilt(argv[1], argv[2], argv[3], atof(argv[4]));
}

void polar_angle_tilt(const char* inFName, const char* outFName, const char* energyFName, const double polarAngleChange){

	// Try to open the file
	std::ifstream iFile;
	iFile.open(inFName);
	if(iFile.is_open()) std::cout<<"reading the mesh file: "<<inFName<<" ..."<<std::endl;
	else{ 
		std::cout<<"Unable to open mesh file: "<<inFName<<"...\n";
		exit(0);
	}

	// Make a temporary file to write everything to
	std::ofstream oFile;
	oFile.open(outFName);

	// Open Energy file for appending
	std::ofstream eFile;
	eFile.open(energyFName, std::fstream::out | std::fstream::app);

	// Find parameters to fix
	std::string token;
	bool dumpQuitBool = false;
	while(std::getline(iFile, token)){
		std::vector<std::string> allTokens;
		std::istringstream ss(token);
		for(std::string ssToVec; ss >> ssToVec;){
		    allTokens.push_back(ssToVec);
		}
		if(allTokens.size() == 0){
			oFile<<"\n";
		} else if(allTokens[0] == "PARAMETER"){ // Fix all the parameters
			if(allTokens[1] == "incl_cen"){
				incl_cen = atof(allTokens[3].c_str());
				oFile<<"PARAMETER incl_cen = "<<incl_cen<<"\n";
			} else if(allTokens[1] == "incl_rad"){
				incl_rad = atof(allTokens[3].c_str());
			    oFile<<"PARAMETER incl_rad = "<<incl_rad<<"\n";
			} else if(allTokens[1] == "incl_ang"){
				incl_ang = atof(allTokens[3].c_str()) + (polarAngleChange)*M_PI/180;
			    oFile<<"PARAMETER incl_ang = "<<incl_ang<<"\n";
			} else if(allTokens[1] == "cont_ang"){
				cont_ang = atof(allTokens[3].c_str());
				oFile<<"PARAMETER cont_ang = "<<cont_ang<<"\n";
			} else if(allTokens[1] == "incl_hght"){
				incl_hght = incl_rad*sin(incl_ang);
			    oFile<<"PARAMETER incl_hght = "<<incl_hght<<"\n";
			} else if(allTokens[1] == "sph_rad"){
				sph_rad = incl_rad/cos((M_PI/2)-cont_ang);
			    oFile<<"PARAMETER sph_rad = "<<sph_rad<<"\n";
			} else if(allTokens[1] == "sphcen_to_inclcen"){
				sphcen_to_inclcen = sph_rad*sin((M_PI/2)-cont_ang);
			    oFile<<"PARAMETER sphcen_to_inclcen = "<<sphcen_to_inclcen<<"\n";
			} else if(allTokens[1] == "sph_hght"){
				sph_hght = incl_hght-cos(cont_ang)*cos(incl_ang)*sph_rad;
			    oFile<<"PARAMETER sph_hght = "<<sph_hght<<"\n";
			} else if(allTokens[1] == "sph_dist"){
				sph_dist = incl_cen-sphcen_to_inclcen*sin(incl_ang);
			    oFile<<"PARAMETER sph_dist = "<<sph_dist<<"\n";
			} else if(allTokens[1] == "incl_thk"){
				incl_thk = 0;//(sphrad-sphcen2procen); //thickness of inclusion (distance from center on plane of contact circle to top of inclusion cap)
			    oFile<<"PARAMETER incl_thk = "<<incl_thk<<"\n";
			} else if(allTokens[1] == "x_plane"){
				x_plane = sin(incl_ang);
			    oFile<<"PARAMETER x_plane = "<<x_plane<<"\n";
			} else if(allTokens[1] == "z_plane"){
				z_plane = cos(incl_ang);
			    oFile<<"PARAMETER z_plane = "<<z_plane<<"\n";
			}
		} else if(allTokens[0].find('_') != -1){
			if(allTokens[0].find('_', allTokens[0].find('_')+1) != -1){
				if(allTokens[0].substr(allTokens[0].find('_', allTokens[0].find('_')+1)+1) == "predicted"){
					continue;
				}else{
					oFile<<token<<"\n";
				}
			}else{
				if(allTokens[0].substr(allTokens[0].find('_')+1) == "predicted"){
					continue;
				}else{
					oFile<<token<<"\n";
				}
			}
		} else if(allTokens.size() > 3){
			if(allTokens[1] == "Total" && allTokens[2] == "energy:"){

				eFile<<allTokens[3]<<"\n";
			} else{
				oFile<<token<<"\n";
			}
		} else{
			if(token == "dump") dumpQuitBool = true;
			oFile<<token<<"\n";
		}
	}
	if(!dumpQuitBool) oFile<<"E\nE\ndump\nq";
	iFile.close();
	oFile.close();
}