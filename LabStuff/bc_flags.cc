#include "HalfEdgeMesh.h"
#include <math.h>

//#define PI 3.14159265 COMMENTED OUT FOR USING RANDOM NUMBER GENERATOR

double angle(double px, double py){
	double theta=0.0;
	if(py>=0.0){
		if(px>=0.0) theta = atan(py/px);
		if(px<0.0) theta = PI+atan(py/px);
	}
	if(py<0.0){
		if(px>=0.0) theta = 2.0*PI + atan(py/px);
		if(px<0.0) theta = PI+atan(py/px);
	}
	return theta;
}

int main(int argc, char** argv){
	
	//./bc_flags 10 square_d10_sub2.vtk square_d10_sub2_bc.vtk
	// reads distance input output
	
	if(argc != 8){
    	std::cout << "**Not enough arguments. Remember to put input file name, output file name, outer circle type (half or full),";
    	std::cout << " x value of inclusion center, contact angle, polar angle, and outer edge radius**" << std::endl;
    	return 1;
  	}

	HalfEdgeMesh mesh1(argv[1]);

	//mesh1.output();

	//double pack = atof(argv[3]);
	//int NN = 15;
	double pack = 1.26;
	int NN = 15; //NOT SURE WHAT THIS DOES, SEEMED IRRELEVANT
	
	double r = 2.27;
	double epsilon = 0.22;
	double radius = r*(1.+epsilon);
	double side = 2.0*radius*sin(36.0*PI/180.);
	double apothem = radius*cos(36.0*PI/180.);
	double axx = radius*(5.+sqrt(5.))/4.;
	double ayy = radius*sqrt( (9./8.)*(5.+sqrt(5.)) );

	double dx = pack*(radius*((1./8.) + (3./8.)*sqrt(5.)));
	double dy = pack*(0.5*ayy);
	double disp = pack*(radius+apothem);	
	double radx =0.5*(NN-1)*disp+11.0;
	double rady = 0.5*(NN-1)*dy+11.0;
	//double corner =0.5*dist;
	
	double Lx = 20., Ly=20.;
	//double cx = 0.5*dist, cy = 0.5*sqrt(3)*dist;
	
	std::vector<int> boundary_conditions;
	
	int nnodes = mesh1.get_n_nodes();
	double xx, yy;
	double eps = 0.01;
	//double rad = 38.0; 
	double sym=3.0;
	double ee = 0.2;
	//double rad2 = 2.0*2.0+ eps;
	for(int ii = 0; ii<nnodes; ii++){
		if(mesh1.checkBoundaryNode(ii)){
			xx = mesh1.getNodePosition(ii, 0);
			yy = mesh1.getNodePosition(ii, 1);
			//std::cout<<"Node: "<<ii<<" "<<xx<<" "<<yy<<"\n"; 
			
			//if(sqrt(xx*xx+yy*yy)<rad) boundary_conditions.push_back(7);
			if(xx<radx && xx>-radx && yy<rady && yy>-rady) boundary_conditions.push_back(7);

			else boundary_conditions.push_back(0);
		}
		
		else{
			boundary_conditions.push_back(0);
		}
		
	}

	mesh1.write_vtk_bc(argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), boundary_conditions);

	return 0;
}
