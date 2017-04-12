// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                    (C) 2007 All Rights Reserved
//
//----------------------------------------------------------------------


#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include <map>
#include <cmath>
#include <string>
#include "ran.h"

  struct HalfEdge;
  struct Vertex;
  struct Face;

  //! A HalfEdge is the Left half of an edge, pointing CCW around a Face.
  struct HalfEdge {
    //! Default constructor
    HalfEdge() {
      opposite = next = prev = 0;
      face = 0; 
      vertex = 0;
      othervertex = 0; // Line added for SE
      id = -1;
    }

    //! Adjacent HalfEdge
    HalfEdge * opposite;
    //! HalfEdge Counter-Clockwise around Face on left
    HalfEdge * next;
    //! HalfEdge Clockwise around Face on left (stored for convenience)
    HalfEdge * prev;
    //! Face on left
    Face * face;
    //! vertex at the end of the HalfEdge
    Vertex * vertex;

    // Added for SE, need both vertices of an edge
    Vertex * othervertex;

    //! HalfEdge id
    int id;
  };
  
  struct Face {
    //! Default constructor
    Face() { halfEdges[0] = halfEdges[1] = halfEdges[2] = 0; id=-1; }

    //! HalfEdges looping Counter-Clockwise around Face
    HalfEdge * halfEdges[3];

    //! Face id
    int id;
  };

  struct Vertex {
    //! Default constructor
    Vertex() { id = -1; boundary = false; }

    //! vertex id
    int id;

    //! is this a boundary vertex?
    bool boundary;

    //! list of incident HalfEdges (pointing towards the vertex)
    std::vector< HalfEdge* > halfEdges;
  };


	typedef struct TriangleConnectivity {
		int tc[3];
	} TriangleConnectivity;
	typedef struct Position {
		double rr[3];
	} Position;
	
  //! HalfEdge data structure to store a triangle mesh
  /*!  This class is built to convert a connectivity array into a
    HalfEdge data structure.
  */
  struct HalfEdgeMesh {

    std::vector< HalfEdge* > halfEdges;
    std::vector< Face* > faces;
    std::vector< Vertex* > vertices;

    std::vector< std::vector< HalfEdge* > > boundaryLoops;

    //! Connectivity of a single triangle
    //typedef tvmet::Vector<int,3> TriangleConnectivity;
   // typedef int[3] TriangleConnectivity;

    //! Container of connectivities for a mesh of triangles
    typedef std::vector<TriangleConnectivity> ConnectivityContainer;
    typedef std::vector<Position> PositionContainer;

    // construct HalfEdge data structure from a connectivity array
    HalfEdgeMesh(const char* ifile);
    
    //
    void construct(const ConnectivityContainer & connectivities, int nVertices);
	PositionContainer _pos;
	ConnectivityContainer _connects;
	int n_nodes;
	int n_faces;
  //
  //
  // A bunch of stuff for write_vtk_bc
  std::ofstream out_to_vtk;
  //MAIN DEFINITIONS FOR CONTROLLING INCLUSION
  double prodist;
  double prorad;
  double proang;
  double contang;
  double prohght;
  // Defining values for the 'sphere' that a portion of will be used to represent the inclusion
  double sphrad;
  double sphcen2procen;
  double sphhght;
  double sphdist;
  double prothk;
  // xplane*(x-prodist) + zplane*(z-prohght) = 0
  double xplane;
  double zplane;
  // For outer edge of "membrane"
  double outerEdge;
  //
	//
  void VertexFormulation(std::set<int>& vertex_list, std::map<int, std::vector<int> >& vertex_constraints,
    double vv, bool halfCircle, Ran* rng, double noiseAmplitude);
  void EdgeFormulation(int i, double vv1, double vv2, std::map <std::pair<int,int>, int>& edge_list,
    std::map <int,std::vector<int> >& edge_constraints, std::map<int, std::vector<int> >& vertex_constraints);
  //
	void write_vtk_bc(const char* ofname, const std::string circleType, const double inclusionCenterX,
    const double contactAngle, const double polarAngle, const double outerEdgeRad, std::vector<int> bcs);
	void output();
	bool checkBoundaryNode(int idx) {return vertices[idx]->boundary;}
	int get_n_nodes(){return n_nodes;}
	double getNodePosition(int idx, int coord){return _pos[idx].rr[coord];}
	//

	
    // Delete structs
    virtual ~HalfEdgeMesh() {
      for(int f=0; f<faces.size(); f++) delete faces[f];
      for(int h=0; h<halfEdges.size(); h++) delete halfEdges[h];
      for(int v=0; v<vertices.size(); v++) delete vertices[v];
    }

  };

void HalfEdgeMesh::construct(const ConnectivityContainer & connectivities, 
			   int nVertices)
{

  std::cout << "Building HalfEdgeMesh..." << std::endl;

  // size containers
  faces.resize( connectivities.size() );
  for(int f=0; f<faces.size(); f++) faces[f] = new Face();
  halfEdges.resize( 3*connectivities.size() );
  for(int h=0; h<halfEdges.size(); h++) halfEdges[h] = new HalfEdge();
  vertices.resize(nVertices);
  for(int v=0; v<vertices.size(); v++) vertices[v] = new Vertex();

  // set vertex ids      
  for(int i=0; i<nVertices; i++) vertices[i]->id = i;

  // Link HalfEdges, Vertices, and Faces for each element; 
  // opposites remain undefined
  for(int f=0; f<connectivities.size(); f++) {
    Face * F = faces[f];
    F->id = f;
    for(int i=0; i<3; i++) {
      HalfEdge * HE = halfEdges[3*f+i];
      HE->id = 3*f+i;
      //int id=connectivities[f](i);
      int id=connectivities[f].tc[i];
      Vertex * V = vertices[id];

      // set vertex from connectivities
      HE->vertex = V;
      V->halfEdges.push_back( HE );

      // link HalfEdges to Face
      HE->face = F;
      F->halfEdges[i] = HE;

      // identify next HalfEdge (i+1 mod 3)
      int j = (i==2) ? 0 : i+1;
      HE->next = halfEdges[3*f+j];

      // identify prev HalfEdge (i-1 mod 3)
      j = (i==0) ? 2 : i-1;
      HE->prev = halfEdges[3*f+j];
    }
  }
  std::cout << "Linked " 
	    << halfEdges.size() << " half-edges, "
	    << faces.size() << " faces, and "
	    << vertices.size() << " vertices. "
	    << std::endl;

  std::cout << "Determining half-edge opposites." << std::endl;
  // For each  HalfEdge, find its opposite
  for(int h=0; h<halfEdges.size(); h++) {
    HalfEdge * H = halfEdges[h];
    // check if next is already set
    if( H->opposite != 0 ) continue;

    // if not, then look at all HalfEdges incident to my vertex and
    // check to see if they duplicate me (i.e., orientations are not
    // consistent) or if their next is my opposite

    for(int hv=0; hv<H->vertex->halfEdges.size(); hv++) {
      HalfEdge * Hv = H->vertex->halfEdges[hv];

      if( H == Hv ) continue;
    
      // Hv is a duplicate of H if both previous HalfEdges point to
      // the same vertex
      if( H->prev->vertex->id == Hv->prev->vertex->id ) {
	std::cout << "Error:  half-edges " << H->id << " and " << Hv->id
		  << " are duplicates.  " << std::endl
		  << "Adjacent faces " << H->face->id
		  << " and " << Hv->face->id 
		  << " must have opposite orientations." << std::endl;
	exit(0);
      }

      // Hv's next is my opposite if it points to the vertex of
      // my previous HalfEdge
      if( H->prev->vertex->id == Hv->next->vertex->id ) {
	H->opposite = Hv->next;
	break;
      }
    }

    if(H->opposite == 0 ) {
//       std::cout << "  Warning."
// 		<< "  Couldn't find HalfEdge opposite to HalfEdge " << H->id 
// 		<< ".  This HalfEdge is on the boundary."
// 		<< std::endl;

      H->vertex->boundary = true;

//       for(int hv=0; hv<H->vertex->halfEdges.size(); hv++) {
// 	HalfEdge * Hv = H->vertex->halfEdges[hv];
// 	std::cout << "h: "<< Hv->id << std::endl
// 		  << "v: "<< Hv->vertex->id << std::endl
// 		  << "n: "<< Hv->next->id << std::endl
// 		  << "nv: "<< Hv->next->vertex->id << std::endl
// 		  << "p: "<< Hv->prev->id << std::endl
// 		  << "pv: "<< Hv->prev->vertex->id << std::endl;
// 	std::cout << std::endl;
//       }
//      exit(0);
    }
	
  }

  // Find boundary loops. Identify a boundary, and walk around it starting
  // from one edge
  
  std::vector<int> edgeLoopLookup(halfEdges.size(), -1);

  // std::vector< std::vector< HalfEdge* > > boundaryLoops;

  for(int h=0; h<halfEdges.size(); h++) {

    //  find a boundary edge not previously added to a loop
    HalfEdge * H = halfEdges[h];
    if( H->opposite == 0 && edgeLoopLookup[H->id] == -1 ) {

      // found one, now walk around boundary, storing boundary edges
      // in order

      std::vector< HalfEdge * > loop;

      HalfEdge * Hstart=H;
      do {
	loop.push_back( H );
	edgeLoopLookup[H->id] = boundaryLoops.size();

	// find next boundary edge (CCW around the boundary) by
	// walking around H's vertex CW
	H = H->next;
	while ( H->opposite != 0 ) {
	  H = H->opposite->next;
	}
      } while ( H != Hstart );

      boundaryLoops.push_back( loop );
      
    }
  }

//   std::cout << "Identified " << boundaryLoops.size() << " boundary loops." 
// 	    << std::endl;
//   for(int L=0; L<boundaryLoops.size(); L++) {
//     std::cout << "Loop " << L << " has " 
// 	      << boundaryLoops[L].size() << " edges: ";
//     for(int e=0; e<boundaryLoops[L].size(); e++) {
//       std::cout << std::setw(10) << boundaryLoops[L][e]->id;
//     }
//     std::cout << std::endl;
//   }

  std::cout << "HalfEdgeMesh built." << std::endl;
	  return;
}

HalfEdgeMesh::HalfEdgeMesh(const char* ifile){
	const char* ifilename; ifilename = ifile;
	std::ifstream vtkstream;
	vtkstream.open(ifilename);
	if(vtkstream.is_open()) std::cout<<"reading the mesh file: "<<ifilename<<" ..."<<std::endl;
	else{ 
		std::cout<<"Unable to open mesh file\n";
		exit(0);
	}
	
	double xx, yy, zz;
	int i1,i2,i3;
	unsigned int idxc=0;
	
	
	
	std::string token; // Read the first string.  It should say "OFF"
	while( token != "POINTS" ) vtkstream >> token;
	vtkstream >> n_nodes; //std::cout<<"nodes: "<<n_nodes<<"\n";
	_pos.resize(n_nodes);
	
	
	// skip number type and read points
	vtkstream >> token;  //std::cout<<"token: "<<token<<"\n";  
	for(int i = 0; i < n_nodes; i++) {
		vtkstream >> xx >> yy >> zz; //std::cout<<"token: "<<xx<<" "<<yy<<" "<<zz<<"\n";  
		//MeshNode* nd = new MeshNode(idxc, xx, yy, zz);
		//idxc++;
		_pos[i].rr[0] = xx; _pos[i].rr[1] = yy; _pos[i].rr[2] = zz; 
		//_nodes[i]=nd;
	}
	
	// Read in triangle connectivities
	idxc=0;
	while( token != "POLYGONS" ) vtkstream >> token;
	vtkstream >> n_faces; //std::cout<<"token: "<<nfaces<<"\n";  	
	vtkstream >> token;
	_connects.resize(n_faces);
	for (int i = 0; i < n_faces; i++){
		vtkstream >> i1; //std::cout<<"token: "<<i1<<"\n"; 
		if(i1 != 3){
			std::cout << "Are you sure that this is a triangular mesh?\n" << std::endl;
			exit(0);
		}
		vtkstream >> i1 >> i2 >> i3;
		_connects[i].tc[0] =i1; _connects[i].tc[1] =i2; _connects[i].tc[2] =i3;
	}
  

	
	std::cout<<"Mesh has "<<n_nodes<<" nodes and "<<n_faces<<" faces.\n";
	
	construct(_connects, n_nodes);
}

void HalfEdgeMesh::output(){
	for(int i=0; i<n_nodes; i++){
		double xx = _pos[i].rr[0];
		double yy = _pos[i].rr[1];
		double zz = _pos[i].rr[2];
		std::cout<<"Node: "<<i<<" "<<xx<<" "<<yy<<" "<<zz<<"\n";
		std::cout<<"bc: "<<vertices[i]->boundary<<"\n";
	}
}

void HalfEdgeMesh::VertexFormulation(std::set<int>& vertex_list, std::map<int, std::vector<int> >& vertex_constraints,
  double vv, bool halfCircle, Ran* rng, double noiseAmplitude){

  if(vertex_list.find(vv) == vertex_list.end()){ // If vertex isn't in list yet, add to list and add to outfile
    vertex_list.insert(vv);
    double xx = _pos[vv].rr[0];
    double yy = _pos[vv].rr[1];
    double zz = _pos[vv].rr[2];
    
    std::string constraintInfo = "";
    if(std::pow((xx-prodist)/cos(proang),2)+std::pow(yy,2)-std::pow(prorad,2) <= 0.002){ // Check if node is within the (maybe oval) disk of radius r
    
      constraintInfo += "constraint 1";
        zz = sqrt(std::pow(sphrad,2)-std::pow(xx-sphdist,2)-std::pow(yy,2)) + sphhght;
      if(std::pow((xx-prodist)/cos(proang),2)+std::pow(yy,2)-std::pow(prorad,2) < -0.002){ // z value should be projected up to be on the sphere
//        zz = sqrt(std::pow(sphrad,2)-std::pow(xx-sphdist,2)-std::pow(yy,2)) + sphhght;
        std::vector<int> constraints;
        constraints.push_back(1);
        vertex_constraints.insert(std::make_pair<int,std::vector<int> >(vv,constraints));
        constraintInfo += " fixed"; //THIS MAY NEED TO BE ADDED, WE'LL SEE
      } else{
        constraintInfo += ", 2";
        std::vector<int> constraints;
        constraints.push_back(1); constraints.push_back(2);
        vertex_constraints.insert(std::make_pair<int,std::vector<int> >(vv,constraints));
        constraintInfo += " fixed"; //THIS MAY NEED TO BE ADDED, WE'LL SEE
      }
    } else{
      if(halfCircle){
        if(fabs(std::pow(xx,2) + std::pow(yy,2) - std::pow(outerEdge,2)) <= 0.002){ //Vertex is on outer edge of membrane, add to contraint 5
          constraintInfo += "constraint 5";
          std::vector<int> constraints;
          constraints.push_back(5);
          if(fabs(xx)<=0.002){ // Vertex is on x = 0 plane constraint
            constraintInfo += ", 6";
            constraints.push_back(6);
          }
          vertex_constraints.insert(std::make_pair<int,std::vector<int> >(vv,constraints));
        } else if(fabs(xx)<=0.002){ // Vertex is on x = 0 plane constraint
          constraintInfo += "constraint 6";
          std::vector<int> constraints;
          constraints.push_back(6);
          vertex_constraints.insert(std::make_pair<int,std::vector<int> >(vv,constraints));
        }
      } else{
        if(fabs(std::pow(xx,2) + std::pow(yy,2) - std::pow(outerEdge,2)) <= 0.002){ //Vertex is on outer edge of membrane, add to contraint 5
          constraintInfo += "constraint 5";
          std::vector<int> constraints;
          constraints.push_back(5);
          vertex_constraints.insert(std::make_pair<int,std::vector<int> >(vv,constraints));
        }
      }
    }

    // ADDING NOISE
    if(constraintInfo == ""){
      zz += noiseAmplitude*(.5-rng->doub());
    }

    out_to_vtk << vv+1 <<" "<< xx <<" "<< yy <<" "<< zz <<" "<<constraintInfo<<"\n"; // change for SE syntax
  }
}

void HalfEdgeMesh::EdgeFormulation(int i, double vv1, double vv2, std::map <std::pair<int,int>, int>& edge_list,
  std::map <int,std::vector<int> >& edge_constraints, std::map<int, std::vector<int> >& vertex_constraints){

  // Make pairs for either ordering of the vertices, don't want to add edges in reverse direction
    std::pair<int,int> e1_1 = std::make_pair<int,int>(vv1,vv2);
    std::pair<int,int> e1_2 = std::make_pair<int,int>(vv2,vv1);
    if(edge_list.find(e1_1) == edge_list.end()
     && edge_list.find(e1_2) == edge_list.end()){ // If neither direction of edge is in list yet, add to list and add to outfile
      int e1num = i;
      edge_list.insert(std::make_pair(e1_1,e1num));

      std::string eConstraintInfo = "";
      if(vertex_constraints.find(vv1)!=vertex_constraints.end() && vertex_constraints.find(vv2)!=vertex_constraints.end()){
        std::vector<int> vv1constraints = vertex_constraints.find(vv1)->second;
        std::vector<int> vv2constraints = vertex_constraints.find(vv2)->second;

        if(vv1constraints[0]==1 && vv2constraints[0]==1){
          eConstraintInfo += "constraint 1";
          if(vv1constraints[1]==2 && vv2constraints[1]==2){
            eConstraintInfo += ", 2";
          }
        } else if(vv1constraints[0]==5 && vv2constraints[0]==5){
          eConstraintInfo += "constraint 5";
          if((vv1constraints[0]==6 && vv2constraints[0]==6)){
            eConstraintInfo += ", 6";
          }
        } else if((vv1constraints[0]==6 && vv2constraints[0]==6)){
          eConstraintInfo += "constraint 6";
        }
      }

      out_to_vtk << e1num <<" "<< e1_1.first+1 <<" "<< e1_1.second+1 <<" "<<eConstraintInfo<<"\n";
      if(e1_1.first+1 == 1 || e1_1.second+1 == 1){
        out_to_vtk << "//WEHERE" << "\n";
      }
    }
}
  
void HalfEdgeMesh::write_vtk_bc(const char* ofname, const std::string circleType, const double inclusionCenterX, const double contactAngle,
  const double polarAngle, const double outerEdgeRad, std::vector<int> bcs){ // CHANGES FOR SE BEING MADE IN HERE MOSTLY
	
  out_to_vtk.open(ofname);
  bool halfCircle = true;
  std::string halfType = "half";
  for(int i = 0; i < circleType.size(); i++){
    if(tolower(circleType[i]) != halfType[i]){
      halfCircle = false;
    }
  }

  // Change header to work in SE

  // Setting variables for SE header

  //MAIN DEFINITIONS FOR CONTROLLING INCLUSION
  prodist = inclusionCenterX; // Distance of center of inclusion from edge of semicircle
  prorad = 1; // Radius of the circle that will represent the inclusion
  proang = (polarAngle*M_PI/180); // Angle between z-axis and 'polar' axis of inclusion
  contang = (contactAngle*M_PI/180); // Angle of contact between inclusion and "membrane"
  prohght = (prorad*sin(proang)); //height of inclusion center

  // Defining values for the 'sphere' that a portion of will be used to represent the inclusion
  sphrad = (prorad/cos((M_PI/2)-contang)); // Radius of sphere
  sphcen2procen = (sphrad*sin((M_PI/2)-contang)); //distance between center of sphere and center of 'inclusion cap'
  sphhght = (prohght-cos(contang)*cos(proang)*sphrad); //height of sphere center from x-y-plane
  sphdist = (prodist-sphcen2procen*sin(proang)); //distance of sphere's center from x-z plane
  prothk  = 0;//(sphrad-sphcen2procen); //thickness of inclusion (distance from center on plane of contact circle to top of inclusion cap)

  // xplane*(x-prodist) + zplane*(z-prohght) = 0
  xplane = (sin(proang)); // For equation of plane that holds the contact circle between inclusion and memberane
  zplane = (cos(proang)); // For equation of plane that holds the contact circle between inclusion and memberane

  // For outer edge of "membrane"
  outerEdge = outerEdgeRad;


  out_to_vtk<<"//Surface Evolver file generated from vtk file from GMSH\n"
    <<"keep_originals //makes life easier\n\n";

  out_to_vtk<<"//Curvature energy definition\n"
    <<"quantity starsq energy method star_sq_mean_curvature\n\n"

    <<"//MAIN DEFINITIONS FOR CONTROLLING INCLUSION\n"
    <<"#define prodist ("<<prodist<<") //distance of center of protein from edge of semicircle\n"
    <<"#define prorad  ("<<prorad<<") //radius of the circle that will represent the protein\n"
    <<"#define proang  ("<<proang<<") //angle between z-axis and 'polar' axis of protein\n"
    <<"#define contang ("<<contang<<") //angle of contact between protein and membrane\n"
    <<"#define prohght ("<<prohght<<") //height of protein center\n\n"

    <<"//defining values for the 'sphere' that a portion of will be used to represent the protein\n"
    <<"#define sphrad  ("<<sphrad<<") //radius of sphere\n"
    <<"#define sphcen2procen ("<<sphcen2procen<<") //distance between center of sphere and center of 'protein cap'\n"
    <<"#define sphhght ("<<sphhght<<") //height of sphere center from x-y-plane\n"
    <<"#define sphdist ("<<sphdist<<") //distance of sphere's center from x-z plane\n"
    <<"#define prothk  ("<<prothk<<") //thickness of protein (distance from center on plane of contact circle to top of protein cap)\n\n"

    <<"//defining values/equations for the plane that holds the contact circle between protein and membrane\n"
    <<"//xplane*(x-prodist) + zplane*(z-prohght) = 0\n"
    <<"#define xplane ("<<xplane<<")\n"
    <<"#define zplane ("<<zplane<<")\n\n"

    // Trying something new for constraints to make polar angle adjustment easier:
    <<"//Making polar angle adjustment easier (hopefully)\n"
    <<"parameter incl_cen = prodist\n"
    <<"parameter incl_rad = prorad\n"
    <<"parameter incl_ang = proang\n"
    <<"parameter cont_ang = contang\n"
    <<"parameter incl_hght = prohght\n\n"

    <<"parameter sph_rad = sphrad\n"
    <<"parameter sphcen_to_inclcen = sphcen2procen\n"
    <<"parameter sph_hght = sphhght\n"
    <<"parameter sph_dist = sphdist\n"
    <<"parameter incl_thk = prothk\n"

    <<"parameter x_plane = xplane\n"
    <<"parameter z_plane = zplane\n\n"


    <<"//Constraints for Protein\n\n"

    <<"constraint 1 //'sphere' for which the cap represents the protein\n"
    //<<"formula: 0=0\n\n"
    <<"formula: (x-sph_dist)^2 + (y)^2 + (z-sph_hght)^2 = (sph_rad)^2\n\n"

    <<"constraint 2 //plane that the base of protein is constrainted to\n"
    //<<"formula: 0 = 0\n\n"
    <<"formula: x_plane*(x-incl_cen) + z_plane*(z-incl_hght) = 0\n\n"

    <<"//any other constraints\n"
    <<"//constraints for outer massive membrane\n"
    <<"constraint 5\n"
    <<"formula: x^2 + y^2 = "<<outerEdge<<"^2\n\n";

    if(halfCircle){
      out_to_vtk<<"constraint 6\n"
      <<"formula: x = 0 \n\n";
    }


  // Make a set to make sure no duplicate vertices are being made
  std::set<int> vertex_list;
  std::map<int, std::vector<int> > vertex_constraints;

  // Initialization for using random numbers for noise
  Ran* rng = new Ran(1);
  double noiseAmplitude = 0;

  out_to_vtk<<"\nVertices\n"; // SE syntax
	for(int i=0; i<n_faces; i++){
    double v1 = _connects[i].tc[0];
    HalfEdgeMesh::VertexFormulation(vertex_list, vertex_constraints, v1, halfCircle, rng, noiseAmplitude);

    double v2 = _connects[i].tc[1];
    HalfEdgeMesh::VertexFormulation(vertex_list, vertex_constraints, v2, halfCircle, rng, noiseAmplitude);
    
    double v3 = _connects[i].tc[2];
    HalfEdgeMesh::VertexFormulation(vertex_list, vertex_constraints, v3, halfCircle, rng, noiseAmplitude);
	}

  //Make map to make sure I'm not making duplicate edges, and to store edge number for any set of vertices,
  //and keep track of which face it was added for (each edge has 2 faces, will need to know which it was added for
  //in order to create faces with properly facing normals)
  std::map <std::pair<int,int>, int> edge_list;
  std::map <int,std::vector<int> > edge_constraints;

  out_to_vtk<<"\nEdges\n"; // SE syntax
  for(int i=0; i<n_faces; i++){ // Figure out how to add all the edges correctly based off of vertices
    double v1 = _connects[i].tc[0];
    double v2 = _connects[i].tc[1];
    double v3 = _connects[i].tc[2];

    HalfEdgeMesh::EdgeFormulation(3*i+1, v1, v2, edge_list, edge_constraints, vertex_constraints);
    HalfEdgeMesh::EdgeFormulation(3*i+2, v2, v3, edge_list, edge_constraints, vertex_constraints);
    HalfEdgeMesh::EdgeFormulation(3*i+3, v3, v1, edge_list, edge_constraints, vertex_constraints);
  }

	out_to_vtk<<"\nFaces\n"; // SE syntax
	for(int i=0; i<n_faces; i++){ // Add all the faces, ensuring correct edge orientation and order
		// Get the vertices for finding the edge numbers, as faces are built using edge numbers in SE
    double v1 = _connects[i].tc[0];
    double v2 = _connects[i].tc[1];
    double v3 = _connects[i].tc[2];

    // Make pairs of vertices and search edge map to find the corresponding edge number for either direction
    std::pair<int,int> e1_1 = std::make_pair(v1,v2);
    std::pair<int,int> e1_2 = std::make_pair(v2,v1);
    std::pair<int,int> e2_1 = std::make_pair(v2,v3);
    std::pair<int,int> e2_2 = std::make_pair(v3,v2);
    std::pair<int,int> e3_1 = std::make_pair(v3,v1);
    std::pair<int,int> e3_2 = std::make_pair(v1,v3);
    int e1num;
    int e2num;
    int e3num;
    if(edge_list.find(e1_1) == edge_list.end()){ // Then other pair is in the list and goes the "wrong" way
      e1num = -1*(edge_list.find(e1_2)->second);
    } else{ // First pair is in the list
      e1num = edge_list.find(e1_1)->second;
    }

    if(edge_list.find(e2_1) == edge_list.end()){ // Then other pair is in the list and goes the "wrong" way
      e2num = -1*(edge_list.find(e2_2)->second);   
    } else{
      e2num = edge_list.find(e2_1)->second;
    }

    if(edge_list.find(e3_1) == edge_list.end()){
      e3num = -1*(edge_list.find(e3_2)->second);
    } else{
      e3num = edge_list.find(e3_1)->second;
    }

		out_to_vtk << i+1 <<" "<< e1num <<" "<< e2num <<" "<< e3num <</*" "<<constraintinfo<<*/"\n";
	}

  out_to_vtk << "\n\n"
    <<"read\n\n"

    <<"r :::= { 'r'; set vertex starsq where not on_constraint 1; set vertex starsq where on_constraint 2; 'e' }\n\n"

    <<"u :::= { 'u'; set vertex starsq where not on_constraint 1; set vertex starsq where on_constraint 2; 'e' }\n\n"

    <<"E :::= {u; u; hessian_seek; hessian_seek; hessian_seek; u; u;}\n\n"

    <<"SUPPRESS_WARNING 1792\n\n"

    <<"set vertex starsq where not on_constraint 1\n"
    <<"set vertex starsq where on_constraint 2\n\n"

    <<"set face color red where on_constraint 1\n"
    <<"set edge color 13 where on_constraint 2\n"
    <<"set facet tension 0\n";
}

