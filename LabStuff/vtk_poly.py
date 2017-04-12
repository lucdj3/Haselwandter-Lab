#! /usr/bin/python
import sys

def vtk_poly(ifn, ofn):
	
	ifile = open(ifn,'r')
	ofile = open(ofn, 'w')
	
	lines = ifile.readlines()

	nNodes = int(lines[4].split()[1])

	
	ofile.write("# vtk DataFile Version 3.0 \nvtk output \nASCII \nDATASET POLYDATA\n")
	ofile.write("POINTS " + str(nNodes) + " double\n")
	for ln in lines[5:5+nNodes]:
		ofile.write(ln)
	
	nFacesAll = int(lines[6+nNodes].split()[1])
	triangles = []
	for ln in lines[7+nNodes: 7+nNodes+nFacesAll]:
		sp = ln.split()
		if(sp[0]=='3'):
			triangles.append(ln)
	
	nFaces = len(triangles)
	ofile.write("POLYGONS "+str(nFaces)+" "+ str(nFaces*4)+"\n")
	for tri in triangles:
		ofile.write(tri)
	

if(len(sys.argv)!=3):
	print "Usage: \n\t"+sys.argv[0]+" [inputfile] [outputfile]"
	exit()
else:
	vtk_poly(sys.argv[1], sys.argv[2])
