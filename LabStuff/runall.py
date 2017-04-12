import subprocess as sproc
import sys

## Main function of program ##
def runall(GMSHfn, SEfn, circle_type, inclusion_center_x, contact_angle, polar_angle, outer_radius):
	# Run GMSH with the .geo file, creating a 2-d mesh and outputting as a vtk file unlikely to overwrite another vtk
	sproc.call("gmsh " + GMSHfn + " -2 -o zAz_temp_192038_1.vtk", shell=True)

	# Run vtk_poly.py with the .vtk file from GMSH preparing for bc_flags, outputting as a vtk file unlikely to overwrite another vtk
	sproc.call("python vtk_poly.py zAz_temp_192038_1.vtk zAz_temp_192038_2.vtk", shell=True)

	# Optional line to compile the c++ code before calling, given that it doesn't auto-compile like python does
	sproc.call("g++ -O3 bc_flags.cc -o bc_flags", shell=True)

	# Run bc_flags with the .vtk file from vtk_poly and the necessary parameters, write out to user specified file name
	SEfn = "fromGMSH/" + SEfn
	sproc.call("./bc_flags zAz_temp_192038_2.vtk " + SEfn + " " + circle_type + " " + inclusion_center_x + " "
		+ contact_angle + " " + polar_angle + " "+ outer_radius, shell=True)

	# Clean-up
	sproc.call("rm zAz_temp_192038_1.vtk", shell=True)
	sproc.call("rm zAz_temp_192038_2.vtk", shell=True)

	# Open the file in surface evolver
	# sproc.call("evolver " + SEfn, shell=True)




## Check that program is run correctly, let user know how to run if not ##
if(len(sys.argv) != 8):
	print "Usage: \n\t"
	print sys.argv[0]+" [inputfile] [outputfile] [outer circle type (half or full)] [x value of inclusion center] [contact angle of inclusion] [polar angle of inclusion] [outer radius]"
	exit()
else:
	runall(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])