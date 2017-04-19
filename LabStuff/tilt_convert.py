import subprocess as sproc
import sys
import os

tilt_size = 0.05
rec_fn = "temp_blah@#iI2#"
evolv_cmd = "../../Evolver/evolver.exe"
devnull = open(os.devnull, 'w')


## Main function of program ##
## This program will run polar_change with small angle changes and then in-surface-evolver minimizations ##
## to reach the requested polar angle change ##
def polar_change(in_fn, out_fn, e_fn, polar_angle_change):

	in_fn = "fromGMSH/" + in_fn
	out_fn = "mesh/" + out_fn
	e_fn = "energy/" + e_fn

	# Optional line to compile the c++ code before calling, given that it doesn't auto-compile like python does
	#sproc.call("g++ -O3 polar_change.cc -o polar_change", shell=True)

	if(float(polar_angle_change) > tilt_size):
		loop_runs = int((float(polar_angle_change)//tilt_size))
		final_change = float(polar_angle_change)%tilt_size

		sproc.call("./polar_change " + in_fn + " " + rec_fn + str(0) + ".txt" + " " + e_fn + " " + str(tilt_size), shell=True)
		sproc.call(evolv_cmd + " " + rec_fn + str(0) + ".txt", shell=True, stdout = devnull)

		for i in range(loop_runs):
			sproc.call("./polar_change " + rec_fn + str(i) + ".txt.dmp" + " " + rec_fn + str(i+1) + ".txt"
				+ " " + e_fn + " " + str(tilt_size), shell=True)
			sproc.call("rm " + rec_fn + str(i) + ".txt", shell=True)
			sproc.call("rm " + rec_fn + str(i) + ".txt.dmp", shell=True)
			sproc.call(evolv_cmd + " " + rec_fn + str(i+1) + ".txt", shell=True, stdout = devnull)

		sproc.call("./polar_change " + rec_fn + str(loop_runs) + ".txt.dmp" + " " + out_fn + " " + e_fn
			+ " " + str(final_change), shell=True)
		sproc.call("rm " + rec_fn + str(loop_runs) + ".txt", shell=True)
		sproc.call("rm " + rec_fn + str(loop_runs) + ".txt.dmp", shell=True)

	else:
		sproc.call("./polar_change " + in_fn + " " + out_fn + " " + e_fn + " " + polar_angle_change, shell=True)

	# Open the file in surface evolver
	# sproc.call("evolver " + SEfn, shell=True)




## Check that program is run correctly, let user know how to run if not ##
if(len(sys.argv) != 5):
	print "Usage: \n\t"
	print sys.argv[0]+" [inputfile] [outputfile] [energyfile] [polar angle change]"
	exit()
else:
	polar_change(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])