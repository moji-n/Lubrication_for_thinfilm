# Lubrication_thinfilm

This is a prgoram that solves the non-linear evolution of spreading a thin film driven by curvature, gravity and disjoining pressure.

Input:
The input parapmeters can be found in the /input/ directory,
  . disjoining_parametres.txt  contains the parameters for calculating the disjoining pressure 
  . parameters.txt             contains the general parameters of the model such as geometry, time-stepping, inital conditions, choosing different mechanicns ...
  . timesteps.txt              contains the times instances which the output file is created
  
Output:

The output of the code describes the film thickness versus location at seperate file for each time instance located in the directory /output/.

Files:
"Lubrication_main.m"      : This is the main M file , where it reads the input parameters, initilize the model, make the mesh, and solve the linear set from the discretized differential equation.
"Disjoin_der.m"           : This function calculates the derivate of different components of disjoining pressure.
"interpolate_Mobility.m"  : This function interpolate the mobility term at the middle of the nodes.
"interpolate_refine.m"    : This function interpolate the film thickness at the created nodes after a change in the mesh.
"populate_matrixes.m"     : This function populate the coefficient matrix of the discretized differential equation.
"grid_gen_refin.m"        : This function generate the mesh and refines it close to the contact point.

"h_min_plot.m"            : A post-processing code, that plots the variation of film thickness at the contact point versus time from the file "output_general.txt"
"plot_lub_results.m"      : A post-processing code, that plot the spreading radius and film profile versus time.
