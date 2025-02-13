%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The programm "steady_travelling_waves" has been written by Matteo Antuono
and provides the analytical solution described in:

  Antuono (2022), "Steady periodic waves over a planar seabed: a global characterization",
                   J. Fluid Mech. (2022), vol. 947, A18, doi:10.1017/jfm.2022.643

NOTE: the OpenMP communication protocol has been used to speed-up the computation 
of series convolutions in the modal solution.


NOTE: the wave is assumed to propagate along DECREASING values of the x-axis.
To invert the motion, it is sufficient to change sign to the velocity field
(namely, to both horizontal and vertical components)



----------------------------------------------------------------------------------------------
The INPUT file is "ST_waves.dat"

In "ST_waves.dat" it is possible to select the following choices:
(o) assign the wave length (by setting 'flag_lp=0')
(o) assign the wave period (by setting 'flag_lp=1')


In the case 'flag_lp=0' (assignment of wavelength), it is possible to:
  (*) REMOVE the wave drift by setting 'flag_no_drift=.TRUE.'
  (*) ADD a superimposed wave current by assigning a non-null value to 'v_current'.
      A positive value indicates a current in the same direction of wave propagation.
NOTE: the above options ARE NOT AVAILABLE for the case 'flag_lp=1' (assignment of wave period)


For each of the above choices, it is possible to:
(o) save the solution at a fixed time (t=0) as a function of space (by setting 'flag_st=0')
(o) save the solution at a fixed position (x=0) as a function of time (by setting 'flag_st=1')





---------------------------------------------------------------------------------------------
The OUTUT files are:

(o) "filename_solution_FS.dat"  
     This file contains the solution along the free surface, that is:
     *) the free surface quote, 'eta'
     *) the velocity components at the free surface, 'u_fs' and 'w_fs'

(o) "filename_solution_3D.dat"
     This file contains the solution in the fluid bulk, that is:
     *) the horizontal component of the velocity field, 'u'
     *) the vertical component of the velocity field, 'w'
     *) the pressure field, 'p' 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
