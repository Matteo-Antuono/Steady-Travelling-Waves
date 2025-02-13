!!######################################################################
!-----------------------------------------------------------------------
      MODULE type_set

      type wave_mode
        real(8) :: freq 
        real(8) :: kwave
        real(8) :: amp 
        real(8) :: phi 
      end type      
        
      type wave_solution
        real(8),allocatable,dimension(:)    :: xs    !! xs(phi)=x(psi_s,phi) // coordinate x along the free surface 
        real(8),allocatable,dimension(:)    :: eta   !! eta(phi)         // free surface elevation
        real(8),allocatable,dimension(:)    :: dxs   !! d(xs)/d(phi)     // derivative of xs(phi)
        real(8),allocatable,dimension(:)    :: deta  !! d(xs)/d(phi)     // derivative of xs(phi)
        real(8),allocatable,dimension(:)    :: u_fs  !! horizontal velocity at the FS
        real(8),allocatable,dimension(:)    :: w_fs  !! vertical velocity at the FS
        real(8),allocatable,dimension(:,:)  :: u     !! horizontal velocity
        real(8),allocatable,dimension(:,:)  :: w     !! vertical velocity
        real(8),allocatable,dimension(:,:)  :: p     !! pressure
      end type
  
      type wave_solution_modes
        real(8),allocatable,dimension(:) :: mode
      end type

      END MODULE type_set
!-----------------------------------------------------------------------      
!#######################################################################
!-
!!######################################################################
!=======================================================================
!-----------------------------------------------------------------------
      MODULE global_variables
           
      USE  type_set  

      character(30)          :: filename,filename_app      
      real(8),dimension(15)  :: par_sol  !! parameters of the solution   

      !! solution & modes
      TYPE(wave_solution)       :: wave_sol
      TYPE(wave_solution_modes) :: wave_modes
           
      
      END MODULE global_variables
!-----------------------------------------------------------------------
!!!#####################################################################
!--
!!######################################################################
!=======================================================================
!-----------------------------------------------------------------------
      MODULE global_parameters
           
      real(8),parameter    :: pi=ACOS(-1.d0)
      real(8),parameter    :: grav=9.81d0     
      real(8),parameter    :: toll_fs=1.E-9
      integer(4),parameter :: iter_max=1000
      
      END MODULE global_parameters
!-----------------------------------------------------------------------
!!!#####################################################################
!-
!!######################################################################
!!----------------------------------------------------------------------  
      PROGRAM steady_travelling_waves
 
      ! REFERENCE for the analytic solution:
      !   Antuono M., Steady periodic waves over a planar seabed: a global characterization, 
      !               J. Fluid Mech. (2022), vol. 947, A18, doi:10.1017/jfm.2022.643 
 
      ! COMPILING OPTIONS: 
      ! COMPILING (DBG-OPENMP) : gfortran -g -fdefault-real-8 -fdefault-double-8 -O0 -Wall -Wextra -pedantic -fbounds-check -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -fopenmp -o ST_waves steady_travelling_waves.f90   
      ! COMPILING (OPT-OPENMP) : gfortran -O3 -fdefault-real-8 -fdefault-double-8 -fopenmp -o ST_waves steady_travelling_waves.f90  
      
      USE global_parameters
      USE global_variables
      
      IMPLICIT NONE
      
      integer(4)    :: flag_lp,flag_st
      real(8)       :: h0,a0,T_per,Length
      integer(4)    :: nsample,nsample_z 
      real(8)       :: v_current
      logical       :: flag_no_drift
      !----------
      real(8)           :: eps_bar,mu_bar  !! barred "epsilon" and "mu" as defined in Antuono (2022)
      real(8)           :: omega_bar       !! dimensionless circular frequency (with shallow water scaling)
      real(8)           :: eps,mu          !! parameters of the scaling proposed in Antuono (2022)
      real(8)           :: H_bot           !! seabed quote in the (phi,psi)-plane (i.e. hodograph space)
      !---------
      real(8)           :: z0,kwave,omega            !! dimensional parameters
      real(8)           :: U_drift,cel_earth,cel_an  !! dimensional velocities
      real(8)           :: err_omega                 
      integer(4)        :: cont 
      !---------
      real(8),external  :: disp_rel,inverse_disp_rel
      !----------------------------------------------------------------------------------------------------
      
      
      write(*,*)
      write(*,*) 
      write(*,*) '==============================================================================='
      write(*,*) 
      write(*,*) '        STEADY PERIODIC WAVES over a planar seabed - Antuono (2022) '
      write(*,*) '               written by Matteo Antuono, 28th November 2024        ' 
      write(*,*) 
      write(*,*) '-------------------------------------------------------------------------------'
      
      
      !! READING data
      open(777,file='ST_waves.dat',status='unknown')
        read(777,*) 
        read(777,*)         
        read(777,*) 
 
        read(777,*) h0                
        read(777,*) a0                
        read(777,*) flag_lp        !! '0' = assignement of wave length   OR   '1' = assignement of wave period 
        if (flag_lp.EQ.0) THEN
          read(777,*) Length   
        else
          read(777,*) T_per      
        endif      
        read(777,*) flag_st        !! '0' = solution in space   OR   '1' = solution in time     
        read(777,*) nsample            
        read(777,*) v_current            
        read(777,*) flag_no_drift  !! 'TRUE' if you want to eliminate the Drift velocity from the solution 
        read(777,*) nsample_z 
        read(777,*) filename              
      close(777)
      !! nonlinearity parameter using shallow water scales
      eps_bar = a0/h0



      !! Writing LOG file
      filename_app = trim(ADJUSTL(trim(ADJUSTL(filename))//'_LOG.dat'))
      open(777,file=filename_app,status='unknown')
        write(777,*) 'LOG file'
        write(777,*)
        write(777,*) 'INPUT DATA'
        write(777,*) 'h0: ', h0                  
        write(777,*) 'a0: ', a0        
        write(777,*) 'flag_lp: ', flag_lp     
        if (flag_lp.EQ.0) THEN
          write(777,*) 'Length: ', Length   
        else
          write(777,*) 'Period: ',T_per      
        endif           
        write(777,*) 'flag_st: ', flag_st
        write(777,*) 'nsample: ', nsample
        write(777,*) 'v_current: ', v_current
        write(777,*) 'flag_no_drift: ', flag_no_drift
        write(777,*) 'nsample_z: ', nsample_z
        write(777,*) 'filename: ', filename
        
        write(777,*)
        write(777,*)
        if (flag_lp.EQ.0) THEN
          write(777,*)  'The wavelength is assigned => compute the wave period/circular frequency'
          write(777,*)
        else
          write(777,*)  'The wave period is assigned => compute the wavelength'  
          write(777,*)  'starting ITERATIVE PROCEDURE'
          write(777,*)
        endif         
 
        
                        
      IF (flag_lp.EQ.0) THEN    !! The wavelength is assigned => compute the wave period/circular frequency

        mu_bar = 2.d0*pi*h0/Length  !! barred dispersive parameter as in Antuono (2022)
        mu     = dtanh(mu_bar)      !! dispersive parameter as in Antuono (2022)
        H_bot  = mu_bar/mu          !! seabed quote in the (phi,psi)-plane
        eps    = H_bot*eps_bar      !! steepness/nonlinearity parameter as in Antuono (2022)
        z0     = h0/H_bot           !! vertical length scale

        !! Solution (computing modes and other wave parameters)    
        call Antuono_waves_init(mu,eps,h_bot,h0,v_current,nsample,flag_no_drift)
         
        cel_an  =       sqrt(grav*z0) * par_sol(7)  !! dimensional celerity (analytical solution)
        U_drift = eps * sqrt(grav*z0) * par_sol(8)  !! dimensional drift velocity 
        IF (flag_no_drift) THEN                                    
           cel_earth = cel_an + v_current - abs(U_drift)   !! celerity in the earth reference frame (WITHOUT drift) 
        ELSE   
           cel_earth = cel_an + v_current                  !! celerity in the earth reference frame (WITH drift)
        ENDIF                                                      
        kwave = mu/z0
        omega = cel_earth*kwave

        !! DIMENSIONAL parameters
        write(*,*) 'dimensional parameters'
        write(*,*) 'Wave Period: T = ', (2.0*pi/omega), ' sec'
        write(*,*) 'Wave Length: L = ', Length, ' m'        
        write(*,*) '-------------------------------------------------------------------------------'
        write(*,*) 'dimensionless parameters'
        write(*,*) 'eps     ', eps,      '  mu     ', mu
        write(*,*) 'eps_bar ', eps_bar,  '  mu_bar ', mu_bar
        write(*,*) 'wave steepness:  sigma ', eps_bar*mu_bar
        write(*,*) '-------------------------------------------------------------------------------'  
        
        write(777,*) 'dimensional parameters'
        write(777,*) 'Wave Period: T = ', (2.0*pi/omega), ' sec'
        write(777,*) 'Wave Length: L = ', Length, ' m'        
        write(777,*) '-------------------------------------------------------------------------------'
        write(777,*) 'dimensionless parameters'
        write(777,*) 'eps     ', eps,      '  mu     ', mu
        write(777,*) 'eps_bar ', eps_bar,  '  mu_bar ', mu_bar
        write(777,*) 'wave steepness:  sigma ', eps_bar*mu_bar        
        write(777,*) '-------------------------------------------------------------------------------'         
            
      
      ELSE  !! wave period is assigned => compute the wavelength
      
        !! OSS: IN THIS CASE (i.e. wave period assignement) the inversion of the dispersion relation  
        !! in the presence of wave current (1) and/or subtraction of the drift velocity (2) is not available. 
        !! Below the motivations for these issues:
        !!
        !! (1) the subroutine "inverse_disp_rel" should be modified to include the wave current 
        !!     after proper upper and lower bound are provided for the bisection method.
        !! (2) the Drift velocity is a function of mu_bar but, unfortunately, there is no 
        !!     regression formula for this dependence at the moment. After such a dependence is provided,
        !!     the subroutine "inverse_disp_rel" should be modified accordingly
         
 
        if (flag_no_drift) then
          write(*,*) 'ERR: the inversion of the disperive relation is not available with the option "NO DRIFT"'
          write(*,*) '========================================================================================'      
          write(*,*)  
          stop
        endif  
        if (ABS(v_current).GT.1.E-12) then
          write(*,*) 'ERR: the inversion of the disperive relation is not available when a current is present'
          write(*,*) '=======================================================================================' 
          write(*,*)    
          stop
        endif
        omega_bar = 2.d0*pi/T_per
        omega_bar = omega_bar*sqrt(h0/grav)              !! dimensionless circular frequency using the shallow water scales
        mu_bar    = inverse_disp_rel(omega_bar,eps_bar)  !! wavelenght estimated through the approximate dispersive relation 
      
        mu     = dtanh(mu_bar)      !! dispersive parameter as in Antuono (2022)
        H_bot  = mu_bar/mu          !! seabed quote in the (phi,psi)-plane
        eps    = H_bot*eps_bar      !! steepness/nonlinearity parameter as in Antuono (2022)
        z0     = h0/H_bot           !! vertical length scale        

        !! Solution (computing modes and other wave parameters)  
        call Antuono_waves_init(mu,eps,h_bot,h0,v_current,nsample,flag_no_drift)
         
        cel_an    = sqrt(grav*z0) * par_sol(7)  !! dimensional celerity (analytical solution)
        cel_earth = cel_an
        kwave  = mu/z0
        omega  = cel_earth*kwave    !! dimensional circular frequency using scales defined in Antuono (2022)

     
        !!-----------------------------------------------------------------
        !! Iterative procedure to improve the estimation of the wavelength     
              
        err_omega = (2.d0*pi/T_per)/omega - 1.d0 !! relative error with respect to the assigned wave period/circular frequency
             
        cont=1
        DO WHILE ((cont.LE.10).AND.(ABS(err_omega).GT.1.E-6)) 
          
               
          write(777,*) 'ASSIGNED wave period: T = ', T_per, ' sec'
          write(777,*) 'COMPUTED wave period: T = ', 2.d0*pi/omega, ' sec'
          write(777,*) '>>> relative error in circular frequency estimation:  ', err_omega
          write(777,*) '-------------------------------------------------------------------------------'
          write(777,*) 'Wave Length: L = ', 2.d0*pi/kwave, ' m'        
          write(777,*) '-------------------------------------------------------------------------------'
          write(777,*) 'dimensionless parameters'
          write(777,*) 'eps     ', eps,      '  mu     ', mu
          write(777,*) 'eps_bar ', eps_bar,  '  mu_bar ', mu_bar
          write(777,*) '-------------------------------------------------------------------------------'            
               
        
          !! deallocating the previous solution
          deallocate(wave_modes%mode)
          
          !! correction of the estimated "mu_bar". NOTE: d(omega)/d(mu) > 0
          mu_bar = mu_bar/( omega/(2.d0*pi/T_per) )  
          
          mu     = dtanh(mu_bar)      !! dispersive parameter as in Antuono (2022)
          H_bot  = mu_bar/mu          !! seabed quote in the (phi,psi)-plane
          eps    = H_bot*eps_bar      !! steepness/nonlinearity parameter as in Antuono (2022)
          z0     = h0/H_bot           !! vertical length scale        
  
          !! Solution (computing modes and other wave parameters)     
          call Antuono_waves_init(mu,eps,h_bot,h0,v_current,nsample,flag_no_drift)
           
          cel_an    = sqrt(grav*z0) * par_sol(7)  !! dimensional celerity (analytical solution)
          cel_earth = cel_an
          kwave     = mu/z0
          omega     = cel_earth*kwave    !! dimensional circular frequency using scales defined in Antuono (2022)
  
          err_omega = (2.d0*pi/T_per)/omega - 1.d0  !! updated relative error    
        
          cont = cont+1
          
        ENDDO
 
        !! Parameters
        write(*,*) 'dimensional parameters'         
        write(*,*) 'ASSIGNED wave period: T = ', T_per, ' sec'
        write(*,*) 'COMPUTED wave period: T = ', 2.d0*pi/omega, ' sec'
        write(*,*) '>>> relative error in circular frequency estimation:  ', err_omega, '  cont = ', cont
        write(*,*) '-------------------------------------------------------------------------------'
        write(*,*) 'Wave Length: L = ', 2.d0*pi/kwave, ' m'        
        write(*,*) '-------------------------------------------------------------------------------'
        write(*,*) 'dimensionless parameters'
        write(*,*) 'eps     ', eps,      '  mu     ', mu
        write(*,*) 'eps_bar ', eps_bar,  '  mu_bar ', mu_bar
        write(*,*) '-------------------------------------------------------------------------------'   
        
        write(777,*) 'dimensional parameters'         
        write(777,*) 'ASSIGNED wave period: T = ', T_per, ' sec'
        write(777,*) 'COMPUTED wave period: T = ', 2.d0*pi/omega, ' sec'
        write(777,*) '>>> relative error in circular frequency estimation:  ', err_omega, '  cont = ', cont
        write(777,*) '-------------------------------------------------------------------------------'
        write(777,*) 'Wave Length: L = ', 2.d0*pi/kwave, ' m'        
        write(777,*) '-------------------------------------------------------------------------------'
        write(777,*) 'dimensionless parameters'
        write(777,*) 'eps     ', eps,      '  mu     ', mu
        write(777,*) 'eps_bar ', eps_bar,  '  mu_bar ', mu_bar
        write(777,*) '-------------------------------------------------------------------------------'        
            
      ENDIF
      !! closing LOG file
      close(777)

     
      !! saving 
      par_sol(13) = 2.d0*pi/omega  !! computed wave period
      par_sol(14) = 2.d0*pi/kwave  !! computed wave length
      par_sol(15) = omega/kwave    !! computed wave celerity (possibly accounting for drift and current corrections)


      !! BUILDING the SOLUTIONS at the Free-Surface and inside the fluid bulk
      !! then SAVING the solutions
      call solution_xz(flag_st,nsample_z)

 
      STOP
      END PROGRAM 
!-----------------------------------------------------------------------
!=======================================================================
!!-
!!=======================================================================
!!-----------------------------------------------------------------------
      SUBROUTINE Antuono_waves_init(mu,eps,h_bot,h0,v_current,nsample,flag_no_drift)
      
      ! Iterative scheme for the solution of travelling permanent waves over arbitrary depths

      USE global_parameters
      USE global_variables
      USE omp_lib

      IMPLICIT NONE
      
      real(8),intent(in)                :: mu,eps,H_bot,h0
      integer(4)                        :: n,n_modes,n_iter
      integer(4),parameter              :: n_iter_max=5000
      real(8),parameter                 :: tol_iter=1.E-9  
      real(8)                           :: Ud,chi_mod
      real(8)                           :: theta_mode_app,Flim 
      real(8),allocatable,dimension(:)  :: theta_mode,theta_mode_prev
      real(8),allocatable,dimension(:)  :: C_theta_mode,Cm1_theta_mode,D_theta_mode
      real(8),allocatable,dimension(:)  :: term1_mode,term2_mode
      real(8),allocatable,dimension(:)  :: S_term_mode,CS_term_mode
      real(8),allocatable,dimension(:)  :: Z_term_mode,CZ_term_mode      
      real(8),dimension(n_iter_max)     :: chi,beta,psi_s,eta_0,nu
      real(8)                           :: Cn_theta,Cn_chi,Cn_beta,Cn_eta0
      real(8),dimension(n_iter_max)     :: Cn_iter,qn_iter
      real(8)                           :: conv_iter,par,res_norm
      real(8)                           :: num_app1,den_app1,theta_max
      real(8)                           :: start_loc,finish_loc,Froude    
                                 
      integer(4),intent(in)  :: nsample        !! sample point for the solution
      real(8)                :: z0  

      logical,intent(in)     :: flag_no_drift  !! TRUE if we want to eliminate the Drift velocity
      real(8),intent(in)     :: v_current      !! DIMENSIONAL velocity of an eventual current (NOTE: a positive sign indicates a current in the same directin of the wave) 

 
      
 
      !! Qualitative estimate of the number of modes of the solution     
      n_modes = nint( (eps/0.4d0) * (500.d0)/mu ) 

      
      !! Initialization 
      ALLOCATE(theta_mode(-n_modes:n_modes),theta_mode_prev(-n_modes:n_modes))
      ALLOCATE(S_term_mode(-n_modes:n_modes),Z_term_mode(-n_modes:n_modes))   
      ALLOCATE(C_theta_mode(-n_modes:n_modes),Cm1_theta_mode(-n_modes:n_modes))
      ALLOCATE(D_theta_mode(-n_modes:n_modes),term1_mode(-n_modes:n_modes),term2_mode(-n_modes:n_modes) )
      ALLOCATE(CS_term_mode(-n_modes:n_modes),CZ_term_mode(-n_modes:n_modes))
      
      theta_mode(:)      = 0.d0     !! modes of "theta"
      theta_mode_prev(:) = 0.d0     !! modes of "theta" at the previous iteration
      D_theta_mode(:)    = 0.d0     !! modes of "dtheta/dlam"   
      C_theta_mode(:)    = 0.d0     !! modes of "C[theta]"
      Cm1_theta_mode(:)  = 0.d0     !! modes of "C^(-1)[theta]"
      S_term_mode(:)     = 0.d0     !! modes of "S = Cm1^2 + dtheta/dlam^2"
      Z_term_mode(:)     = 0.d0     !! modes of "Z = 2*theta*Cm1 + eps*theta*S"
      CS_term_mode(:)    = 0.d0     !! modes of "C[S]"
      CZ_term_mode(:)    = 0.d0     !! modes of "C[Z]"

      chi(:)     = 0.d0  !! chi - preallocation
      beta(:)    = 0.d0  !! beta - preallocation
      Psi_s(:)   = 0.d0  !! Psi_s - preallocation
      eta_0(:)   = 0.d0  !! eta_0 - preallocation
      Cn_iter(:) = 1.d0  !! C_n (convergence) - preallocation
      qn_iter(:) = 1.d0  !! q_n (rate of convergence) - preallocation
      nu(:)      = 0.d0  !! parameter of singularity, chi-2*eps*theta_max
      
      Cn_theta = 0.d0  !! C_n (convergence) for "theta"
      Cn_chi   = 0.d0  !! C_n (convergence) for "chi"
      Cn_beta  = 0.d0  !! C_n (convergence) for "beta"
      Cn_eta0  = 0.d0  !! C_n (convergence) for "eta0"
 
      
      !!------------------------------------------------------------------
      !! INIZIALIZATION     
      !! initial assignment:  theta=cos(mu*Phi)
      theta_mode_prev(-1) = 0.5d0
      theta_mode_prev( 1) = 0.5d0 
      !!------------------------------------------------------------------      
      


      !! ITERATIVE SCHEME
      Psi_s(1)  = H_bot     !! assuming eta_0=0
      par       = Psi_s(1)  !! inner parameter for the operators C and C^(-1)      
      conv_iter = 1000.d0
      n_iter    = 1
      
      start_loc = OMP_GET_WTIME() 
      DO WHILE ((n_iter.LT.n_iter_max).AND.(conv_iter.GT.tol_iter))
      

        !! NB: f = F0 + sum_{n>=1} Fn cos(ny)   equiv   f = sum_{n in Z} f_n exp(iny)   <==>    F_0 = f_0   &   F_n = 2 * f_n
        !! 
        !! NB: there is no need to compute the zero-mode since "theta" has a null mean, i.e. theta_mode_prev(0)=0 
        DO n=1,n_modes
        
          !! theta
          theta_mode_app = theta_mode_prev(n)

          !! C[theta]
          C_theta_mode( n)   = ( dtanh(dble(n)*mu*par)/(dble(n)*mu) ) * theta_mode_app  
          C_theta_mode(-n)   = C_theta_mode(n)
          
          !! C^(-1)[theta]
          Cm1_theta_mode( n) = ( (dble(n)*mu)/dtanh(dble(n)*mu*par) ) * theta_mode_app     
          Cm1_theta_mode(-n) = Cm1_theta_mode(n) 
          
          !! dtheta/dlamda
          D_theta_mode( n)   = -dble(n)*mu*theta_mode_app  !! NOTE: this is a sine series but, later, the derivative is squared (note: the immaginary number 'i' is here omitted, see later)  <<<
          D_theta_mode(-n)   = -D_theta_mode(n)            !! anti-symmetric 
          
        ENDDO        
        theta_max   = SUM(theta_mode_prev)     !sum over ALL the modes (OK)
        
      
        !! Evaluating "S_term" 
        call Discrete_Convolution(n_modes,Cm1_theta_mode,Cm1_theta_mode,term1_mode)
        call Discrete_Convolution(n_modes,D_theta_mode,D_theta_mode,term2_mode)     
        S_term_mode(:) = term1_mode(:) - term2_mode(:)     !! NOTE: "minus sign" => here the immaginary number 'i' from "D_theta_mode" is taken into account again  <<<

        !! Evaluating "Z_term" 
        call Discrete_Convolution(n_modes,theta_mode_prev,Cm1_theta_mode,term1_mode)  !! theta * C^(-1)[theta]
        call Discrete_Convolution(n_modes,theta_mode_prev,S_term_mode,term2_mode)     !! theta * S              
        Z_term_mode(:) = 2.d0*term1_mode(:) + eps*term2_mode(:) 
                   
        !! Evaluating "eta_0"
        eta_0(n_iter) = -eps * term1_mode(0)        
        

        DO n=1,n_modes  
             
          !! C[S]
          CS_term_mode( n) = ( dtanh(dble(n)*mu*par)/(dble(n)*mu) ) * S_term_mode(n)
          CS_term_mode(-n) = CS_term_mode(n)
          
          !! C[Z]
          CZ_term_mode( n) = ( dtanh(dble(n)*mu*par)/(dble(n)*mu) ) * Z_term_mode(n) 
          CZ_term_mode(-n) = CZ_term_mode(n)
        
        ENDDO
        CS_term_mode(0) = par * S_term_mode(0) 
        CZ_term_mode(0) = par * Z_term_mode(0)
        
          
        !! Evaluating "nu" (OK)
        !! NB: if we consider  f = sum_{n>=0} Fn cos(ny)   =>    [[ f ]] = 2 * sum_{n odd > 0} F_n  
        !! NB: f = F0 + sum_{n>=1} Fn cos(ny)   equiv   f = sum_{n in Z} f_n exp(iny)   <==>    F_0 = f_0   &   F_n = 2 * f_n
        !! ==> the computed summations are 1/4 of their expected value           
        num_app1 = 0.d0
        den_app1 = 1.d0
        DO n=1,n_modes,2
          num_app1 = num_app1 + ( C_theta_mode(n) + eps * CZ_term_mode(n) )
          den_app1 = den_app1 + eps * CS_term_mode(n)
        ENDDO
        chi(n_iter) = 2.d0*(num_app1/den_app1)
        nu(n_iter)  = chi(n_iter) - 2.d0*eps*theta_max



        !! WELL-POSEDNESS 
        !! Computed using "theta_prev" (correct choice)
        IF (nu(n_iter).LE.0.d0) THEN
          write(*,*) '////  SINGULAR solution  ////'
          write(*,*) 'nu:        ', nu(n_iter)
          write(*,*) 'chi:       ', chi(n_iter)          
          write(*,*) 'iteration: ', n_iter
          write(*,*) 'theta_max: ', theta_max
          write(*,*) '-----------------------------------------------------'
          STOP
        ENDIF        

            
        
        !! Evaluating "beta"
        beta(n_iter) = Z_term_mode(0) - 0.5d0 * chi(n_iter) * S_term_mode(0)       
   
        !! Evaluating "theta"
        theta_mode(:) = ( C_theta_mode(:) + eps * CZ_term_mode(:) )/chi(n_iter) - 0.5d0 * eps * CS_term_mode(:)      
        theta_mode(0) = 0.d0     !!NB: this is equivalent to subtract "eps*beta/chi" from the former equation
           
        !! Evaluating "Psi_s"    
        Psi_s(n_iter) = H_bot + eps * eta_0(n_iter)   
        par           = Psi_s(n_iter)  !! to be used in the next iteration

        
        !!-------------------------------------------------------------
        !! CONVERGENCE 
        
        !! L2 norm
        den_app1  = 0.d0
        conv_iter = 0.d0
        DO n=1,n_modes
          conv_iter = conv_iter + ( theta_mode(n) - theta_mode_prev(n) )**2  
          den_app1  = den_app1  + ( theta_mode(n) )**2                        
        ENDDO      
        Cn_theta = DSQRT(conv_iter/den_app1)

        IF (n_iter.GT.2) THEN
          Cn_chi   = DABS(chi(n_iter)  -  chi(n_iter-1))/DABS(chi(n_iter))
          Cn_beta  = DABS(beta(n_iter) - beta(n_iter-1))/DABS(beta(n_iter))
          Cn_eta0  = DABS(eta_0(n_iter)-eta_0(n_iter-1))/DABS(eta_0(n_iter))
        ENDIF
        conv_iter       = MAXVAL( (/ Cn_theta, Cn_chi, Cn_beta, Cn_eta0 /) )
        Cn_iter(n_iter) = conv_iter
        IF (n_iter.GT.2) THEN
          qn_iter(n_iter) = DLOG( Cn_iter(n_iter)/Cn_iter(n_iter-1) )/ &
                            DLOG( Cn_iter(n_iter-1)/Cn_iter(n_iter-2) )
        ENDIF 
        !!--------------------------------------------------------------       

  
        !! subsequent iteration
        DO n=1,n_modes
          theta_mode_app      = Flim( theta_mode(n) )
          theta_mode( n)      = theta_mode_app
          theta_mode(-n)      = theta_mode_app
          theta_mode_prev( n) = theta_mode_app   
          theta_mode_prev(-n) = theta_mode_app 
        ENDDO
        
        n_iter = n_iter+1
        
 
      ENDDO
      finish_loc = OMP_GET_WTIME() 
 
 
      !! Froude number
      Froude = DSQRT( chi(n_iter-1)-2.d0*beta(n_iter-1)*eps**2 )   
      !! Drift velocity
      Ud = Froude*eta_0(n_iter-1)/H_bot  
      !! vertical length scale of the wave propagation phenomenon
      z0 = h0/H_bot 
      
      !! Additional variable to be used in the inflow subroutine
      !! chi_mod = F^2 + 2 eps^2 B = F^2 + 2 eps^2 ( beta + eta0/eps ) 
      !!         = F^2 + 2 eps^2 beta + 2 eps eta0 = chi + 2 eps eta0
      chi_mod = chi(n_iter-1) + 2.d0 * eps * eta_0(n_iter-1) 
 
      par_sol(1)  = dble(nsample)
      par_sol(2)  = mu
      par_sol(3)  = eps 
      par_sol(4)  = H_bot 
      par_sol(5)  = z0     
      par_sol(6)  = chi_mod 
      par_sol(7)  = Froude  
      par_sol(8)  = Ud  !! NOTE: Ud<0 because eta0<0
      par_sol(9)  = par !! Psi_s   !! to be used for 3D
      par_sol(10) = dble(n_modes)  !! to be used for 3D
      par_sol(11) = v_current      !! eventual current (dimensional)
      IF (flag_no_drift) THEN
        par_sol(12) = 1.d0  !! to eliminate the Drift velocity
      ELSE 
        par_sol(12) = 0.d0  
      ENDIF
      
      

      write(777,*) '====================================================='
      write(777,*) 'total computational time:  ', finish_loc-start_loc     
      write(777,*) 'number of modes:       ', n_modes
      write(777,*) 'number of iterations:  ', n_iter-1
      write(777,*) 'convergence:    ', conv_iter
      write(777,*) 'Well Posedness: ', nu(n_iter-1)   


      !! Energy condition
      !! Percentage of Energy in the last tenth of modes
      res_norm = 0.d0
      den_app1 = 0.d0
      DO n=1,n_modes
        den_app1 = den_app1 + theta_mode_prev(n)**2 
        IF (n.LT.nint(dble(n_modes*9)/10.d0)) cycle      
        res_norm = res_norm + theta_mode_prev(n)**2                      
      ENDDO 
      res_norm = res_norm/den_app1  
      
      write(777,*) 'Average of "eta^2" over one period = ', (2.d0*den_app1 + eta_0(n_iter-1)**2)  
      write(777,*) '      energy condition, res_energy = ', res_norm
      write(777,*) '-----------------------------------------------------'
      IF (flag_no_drift) THEN
        write(777,*) '==> solution WITHOUT DRIFT' 
      ELSE
        write(777,*) '==> solution WITH DRIFT'
      ENDIF  
      write(777,*) '====================================================='
         

 
      !! Saving the modes for eta(Phi)
      ALLOCATE(wave_modes%mode(0:n_modes))
      wave_modes%mode(0) = eta_0(n_iter-1)  
      DO n=1,n_modes 
        wave_modes%mode(n) = 2.d0*theta_mode(n)     
      ENDDO


      !! DE-allocations
      DEALLOCATE(theta_mode,theta_mode_prev)
      DEALLOCATE(S_term_mode,Z_term_mode)   
      DEALLOCATE(C_theta_mode,Cm1_theta_mode)
      DEALLOCATE(D_theta_mode,term1_mode,term2_mode)
      DEALLOCATE(CS_term_mode,CZ_term_mode)     
     
    
      RETURN
      END SUBROUTINE Antuono_waves_init
!!----------------------------------------------------------------------       
!!######################################################################
!
!!======================================================================
!!----------------------------------------------------------------------  
      SUBROUTINE Discrete_Convolution(n_modes,Fn,Gn,Conv) 
      
      ! Discrete Convolution for cosine-series
      
      USE omp_lib 
    
      IMPLICIT NONE   
      
      integer(4),intent(in)                            :: n_modes     
      integer(4)                                       :: i,j
      real(8),dimension(-n_modes:n_modes),intent(in)   :: Fn,Gn
      real(8),dimension(-n_modes:n_modes),intent(out)  :: Conv
      real(8)                                          :: appo,Flim      


      !! Using option "dynamic"
      !$omp parallel  &
      !$omp default(shared)
      
      !$omp do                  &
      !$omp private(i,j,appo)   &    
      !$omp schedule(dynamic,10)        
      DO i=0,n_modes   
        appo = 0.d0
        !! optimized version
        DO j=i-n_modes,n_modes               
          !! simmetric form
          appo = appo + ( Fn(j)*Gn(i-j) + Fn(i-j)*Gn(j) )
        ENDDO  
        appo     = Flim(0.5d0*appo)   
        Conv(i)  = appo
        Conv(-i) = appo  !!symmetry: imposing cosine series      
      ENDDO
      !$omp end do
      !$omp end parallel           
      
         
      RETURN
      ENDSUBROUTINE Discrete_Convolution
!!----------------------------------------------------------------------       
!!======================================================================
!-
!!=======================================================================
!!-----------------------------------------------------------------------
      function Flim(arg)
      !!----------------------------------
      !! if ABS(arg) < bound => Flim = 0
      !! otherwise  Flim = arg
      !!---------------------------------        
      IMPLICIT NONE             

      real(8)            :: Flim,arg
      real(8),parameter  :: bound = 1.E-15  !! NB: it must be smaller than "tol_iter"
      !-----------------------------------
        
      Flim = arg
      IF (DABS(arg).LT.bound) Flim = 0.d0   
       
      RETURN
      ENDfunction Flim 
!!-----------------------------------------------------------------------
!!=======================================================================
!-
!!=======================================================================
!!-----------------------------------------------------------------------
      SUBROUTINE solution_xz(flag_st,nsample_z)

      ! Building the free surface solution:
      ! eta: free surface quote
      ! (u,w): velocity components along x and z

      USE global_parameters
      USE global_variables
      USE omp_lib

      IMPLICIT NONE
 
      integer(4),intent(in) :: flag_st,nsample_z     
      integer(4)            :: nsample,n_modes,i,n,j
      real(8)               :: mu,eps,z0,chi_mod,Froude,H_bot,U_drift,par
  
      logical     :: flag_no_drift   !! TRUE if we want to eliminate the Drift velocity
      real(8)     :: v_current       !! velocity of an eventual current (NOTE: a positive sign indicates a current in the same directin of the wave)       
      real(8)     :: T_in,L_in,C_in  !! period, length and celerity

      real(8),allocatable,dimension(:) :: xi,phis !! composite variable xi = x + c*t // variable phi as function of xi at the free surface 
      real(8)                          :: dxi,nmu,dzi,xx,zz,eta_ref,pp,Bern,scal 
      real(8),dimension(2)             :: vel_sol
      character(1)                     :: char_st
      

      nsample   = nint( par_sol(1) ) 
      mu        = par_sol(2) 
      eps       = par_sol(3) 
      H_bot     = par_sol(4) 
      z0        = par_sol(5) 
      chi_mod   = par_sol(6) 
      Froude    = par_sol(7) 
      U_drift   = par_sol(8)   !! NOTE: U_drift<0
      par       = par_sol(9)   !! Psi_s
      n_modes   = nint(par_sol(10))  
      v_current = par_sol(11)  !!NOTE: "v_current" is a DIMENSIONAL variable
      IF (nint(par_sol(12)).EQ.1) THEN
        flag_no_drift=.TRUE.
      ELSE
        flag_no_drift=.FALSE.
      ENDIF
      T_in = par_sol(13) !wave period
      L_in = par_sol(14) !wave Length
      C_in = par_sol(15) !wave celerity
      
      Bern = 0.5d0*(chi_mod - Froude**2)/eps**2  !! Bernoulli constant (dimensionless)
      
      
      !! allocating the spatial variable 
      allocate(xi(nsample),phis(nsample))
      
      !! sample points for the solution 
      !! NB: for simplicity a uniform grid is adopted
      dxi = L_in/(nsample)
      do i=1,nsample
        xi(i)   = -0.5d0*L_in + dxi*(i-1)
      enddo
      xi   = xi/z0 !! DIMENSIOLESS variables
      phis = xi/z0 !! initialization  phi_s = xi
                  
      !! Inverting the relation xs=xs(phi)
      !! See the equation (D5) in Antuono (2022)
      call inversion_phi_xi_fs(n_modes,nsample,xi,phis,eps,mu,par)
      
      xi = xi*z0  !!DIMENSIONAL variables
      

      !! allocating the solution (sampling at the Free Surface)
      allocate(wave_sol%xs(nsample))
      allocate(wave_sol%eta(nsample))
      allocate(wave_sol%dxs(nsample))
      allocate(wave_sol%deta(nsample))  
      allocate(wave_sol%u_fs(nsample))
      allocate(wave_sol%w_fs(nsample))
  
      wave_sol%xs   = xi
      wave_sol%eta  = 0.d0
      wave_sol%dxs  = 0.d0  !! d(xs)/d(phi)
      wave_sol%deta = 0.d0  !! d(eta)/d(phi)
      wave_sol%u_fs = 0.d0 
      wave_sol%w_fs = 0.d0 
      
      do i=1,nsample
      
        DO n=1,n_modes 
          nmu = dble(n)*mu
          wave_sol%eta(i)  = wave_sol%eta(i)  + wave_modes%mode(n) * dcos(nmu*phis(i))          
          wave_sol%dxs(i)  = wave_sol%dxs(i)  + wave_modes%mode(n)/dtanh(nmu*par) * nmu * dcos(nmu*phis(i))
          wave_sol%deta(i) = wave_sol%deta(i) - wave_modes%mode(n) * nmu * dsin(nmu*phis(i))
        ENDDO         
        wave_sol%eta(i) = wave_modes%mode(0) + wave_sol%eta(i)  ! eta(x)
        wave_sol%dxs(i) = 1.d0 + eps * wave_sol%dxs(i)          ! d(xs)/d(Phi) evaluated at Phi(x)  
        
        wave_sol%u_fs(i) = ( chi_mod - 2.d0*eps*wave_sol%eta(i) ) * wave_sol%dxs(i)  / Froude**2
        wave_sol%u_fs(i) = wave_sol%u_fs(i) - 1.d0  !! earth frame of reference
        
        wave_sol%w_fs(i) = eps * ( chi_mod - 2.d0*eps*wave_sol%eta(i) ) * wave_sol%deta(i) / Froude**2

      enddo
      
      !! dimensional variables
      wave_sol%eta  = eps * z0 * wave_sol%eta
      IF (flag_no_drift) THEN                                    
         wave_sol%u_fs = eps * sqrt(grav*z0) * wave_sol%u_fs + v_current - abs(U_drift)   !! WITHOUT drift
      ELSE   
         wave_sol%u_fs = eps * sqrt(grav*z0) * wave_sol%u_fs + v_current                  !! WITH drift
      ENDIF 
      wave_sol%w_fs = eps * sqrt(grav*z0) * wave_sol%w_fs
      
      
      !! SCALING for SAVING the solution
      !! xi = x + c*t
      IF (flag_st.EQ.0) THEN  !! solution in space (at t=0)
        scal = 1.d0
        char_st = 'x'
      ELSE       !! solution in time (at x=0)
        scal = c_in   
        char_st = 't'
      ENDIF    
      
      filename_app = trim(ADJUSTL(trim(ADJUSTL(filename))//'_solution_FS.dat'))
      open(unit=999,file=filename_app,status='unknown',form='formatted')
      write(999,*) ' VARIABLES= "'//char_st//'", "eta", "u_fs", "w_fs" '
      do i=1,nsample
          write(999,*) wave_sol%xs(i)/scal,wave_sol%eta(i),wave_sol%u_fs(i),wave_sol%w_fs(i)
      enddo
      close(999)
 
      
      
      !!-------------------------------------------------------
      !! allocating the solution (sampling in the fluid bulk)
      allocate(wave_sol%u(nsample,nsample_z))
      allocate(wave_sol%w(nsample,nsample_z)) 
      allocate(wave_sol%p(nsample,nsample_z))     
      wave_sol%u = 0.d0 
      wave_sol%w = 0.d0
      wave_sol%p = 0.d0
      
      
      !! Sample points for the solution in the bulk
      !! NOTE: the z coordinate goes from the seabed bottom to the free surface    
    
      !! Using option "dynamic"
      !$omp parallel  &
      !$omp default(shared)
      
      !$omp do                                         &
      !$omp private(i,j,xx,zz,eta_ref,dzi,pp,vel_sol)  &    
      !$omp schedule(dynamic,2)                   
      do i=1,nsample
        
        eta_ref = wave_sol%eta(i)/z0   !! DIMENSIONLESS FREE SURFACE quote
        
        xx  = xi(i)/z0   !! DIMENSIONLESS POSITION along x
        
        dzi = ( H_bot + eta_ref )/(nsample_z)  !! DIMENSIONLESS STEP along z
        
        do j=1,nsample_z
          
          zz = -H_bot + dzi*(j-1) !! DIMENSIONLESS POSITION along z
        
          call inversion_xz(xx,zz,par,mu,eps,n_modes,vel_sol,eta_ref)  
        
          wave_sol%u(i,j) = vel_sol(1) * eps * sqrt(grav*z0)
          wave_sol%w(i,j) = vel_sol(2) * eps * sqrt(grav*z0)
 
          !! dimensionless pressure
          pp = Bern - zz/(eps**2) - ( vel_sol(1) + 0.5d0*(vel_sol(1)**2 + vel_sol(2)**2) )*(Froude**2)/(eps**2)        
          wave_sol%p(i,j) = (eps**2) * grav * z0 * pp  !! DIMENSIONAL PRESSURE 
         
        enddo
      enddo
      !$omp end do
      !$omp end parallel    
     
    
      !! SAVING the solution in the bulk
      filename_app = trim(ADJUSTL(trim(ADJUSTL(filename))//'_solution_3D.dat'))
      open(unit=888,file=filename_app,status='unknown',form='formatted')
      write(888,*) ' VARIABLES= "'//char_st//'", "z", "u", "w", "p" '
      write(888,*) 'ZONE  T="',0.d0,'" I=',nsample_z+1,' J=',nsample             
      do i=1,nsample  
        eta_ref = wave_sol%eta(i)/z0
        xx      = xi(i)      
        dzi     = ( H_bot + eta_ref )/(nsample_z)  
        do j=1,nsample_z        
          zz  = ( -H_bot + dzi*(j-1) )*z0
          write(888,*) xx/scal, zz, wave_sol%u(i,j), wave_sol%w(i,j), wave_sol%p(i,j)
        enddo    
        write(888,*) xx/scal, wave_sol%eta(i), wave_sol%u_fs(i), wave_sol%w_fs(i), 0.d0
      enddo
      close(888)    
      
      
      write(*,*) ' ==> SOLUTION SAVED! '
      write(*,*) '==============================================================================='
      write(*,*)
      write(*,*)
      

      RETURN 
      ENDSUBROUTINE solution_xz
!!-----------------------------------------------------------------------
!!=======================================================================
!-
!!=======================================================================
!!-----------------------------------------------------------------------
      SUBROUTINE inversion_phi_xi_fs(n_modes,nsample,xi,phis,eps,mu,par)

      !! Inversion of the equation xs(phi) through a fixed point algorithm
      !! using the formula (D5) of Antuono (2022)
      !!
      !! OUTPUT: we obtain  phi = phi(xs)
      !! to be used for the computation of the solutions along the free surface


      USE global_parameters
      USE global_variables
      USE omp_lib

      IMPLICIT NONE
      
      integer(4),intent(in)                    :: n_modes,nsample
      integer(4)                               :: i,n,iter
      real(8),intent(in)                       :: eps,mu,par
      real(8),dimension(nsample),intent(in)    :: xi
      real(8),dimension(nsample),intent(inout) :: phis
      real(8)                                  :: phis_prev
      real(8)                                  :: err_iter,err_iter_tot
      real(8)                                  :: nmu
      
      err_iter_tot = 0.d0
      
      DO i=1,nsample
        
        iter = 1
        
        err_iter = 1.E+3
        
        DO WHILE ((iter.LT.iter_max).AND.(err_iter.GT.toll_fs))

          phis_prev = phis(i)         
          phis(i)   = xi(i) 
          
          DO n=1,n_modes 
            nmu = dble(n)*mu
            phis(i) = phis(i) - eps * wave_modes%mode(n)/dtanh(nmu*par) * dsin(nmu*phis_prev)
          ENDDO      
        
          err_iter = abs(phis(i)-phis_prev)
                   
          iter = iter+1

        ENDDO  
        
        err_iter_tot = err_iter_tot + err_iter
              
      ENDDO


      IF (err_iter_tot.GT.toll_fs*iter_max) write(*,*) 'Problem in inversion_phi_xi_fs'


      RETURN 
      ENDSUBROUTINE inversion_phi_xi_fs
!-----------------------------------------------------------------------
!=======================================================================
!- 
!=======================================================================    
!-----------------------------------------------------------------------
      SUBROUTINE inversion_xz(x0,z0,par,mu,eps,n_modes,vel_sol,eta_ref)
                                                  
      !! Finds the potentials Phi and Psi as functions of (x0,z0)
      !! Then, it computes the solutions
                             
      USE global_variables
      USE omp_lib
            
      IMPLICIT NONE       
      
      integer(4)                       :: iter,n
      integer(4),parameter             :: iter_max=1000
      integer(4),intent(in)            :: n_modes
      real(8),intent(in)               :: x0,z0,par,mu,eps
      real(8),intent(in)               :: eta_ref      
      real(8),dimension(2),intent(out) :: vel_sol
      real(8),parameter                :: tol_exit=1.E-9
      real(8)                          :: err_iter
      real(8)                          :: phi0,psi0
      real(8)                          :: phi0_prev,psi0_prev
      real(8)                          :: dxdphi,dzdphi,nmu
      real(8)                          :: cosh_o_cosh,sinh_o_sinh
      real(8)                          :: Flim
 
      !! initial guess
      phi0_prev = x0
      psi0_prev = min( z0 - eta_ref + par, par)      
      
      iter     = 1
      err_iter = 1.0
 
      DO WHILE ( (err_iter.GT.tol_exit).AND.(iter.LT.iter_max) )
      
        iter=iter+1
        
        phi0 = x0
        psi0 = z0 - eta_ref + par                
 
        DO n=1,n_modes 
        
          nmu  = dble(n)*mu
 
          phi0 = phi0 - eps * wave_modes%mode(n) * &
                 ( cosh_o_cosh(nmu*psi0_prev,nmu*par)/dtanh(nmu*par) ) * dsin(nmu*phi0_prev)
 
            
          !! STABLE formulation close to the free surface:  
          !! note: theta = eta - eta0, Psi_s = H + eps*eta0
          !!      
          !! Psi = z + H - eps*G
          !! Psi = z + H - eps*theta + eps*theta - eps*G =
          !!     = z + H - eps*eta + eps*eta0 + eps*( theta - G) = 
          !!     = z - eps*eta + Psi_s + eps*( theta - G ) 
          !!     =                psi0 + eps*( theta - G )
          !!
          psi0 = psi0 + eps * wave_modes%mode(n) * &
                   ( 1.d0 - sinh_o_sinh(nmu*psi0_prev,nmu*par) ) * dcos(nmu*phi0_prev) 
         
        ENDDO
     
                  
        err_iter = (phi0-phi0_prev)**2 + (psi0-psi0_prev)**2
        err_iter = sqrt( err_iter )/par  !! NOTE: dimensionless value
        
        phi0_prev = phi0
        psi0_prev = min( psi0, par)
          
      ENDDO     
      
      phi0 = Flim(phi0)
      psi0 = Flim(psi0)
      
      
      IF (iter.GE.iter_max) THEN
        write(*,*) 'Convergence problems for evaluating phi_0 and psi_0'
        write(*,*) 'PHYSICAL POINTS:  x0 = ', x0, ' z0 = ',z0
        write(*,*) 'transformed plane: phi =', phi0, ' psi0 = ', psi0
        write(*,*) 'Psi_s = ', par
        write(*,*) 'iteration erro: ', err_iter
        STOP
      ENDIF   
      
      
      !! Velocity solution "u" at point (x0,z0)
      dxdphi = 1.d0
      dzdphi = 0.d0
      DO n=1,n_modes 
      
        nmu  = dble(n)*mu
        
        dxdphi = dxdphi + eps * nmu * wave_modes%mode(n) * &
                  ( cosh_o_cosh(nmu*psi0,nmu*par)/tanh(nmu*par) ) * cos(nmu*phi0)
        
        dzdphi = dzdphi - eps * nmu * wave_modes%mode(n) * &
                                ( sinh_o_sinh(nmu*psi0,nmu*par) ) * sin(nmu*phi0)
 
      ENDDO
      
      dxdphi = Flim(dxdphi)
      dzdphi = Flim(dzdphi)
      
      vel_sol(1) = dxdphi/( dxdphi**2 + dzdphi**2 )  !! NOTE: DIMENSIONLESS value (u)
      vel_sol(1) = vel_sol(1) - 1.d0                 !! Back to the earth (i.e. fixed) frame of reference
      
      vel_sol(2) = dzdphi/( dxdphi**2 + dzdphi**2 )  !! NOTE: DIMENSIONLESS value (w)
      
      
      RETURN
      END SUBROUTINE inversion_xz
!-----------------------------------------------------------------------  
!=======================================================================
!-
!=======================================================================
!-----------------------------------------------------------------------  
      FUNCTION cosh_o_cosh(a,b)
      
      !! numerically well-posed formulation for the function
      !! cosh(a)/cosh(b)  with  0<=a<=b
      !!
      !! funz = exp(-(b-a)) * (1 + exp(-2*a))/ (1 + exp(-2*b))
      !! where (b-a) < 2a < 2b
      
      IMPLICIT NONE
      
      real(8)           :: a,b,cosh_o_cosh
      real(8),parameter :: tol_cut=100.d0
      
      
      IF (ABS(a-b).GT.tol_cut) THEN
      
        cosh_o_cosh = 0.d0

      ELSE
      
        IF (2.d0*ABS(a).GT.tol_cut) THEN
          cosh_o_cosh = exp(a-b) 
        ELSEIF (2.d0*ABS(b).GT.tol_cut) THEN
          cosh_o_cosh = exp(a-b)*( 1.d0 + exp(-2.d0*a))
        ELSE
          cosh_o_cosh = exp(a-b)*( 1.d0 + exp(-2.d0*a))/( 1.d0 + exp(-2.d0*b))
        ENDIF
        
      ENDIF

      
      RETURN
      END FUNCTION
!----------------------------------------------------------------------- 
!======================================================================= 
!-
!=======================================================================
!-----------------------------------------------------------------------  
      FUNCTION sinh_o_sinh(a,b)
      
      !! numerically well-posed formulation for the function
      !! sinh(a)/sinh(b)   with   0<=a<=b
      
      !! funz = exp(-(b-a)) *  (1 - exp(-2*a))/ (1 - exp(-2*b))
      !! where (b-a) < 2a < 2b
      
      IMPLICIT NONE
      
      real(8)            :: a,b,sinh_o_sinh
      real(8), parameter :: tol_cut=100.d0 
      
      
      IF (ABS(a-b).GT.tol_cut) THEN
      
        sinh_o_sinh = 0.d0

      ELSE
      
        IF (2.d0*ABS(a).GT.tol_cut) THEN
          sinh_o_sinh = exp(a-b) 
        ELSEIF (2.d0*ABS(b).GT.tol_cut) THEN
          sinh_o_sinh = exp(a-b)*( 1.d0 - exp(-2.d0*a))
        ELSE
          sinh_o_sinh = exp(a-b)*( 1.d0 - exp(-2.d0*a))/( 1.d0 - exp(-2.d0*b))
        ENDIF
        
      ENDIF

      
      RETURN
      END FUNCTION
!----------------------------------------------------------------------
!=======================================================================
!
!=======================================================================
!-----------------------------------------------------------------------  
      function inverse_disp_rel(omega_bar,eps_bar)
      ! Inversion of the dispersion relation through bisection method
      ! NOTE: using DIMENSIONLESS variables

      IMPLICIT NONE 

      real(8) :: inverse_disp_rel       
      real(8) :: omega_bar
      real(8) :: mu_bar1,mu_1
      real(8) :: mu_bar2,mu_2
      real(8) :: mu_bar3,mu_3
      real(8) :: eps_bar,eps_1,eps_2,eps_3
      real(8) :: eps_max
      real(8) :: eq1,eq2,eq3
      !!
      real(8),parameter    :: toll = 1.E-9
      integer(4),parameter :: cont_max=5000  !max number of iterations
      integer(4)           :: cont
      real(8),external     :: disp_rel
      
      
      ! lower extremum
      mu_bar1  = (omega_bar/1.3d0)**2
      mu_1     = tanh(mu_bar1)
      eps_1    = eps_bar * mu_bar1/mu_1
      !---
      eps_max =  -0.0427 *mu_1**3 + 0.1960 * mu_1**2  - 0.1248*mu_1 + 0.4091
      if (eps_1.GT.eps_max) then
        write(*,*) 'Problem in FUNCTION "inverse_dips_rel"'
        write(*,*) 'P1: eps > eps_MAX '
        stop
      endif
      !---
      eq1 = omega_bar - disp_rel(mu_bar1,eps_bar)
      
      
      ! upper extremum
      mu_bar2 = 0.5d0*( omega_bar**2 + omega_bar * sqrt(omega_bar**2 + 4.0) )
      mu_2    = tanh(mu_bar2)
      eps_2   = eps_bar * mu_bar2/mu_2
      !---
      eps_max =  -0.0427 * mu_2**3 + 0.1960 * mu_2**2  - 0.1248*mu_2 + 0.4091
      if (eps_2.GT.eps_max) then
        write(*,*) 'Problem in FUNCTION "inverse_dips_rel"'
        write(*,*) 'P2: eps > eps_MAX '
        stop
      endif
      !---
      eq2 = omega_bar - disp_rel(mu_bar2,eps_bar)
      
      
      !---------------------------
      if (eq1*eq2.GE.0) then
         write(*,*) 'Something unexpected!'
         write(*,*) 'eq1: ', eq1
         write(*,*) 'eq2: ', eq2 
      endif
      !---------------------------
      
      
      
      mu_bar3 = 0.5*(mu_bar1+mu_bar2)
      mu_3    = tanh(mu_bar3)
      eps_3   = eps_bar * mu_bar3/mu_3
      !---
      eps_max =  -0.0427 * mu_3**3 + 0.1960 *mu_3**2  - 0.1248*mu_3 + 0.4091
      if (eps_3.GT.eps_max) then
        write(*,*) 'Problem in FUNCTION "inverse_dips_rel"'
        write(*,*) 'P3: eps > eps_MAX '
        stop
      endif
      !---
      eq3 = omega_bar - disp_rel(mu_bar3,eps_bar)
      
      
      

      cont=1
      do while ((eq1*eq2.LT.0.d0).AND.(cont.LT.cont_max).AND.(abs(eq3)>toll))
          
          cont=cont+1
          
          mu_bar3 = 0.5*(mu_bar1+mu_bar2)
          mu_3    = tanh(mu_bar3)
          eps_3   = eps_bar * mu_bar3/mu_3
          !---
          eps_max =  -0.0427 * mu_3**3 + 0.1960 * mu_3**2  - 0.1248*mu_3 + 0.4091
          if (eps_3.GT.eps_max) then
            write(*,*) 'Problem in FUNCTION "inverse_dips_rel"'
            write(*,*) 'P4: eps > eps_MAX '
            stop     
          endif
          !---
          eq3 = omega_bar - disp_rel(mu_bar3,eps_bar)    
          
      
          if (eq3*eq1.GT.0.d0) then              
            mu_bar1 = mu_bar3                    
          elseif (eq3*eq2.GT.0) then
            mu_bar2 = mu_bar3              
          endif
          
      end do
      
      
      inverse_disp_rel = mu_bar3
      
      
      return
      END FUNCTION
!!----------------------------------------------------------------------  
!=======================================================================
!
!=======================================================================
!!----------------------------------------------------------------------  
      FUNCTION disp_rel(mu_bar,eps_bar)
      
      ! compute the DIMENSIONLESS circular frequency "omega_bar"
      ! using the regression formula described in Antuono (2022)
 
      IMPLICIT NONE  

      real(8) :: disp_rel 
      real(8) :: mu_bar,eps_bar
      real(8) :: mu,eps,Fr 
      real(8) :: eps_max
      real(8) :: a,b,c,d,e
     
      
      mu  = tanh(mu_bar) 
      eps = (mu_bar/mu) * eps_bar
      
      eps_max =  -0.0427 * (mu**3) + 0.1960 * (mu**2) - 0.1248*mu + 0.4091
      
      IF (eps.GE.eps_max) THEN
        write(*,*) 'Problems in FUNCTION "disp_rel"'
        write(*,*) 'The wave amplitude is too large!'
        write(*,*) 'P5: eps > eps_MAX ' 
      ENDIF
      
      
      !! computing the Froude number // see Antuono (2022) equations (6.6)-(6.9) 
      a = -15.2479d0 * eps**3 + 10.9266d0 * eps**2 -  3.0256d0 * eps + 1.1712d0
      b =  26.9416d0 * eps**3 - 17.1964d0 * eps**2 +  5.6883d0 * eps - 2.2527d0 
      c = -18.2407d0 * eps**3 + 11.7909d0 * eps**2 -  4.4316d0 * eps + 1.4138d0
      !---
      d =  141.0995d0 * eps**4 -119.6816d0 * eps**3 + 38.0355d0 * eps**2 -5.6719d0 * eps + 0.3527d0
      e = -119.3300d0 * eps**4 + 84.4075d0 * eps**3 - 19.0126d0 * eps**2 -0.0593d0 * eps + 0.3693d0
 
      Fr = 1.d0 + eps**2 * ( a + b * mu + c * mu**2 )/( eps + d * mu + e * mu**2 + 1.E-12 ) 
      
      disp_rel = sqrt( mu * mu_bar ) * Fr
       
      
      return
      END FUNCTION
!!----------------------------------------------------------------------  
!=======================================================================
!
!=======================================================================      
!!######################################################################
