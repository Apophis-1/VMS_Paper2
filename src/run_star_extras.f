! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib
      use chem_def

      use const_def

      
      implicit none
      real(dp) :: flag = 0  
      real(dp) :: mdot_check = 0, L_iter = 0, M_iter = 0
      real(dp) :: Gamma_e_switch = 0, Gamma_e_switch_old = 0                 
      real(dp) :: Gamma_e = 0, Gamma_e_old = 0, Gamma = 0, Omega_Gamma = 0  

      integer :: do_iter = 1, check_once = 1 
      integer :: switch_to_V01 = 1         ! boolean to control switches. Begins with true
      integer :: switch_to_VMS = 1         ! boolean to control switches. Begins with true



 
      ! these routines are called by the standard run_star check_model
      contains


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
         s% other_wind => VMS_framework
         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras =.false.
            
      end subroutine extras_controls
      
      ! None of the following functions are called unless you set their
      ! function point in extras_control.
      ! these routines are called by the standard run_star check_model




      subroutine VMS_framework(id, L_phot, M_phot, R_surf, T_phot, w, ierr)

         use crlibm_lib
         type (star_info), pointer :: s
         integer, intent(in) :: id
         real(dp), intent(in) :: L_phot, M_phot, T_phot, R_surf! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

   ! -------------------------------------------------------- Variable decleration ------------------------------------------------------------

         integer :: use_V01, use_VMS 

         real(dp) ::  L1, M1, T1, R1, X, Y, Z, w_vink                     ! surface properties
         real(dp) ::  Z_init, X_c, X_c_init, Y_c_init, Zsolar, const_k          
         real(dp) ::  vesc, vesc_eff, vinf, vinf_fac
         real(dp) ::  M, Gamma_e_iter, R_iter, eta_vink, eta_switch, vcrit_iter, v_sq_iter, Boost_iter, Gamma_iter 
         real(dp) ::  diff_old, diff_new, w1, w2
         real(dp) ::  vcrit, vsurf, v_vcrit_sq, Boost, alpha, Teff_jump
         real(dp) ::  w_all, w_RSG, w_WR, alfa_WR, mass_loss

         integer :: nz, i
         call get_star_ptr(id,s,ierr)
         w = 0
         ierr = 0

         use_V01 = 0
         use_VMS = 0

         Gamma_e_old = Gamma_e
         Gamma_e_switch_old = Gamma_e_switch

         L1 = L_phot                                                    ! surface L, M, T and R
         M1 = M_phot
         T1 = T_phot
         R1 = R_surf

         nz = s% nz
         X = s% xa(s%net_iso(ih1),1)                                    ! surface X
         Y = s% xa(s%net_iso(ihe4),1)                                   ! surface Y
         Z = 1 - (Y + X)                                                ! surface Z
         
         Z_init = s% initial_z                                          ! mdot scaling with Z_fe for Vink recipe
         X_c = s% xa(s%net_iso(ih1),nz)                                 ! central X
         Y_c_init = 0.24 + 2*Z_init
         X_c_init = 1 - Y_c_init - Z_init
         Gamma_e = 10**(-4.813)*(1+X)*(L1/M1)*(Msun/Lsun)               ! electron scattering Eddington parameter

         Zsolar = 0.019
         const_k = (clight*Msun*1d5)/(Lsun*3600d0*24d0*365d0)           ! constant to evaluate eta


   ! ----------------------------------------- Store and retrieve switch properties for He burning --------------------------------------------

         if (s% x_ctrl(1) == 0) then                                  
              open(unit = 1,file = 'LM_switch.txt',status='replace',form='formatted')       
              write(1,*) Gamma_e_switch, L_iter, M_iter, mdot_check, do_iter
    	      close(1)
         else if (s% x_ctrl(1) == 1 .and. check_once == 1) then          ! For helium burning, retrieve the L_switch, and do not perform the iteration at start of step
              open (unit = 2, file = 'LM_switch.txt', status = 'old')
              read(2,*) Gamma_e_switch, L_iter, M_iter, mdot_check, do_iter
              close(2)
              Gamma_e_switch_old = Gamma_e_switch
              Gamma_e_old = Gamma_e
              check_once = 0
         end if



   ! ------------------------------------------- iteration procedure to estimate L_switch -----------------------------------------------------
                 
         call eval_Vink01_wind(w, vinf_fac, Teff_jump, L1, M1)           ! evaluate mdot_vink
         w_vink = w


         if (X_c/X_c_init < 0.999 .and. do_iter == 1) then   

             do i = 1, 3000, 1
                  L_iter = 5d0 + 0.001*i
                  if (L_iter < 6.5) then
                  	call eval_homogeneous_mass_low(M)
                  else
                        call eval_homogeneous_mass_high(M)
                  end if
                  L_iter = exp10_cr(L_iter) * Lsun
                  M_iter = M * Msun 
                  Gamma_e_iter = 10**(-4.813)*(1+X)*(L_iter/M_iter)*(Msun/Lsun)              ! electron scattering gamma_e of iteration to correct M
                  R_iter = pow_cr((L_iter)/(4*pi*T1*T1*T1*T1*boltz_sigma), 0.5d0)
                  
                  vcrit_iter = pow_cr((2d0/3d0)*standard_cgrav*M_iter/R_iter, 0.5d0)/1d5
                  vsurf = s% v_rot_avg_surf/1d5 
                  v_sq_iter = (vsurf/vcrit_iter)**2

                  alpha = 0.52d0
                  Gamma_iter = s% photosphere_opacity * L_iter/(4*pi*standard_cgrav*clight*M_iter) 
                  Boost_iter = (pow_cr((1-Gamma_iter), 1/alpha - 1))/(pow_cr((1-Gamma_iter - (4d0/9d0)*v_sq_iter), 1/alpha - 1))

                  call eval_Vink01_wind(w, vinf_fac, Teff_jump, L_iter, M_iter)
                  mdot_check = w
                  
                  vesc = pow_cr(2d0*standard_cgrav*M_iter/R_iter, 0.5d0)/1d5
                  vesc_eff = pow_cr(2d0*standard_cgrav*M_iter*(1-Gamma_e_iter)/R_iter, 0.5d0)/1d5
                  vinf =  vinf_fac * vesc_eff * pow_cr(Z/Zsolar,0.20d0) 
                  
                  eta_vink = const_k*(mdot_check*vinf)/(L_iter/Lsun)
                  eta_switch = 0.75/(1+(vesc*vesc)/(vinf*vinf)) 
                           
                  if (eta_vink > eta_switch) then                                        ! Switch is reached when the etas are equal      
                       exit
                  end if               
              end do
  
              Gamma_e_switch = (10**-4.813)*(1+X)*(L_iter/M_iter)*(Msun/Lsun)           

         end if

         if (X_c/X_c_init < 0.999) then

              diff_old = Gamma_e_old - Gamma_e_switch_old
              diff_new = Gamma_e - Gamma_e_switch
              
              if (diff_new < 0) then            ! sets what mass loss recipe to use
                  use_V01 = 1
              else if (diff_new > 0) then
                  use_VMS = 1
                  do_iter = 0
              end if
         

              if (diff_old < 0 .and. diff_new > 0) then              ! Additional switch condition to go from V01 to VMS
                  switch_to_VMS = 1
                  switch_to_V01 = 0
              else if (diff_old > 0 .and. diff_new < 0) then         ! Additional switch condition to go from VMS back to V01 
                  switch_to_V01 = 1
                  switch_to_VMS = 0
              end if  
         end if
         

   ! -------------------------------------------------------- Mass loss for different stages --------------------------------------------------

         if (use_V01 == 1 .and. switch_to_V01 == 1) then
             call eval_Vink01_wind(w, vinf_fac, Teff_jump, L1, M1)
             flag = 1
             w_all = w
         else if (use_VMS == 1 .and. switch_to_VMS == 1) then
             call eval_Vink11_wind(w, L1, M1)
             w_all = w
  	     flag = 2
         end if
         
         call eval_de_Jager_wind(w)                   ! de Jager recipe for below 4,000 K
         w_RSG = w
                 
         call eval_sv2020_WR_mdot(w1)                  ! Sander and Vink 2021 recipe for temperatures above 100 kK, (Also has Vink 2017 now, but was not used in the paper)
         call eval_v2017_WR_mdot(w2)  
         w_WR = MAX(w1,w2)  

         


         if (X_c/X_c_init > 0.999) then               ! mass loss during PMS is set to zero
              mass_loss = 0
              flag = 0
         else if (T1 < 9d4 .and. T1 > 4d3) then       ! mass loss is the new VMS implementation in the temp range, regardless of the fuel burning in the core
              mass_loss = w_all
         else if (T1 < 4d3) then
              mass_loss = w_RSG
              flag = 4
         else if (T1 > 1d5) then
              mass_loss = w_WR
              flag = 3
         else if (T1 > 9d4 .and. T1 < 1d5) then
              alfa_WR = (T1 - 9d4) / (1d4)
              mass_loss = alfa_WR * w_WR + (1 - alfa_WR) * w_all
         
         end if

   ! -------------------------------------------------------- For rotation only ---------------------------------------------------------------

         vcrit = pow_cr((2d0/3d0)*standard_cgrav*M1/R1, 0.5d0)/1d5
         vsurf = s% v_rot_avg_surf/1d5 
         v_vcrit_sq = (vsurf/vcrit)**2

         alpha = 0.52d0
         Gamma = s% photosphere_opacity * L1/(4*pi*standard_cgrav*clight*M1) 
         Boost = (pow_cr((1-Gamma), 1/alpha - 1))/(pow_cr((1-Gamma - (4d0/9d0)*v_vcrit_sq), 1/alpha - 1))
         Omega_Gamma = Gamma/(1-(4d0/9d0)*v_vcrit_sq)
         mass_loss = Boost * mass_loss

   ! -------------------------------------------------------- Outputs for terminal ------------------------------------------------------------

         write(*,*) 'flag = ', flag
         !write(*,*) 'use_V01 = ', use_V01
         !write(*,*) 'use_VMS = ', use_VMS
         !write(*,*) 'switch_to_V01 = ', switch_to_V01
         !write(*,*) 'switch_to_VMS = ', switch_to_VMS
         !write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         !write(*,*) 'L switch = ', log10_cr(L_iter/Lsun)
         !write(*,*) 'M switch = ', M_iter/Msun
         !write(*,*) 'mdot switch = ', log10_cr(mdot_check)
         !write(*,*) 'vinf = ', vinf_fac 
         !write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         !write(*,*) 'L model = ', log10_cr(L1/Lsun)
         !write(*,*) 'M model = ', M1/Msun
         !write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         !write(*,*) 'eta_vink = ', eta_vink
         !write(*,*) 'eta_switch = ', eta_switch
         !write(*,*) 'Gamma_e = ', Gamma_e
         !write(*,*) 'Gamma_e at switch = ', Gamma_e_switch
         !write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         !write(*,*) 'Rotation boost = ', Boost
         !write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,*) 'Mdot from V01 = ', log10_cr(w_vink)
         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,*) 'Mdot used in evolution = ', log10_cr(mass_loss)
         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


         w = s% Dutch_scaling_factor * mass_loss







         contains


         subroutine eval_homogeneous_mass_low(M)               ! For logL < 6.5
             real(dp), intent(out) :: M
             real(dp) :: f, logL, g
 
             real(dp), parameter :: F1 = 4.026
             real(dp), parameter :: F2 = 4.277
             real(dp), parameter :: F3 = -1
             real(dp), parameter :: F4 = 25.48
             real(dp), parameter :: F5 = 36.93
             real(dp), parameter :: F6 = -2.792
             real(dp), parameter :: F7 = -3.226
             real(dp), parameter :: F8 = -5.317
             real(dp), parameter :: F9 = 1.648

             !logL = log10_cr(L_iter/Lsun)
             f = F4 + F5*X + F6*pow2(X) + (F7 + F8*X)*L_iter
             g = pow_cr(f, 0.5d0)
             M = (F1 + F2*X + F3*g)/(1+F9*X)
             M = pow_cr(10d0, M)
             
        
         end subroutine eval_homogeneous_mass_low              

         subroutine eval_homogeneous_mass_high(M)              ! For logL > 6.5
             real(dp), intent(out) :: M
             real(dp) :: f, logL, g
 
             real(dp), parameter :: F1 = 10.05
             real(dp), parameter :: F2 = 8.204
             real(dp), parameter :: F3 = -1
             real(dp), parameter :: F4 = 151.7
             real(dp), parameter :: F5 = 254.5
             real(dp), parameter :: F6 = -11.46
             real(dp), parameter :: F7 = -13.16
             real(dp), parameter :: F8 = -31.68
             real(dp), parameter :: F9 = 2.408

             !logL = log10_cr(L_iter/Lsun)
             f = F4 + F5*X + F6*pow2(X) + (F7 + F8*X)*L_iter
             g = pow_cr(f, 0.5d0)
             M = (F1 + F2*X + F3*g)/(1+F9*X)
             M = pow_cr(10d0, M)          
        
         end subroutine eval_homogeneous_mass_high




         subroutine eval_sv2020_WR_mdot(w)
            real(dp), intent(out) :: w
            real(dp) :: log10w, a, c, gamma_eb, log10w_off
            include 'formats'

            a = 2.932 
            c = -0.44*log10_cr(Z_init/Zsolar) + 9.15
            gamma_eb = -0.324*log10_cr(Z_init/Zsolar) + 0.244
            log10w_off = 0.23*log10_cr(Z_init/Zsolar) - 2.61

            log10w = a*log10_cr(-log10_cr(1-Gamma_e)) - 0.3010*pow_cr(Gamma_eb/Gamma_e, c) + log10w_off
            w = exp10_cr(log10w)

         end subroutine eval_sv2020_WR_mdot


         subroutine eval_v2017_WR_mdot(w)           ! This part of the recipe was not used in the paper
            real(dp), intent(out) :: w
            real(dp) :: log10w
            
            include 'formats'
            
            log10w = -13.3 + 1.36 * log10_cr(L1/Lsun) + 0.61 * log10_cr(Z_init/Zsolar)
            w = exp10_cr(log10w)

         end subroutine eval_v2017_WR_mdot


         subroutine eval_de_Jager_wind(w)
            ! de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259.
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = 1.769d0*log10_cr(L1/Lsun) - 1.676d0*log10_cr(T1) - 8.158d0
            w = exp10_cr(log10w)
         end subroutine eval_de_Jager_wind



         subroutine eval_Vink01_wind(w, vinf_fac, Teff_jump, L, M)
            real(dp), intent(inout) :: w, vinf_fac, Teff_jump
            real(dp), intent(in) :: L, M
            real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc_1, vinf_div_vesc_2, vinf_div_vesc, f_z

            f_z = 0.85d0

 
            ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
            Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10_cr(Z_init/Zsolar)))
            dT = 2000d0
            if (T1 > Teff_jump + dT) then
                alfa = 1
            else if (T1 < Teff_jump - dT) then
                alfa = 0
            else
                alfa = (T1 - (Teff_jump - dT)) / (2*dT)
            end if
            

            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc_1 = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc_1*pow_cr(Z/Zsolar,0.20d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10_cr(L/Lsun/1d5) &
                  - 1.313d0*log10_cr(M/Msun/30) &
                  - 1.226d0*log10_cr(vinf_div_vesc/2d0) &
                  + f_z*log10_cr(Z_init/Zsolar)
               w1 = exp10_cr(logMdot)
            else
               vinf_div_vesc_1 = 0
               w1 = 0
            end if

            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc_2 = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc_2*pow_cr(Z/Zsolar,0.20d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10_cr(L/Lsun/1d5) &
                  - 1.339d0*log10_cr(M/Msun/30) &
                  - 1.601d0*log10_cr(vinf_div_vesc/2d0) &
                  + f_z*log10_cr(Z_init/Zsolar)
               w2 = exp10_cr(logMdot)
            else
               vinf_div_vesc_2 = 0
               w2 = 0
            end if

            w = alfa*w1 + (1 - alfa)*w2
            vinf_fac = alfa*vinf_div_vesc_1 + (1 - alfa)*vinf_div_vesc_2

         end subroutine eval_Vink01_wind


         subroutine eval_Vink11_wind(w, L, M)
            real(dp), intent(inout) :: w
            real(dp), intent(in) :: L, M
            real(dp) :: logMdot

            logMdot = &
                  log10_cr(mdot_check) &        
                  + 4.77d0*log10_cr(L/L_iter) & 
                  - 3.99d0*log10_cr(M/M_iter)
            w = exp10_cr(logMdot)

         end subroutine eval_Vink11_wind

      end subroutine VMS_framework



      
      ! None of the following functions are called unless you set their
      ! function point in extras_control.
      
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         !real(dp) :: frac, vct30, vct100
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if




      end function extras_startup
      

      integer function extras_start_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr, nz
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going

         
         nz = s% nz

         if (s% Teff < 7d3 .and. s% xa(s%net_iso(ih1),nz) < 0.7 .and. s% x_ctrl(1) == 0) then
               extras_check_model = terminate
         end if 
         if (s% xa(s%net_iso(ih1),nz) < 0.7 .and. s% time_step < 5 .and. s% x_ctrl(1) == 0) then
               extras_check_model = terminate
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 5
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         real(dp), parameter :: frac = 0.90
         integer :: i
         real(dp) :: edot, edot_partial
 
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
         edot = dot_product(s% dm(1:s% nz), s% eps_nuc(1:s% nz))

         ! the center of the star is at i = s% nz and the surface at i = 1 .
         ! so go from the center outward until 90% of the integrated eps_nuc
         ! is enclosed.  exit and then i will contain the desired cell index.
         edot_partial = 0
         do i = s% nz, 1, -1
            edot_partial = edot_partial + s% dm(i) * s% eps_nuc(i)
            if (edot_partial .ge. (frac * edot)) exit
         end do

         ! note: do NOT add these names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         names(1) = "flag"
         vals(1) = flag

         names(2) = "Gamma_e"
         vals(2) = Gamma_e

         names(3) = "Gamma_e_switch"
         vals(3) = Gamma_e_switch

         names(4) = "Gamma"
         vals(4) = Gamma

         names(5) = "Omega_Gamma"
         vals(5) = Omega_Gamma



      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         if (n /= 1) stop 'data_for_extra_profile_columns'

         
         
      end subroutine data_for_extra_profile_columns

      subroutine how_many_extra_history_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols=0
      end subroutine how_many_extra_history_header_items
      
      subroutine data_for_extra_history_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra history header item
      !set num_cols=1 in how_many_extra_history_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items


      subroutine how_many_extra_profile_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols = 0
      end subroutine how_many_extra_profile_header_items
      
      subroutine data_for_extra_profile_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra profile header item
      !set num_cols=1 in how_many_extra_profile_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg
         num_ints = i
         
         i = 0
         ! call move_dbl
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info




      end module run_star_extras
      
