#include "../include/dprec.fh"
#include "../include/assert.fh"

module qm2_extern_mlmm_module
! ----------------------------------------------------------------------
! Interface for torchtorchani based ML and ML/MM simulations in SANDER
!
! The ML/MM coupling can be performed at the mechanical embedding level
! using fixed charges from the topology file, or by computing QM charges
! for the ML subsystem in vacuo, using PSI4
!
! Author: Jonathan Semelak
!
! Based on qm2_extern_orc_module.F90 and qm2_extern_lio_module.F90
! Date: February 2024
! ----------------------------------------------------------------------

  implicit none

  private
  public :: get_mlmm_forces

  type mlmm_nml_type
     character(len=256) :: charges_functional
     character(len=256) :: charges_basis
     character(len=256) :: charges_model
     character(len=256) :: mlmm_patch_path
     character(len=256) :: psi4_path
     character(len=256) :: python_executable
     integer :: psi4_nthreads
     logical :: polarize_qm_charges
     logical :: distort_qm_energy
     logical :: write_xyz
     logical :: write_forces
     logical :: write_charges
     logical :: use_charges_derivatives
     double precision :: pol_H
     double precision :: pol_He
     double precision :: pol_Li
     double precision :: pol_Be
     double precision :: pol_B
     double precision :: pol_C
     double precision :: pol_N
     double precision :: pol_O
     double precision :: pol_F
     double precision :: pol_Ne
     double precision :: pol_Na
     double precision :: pol_Mg
     double precision :: pol_Al
     double precision :: pol_Si
     double precision :: pol_P
     double precision :: pol_S
     double precision :: pol_Cl
     double precision :: pol_Ar
  end type mlmm_nml_type

contains

  ! ----------------------------------------------------------------
  ! Get ML energy and forces from torchani
  ! ----------------------------------------------------------------
    subroutine get_mlmm_forces( nstep, ntpr_internal , nqmatoms, qmcoords,&
      nclatoms, clcoords, escf, dxyzqm, dxyzcl)

      use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS

      use qmmm_module, only : qmmm_struct, qmmm_nml

      implicit none

      integer, intent(in) :: nqmatoms             ! Number of QM atoms
      integer, intent(in) :: nstep                ! MDIN input
      integer, intent(in) :: ntpr_internal        ! MDIN input
      integer, intent(in) :: nclatoms             ! Number of MM atoms
      _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
      _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in a
      _REAL_, intent(out) :: escf                 ! torchani energy
      _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! torchani forces
      _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! torchani MM forces
      _REAL_ dxyzqm_qmmm(3,nqmatoms)   ! torchani QMMM forces on QM atoms
      _REAL_ alpha(nqmatoms)        ! polarizabilities on QM atoms
      _REAL_ qmcharges(nqmatoms) ! QM atom charges
      _REAL_ dqmcharges(3,nqmatoms,nqmatoms) ! QM atom charges forces
    ! Note that atomic_charge_derivatives is an array of shape [3, num_atoms, num_atoms]
    ! where the element [i, j, k] is the derivative of the **charge on
    ! k-th atom** with respect to the **i-th position of the j-th atom**
      _REAL_ eqmmm ! EQMMM energy
      _REAL_ psi4_escf !psi4 energy for QM region
      character(len=2) :: qm_atomic_symbols(nqmatoms)
      type(mlmm_nml_type), save     :: mlmm_nml
      logical, save                :: first_call = .true.
      logical first_step
      integer i

      ! Setup on first call
      first_step=first_call
      if ( first_call ) then
         first_call = .false.
         write (6,'(/,a,/)') '  >>> RUNNING ML/MM SIMULATION <<<'
         call get_namelist_ml(mlmm_nml)
         call print_namelist(mlmm_nml)
         if (  qmmm_nml%qmmm_int == 0) then
             write(6,*) "ML/MM: ML/MM electrostatic coupling will be calculated by SANDER"
         else if (  qmmm_nml%qmmm_int == 5) then
             write(6,*) "ML/MM: ML/MM electrostatic coupling will be calculated by SANDER"
             write(6,*) "ML/MM: The mechanical embedding approximation will be used"
             write(6,*) "ML/MM: Charges of the ML subsystem will be read from the topology file"
         else if (  qmmm_nml%qmmm_int == 2) then
               write(6,*) "ML/MM: ML/MM electrostatic coupling will be calculated by this interface"
               write(6,*) "ML/MM: The mechanical embedding approximation will be used"
               if (trim(mlmm_nml%charges_model)/='AMBER_CHARGES') then
                  write(6,*) "ML/MM: Charges for the ML subsystem will be computed with PSI4"
                  if (mlmm_nml%use_charges_derivatives) then
                      write(6,*) "ML/MM: Charges derivatives will be obtained numberically"
                  else
                      write(6,*) "ML/MM: Charges derivatives will be not calculated (not recommended)"
                  end if
               else
                 write(6,*) "ML/MM: Charges of the ML subsystem will be read from the topology file"
               endif
               if (mlmm_nml%polarize_qm_charges) write(6,*) "ML/MM: A correction term to account for the polarization of the charges of the ML subsystem will be used"
               if (mlmm_nml%distort_qm_energy) write(6,*) "ML/MM: A correction term to account for the distortion of the ML energy will be used"
         end if
      write(6,*) "ML/MM: Energies are printed in kcal/mol"
      end if !first_call

      call get_atomic_symbols(nqmatoms,qmmm_struct%iqm_atomic_numbers,qm_atomic_symbols)

      if ( ntpr_internal > 0 .and. mod(nstep+1, ntpr_internal) == 0 ) then
         write(6,*) ""
         write(6,*) "------------------------------------------------------------------------------"
      end if

      if (mlmm_nml%polarize_qm_charges .and. qmmm_nml%qmmm_int == 2) then
          call get_atomic_polarizability(nqmatoms,qmmm_struct%iqm_atomic_numbers,alpha,mlmm_nml%pol_H,mlmm_nml%pol_He,mlmm_nml%pol_Li,&
          mlmm_nml%pol_Be,mlmm_nml%pol_B,mlmm_nml%pol_C,mlmm_nml%pol_N,mlmm_nml%pol_O,&
          mlmm_nml%pol_F,mlmm_nml%pol_Ne,mlmm_nml%pol_Na,mlmm_nml%pol_Mg,mlmm_nml%pol_Al,&
          mlmm_nml%pol_Si,mlmm_nml%pol_P,mlmm_nml%pol_S,mlmm_nml%pol_Cl,mlmm_nml%pol_Ar,first_step)
      end if

      dxyzqm=0.d0
      dxyzcl=0.d0
      dxyzqm_qmmm=0.d0

      call get_ml_properties(qmcoords,nqmatoms,qm_atomic_symbols,dxyzqm,escf,mlmm_nml%mlmm_patch_path,mlmm_nml%python_executable)

      if ( ntpr_internal > 0 .and. mod(nstep+1, ntpr_internal) == 0 ) then
          write(6,'(1x,"ML/MM: IN VACUO ML ENERGY = ",f14.4)') escf
      end if

      ! Invert dxyzqm because torchani retrieves forces instead of gradients
      dxyzqm=-dxyzqm

      if (qmmm_nml%qmmm_int == 2) then
          if (trim(mlmm_nml%charges_model)=='AMBER_CHARGES') then
             qmcharges=qmmm_struct%qm_resp_charges/18.2223d0
          else
             call get_psi4_charges(qmcoords,nqmatoms,qmcharges,qm_atomic_symbols,psi4_escf,mlmm_nml%mlmm_patch_path,&
                                   mlmm_nml%psi4_path,mlmm_nml%python_executable,mlmm_nml%charges_model,mlmm_nml%charges_functional,mlmm_nml%charges_basis,&
                                   mlmm_nml%psi4_nthreads)
          endif
          call get_mlmm_ener_forces(nstep,ntpr_internal,nqmatoms,nclatoms,dxyzqm,qmcoords,clcoords,dxyzcl,dxyzqm_qmmm,&
                                   qmmm_struct%iqm_atomic_numbers,qmcharges,dqmcharges,&
                                   mlmm_nml%polarize_qm_charges,mlmm_nml%distort_qm_energy,alpha,eqmmm,mlmm_nml%use_charges_derivatives,&
                                   qm_atomic_symbols,mlmm_nml%mlmm_patch_path,mlmm_nml%psi4_path,mlmm_nml%python_executable,&
                                   mlmm_nml%charges_model,mlmm_nml%charges_functional,mlmm_nml%charges_basis,mlmm_nml%psi4_nthreads)

          if (mlmm_nml%write_charges) then
             if (first_step) then
                 open(unit=271921, file='charges_qm_region.dat', status='REPLACE')
             else
                 open(unit=271921, file='charges_qm_region.dat', status='OLD', position='APPEND')
             endif
             do i = 1, nqmatoms
                 write(271921, '(F15.7)') qmcharges(i)
             end do
             write(271921, *)
             close(271921)
          endif
       end if

       if (mlmm_nml%write_forces) then
           if (first_step) then
               open(unit=271920, file='forces_qm_region.dat', status='REPLACE')
           else
               open(unit=271920, file='forces_qm_region.dat', status='OLD', position='APPEND')
           endif
           do i = 1, nqmatoms
               write(271920, '(3F15.7)') -dxyzqm(1, i),-dxyzqm(2, i),-dxyzqm(3, i)
           end do
           write(271920, *)
           close(271920)
           if (nclatoms .ge. 1) then
               if (first_step) then
                   open(unit=271923, file='forces_qmmm_mm_region.dat', status='REPLACE')
                   open(unit=271927, file='forces_qmmm_qm_region.dat', status='REPLACE')
               else
                   open(unit=271923, file='forces_qmmm_mm_region.dat', status='OLD', position='APPEND')
                   open(unit=271927, file='forces_qmmm_qm_region.dat', status='OLD', position='APPEND')
               endif
               do i = 1, nclatoms
                   write(271923, '(3F15.7)') -dxyzcl(1, i),-dxyzcl(2, i),-dxyzcl(3, i)
               end do
               write(271923, *)
               close(271923)
               do i = 1, nqmatoms
                  write(271927, '(3F15.7)') -dxyzqm_qmmm(1, i),-dxyzqm_qmmm(2, i),-dxyzqm_qmmm(3, i)
               end do
               write(271927, *)
               close(271927)
           end if
       endif

       if (mlmm_nml%write_xyz) then
           if (first_step) then
               open(unit=271922, file='qm_region.xyz', status='REPLACE')
           else
               open(unit=271922, file='qm_region.xyz', status='OLD', position='APPEND')
           endif
           write(271922,'(I8)') nqmatoms
           write(271922,*) "ML/MM: .xyz trajectory file"
           do i = 1, nqmatoms
               write(271922, '(I3,3F15.7)') qmmm_struct%iqm_atomic_numbers(i), qmcoords(1, i),qmcoords(2, i),qmcoords(3, i)
           end do
           close(271922)
       endif

    if ( ntpr_internal > 0 .and. mod(nstep+1, ntpr_internal) == 0 ) then
         write(6,*) "------------------------------------------------------------------------------"
     end if

    end subroutine get_mlmm_forces

  ! ----------------------------------------------------------------
  ! Read mlmm namelist values from file mdin,
  ! use default values if none are present.
  ! ----------------------------------------------------------------

    subroutine get_namelist_ml(mlmm_nml)
      implicit none
      type(mlmm_nml_type), intent(out) :: mlmm_nml

      character(len=256) :: charges_functional
      character(len=256) :: charges_basis
      character(len=256) :: charges_model
      character(len=256) :: mlmm_patch_path
      character(len=256) :: psi4_path
      character(len=256) :: python_executable
      integer :: psi4_nthreads
      logical :: polarize_qm_charges
      logical :: distort_qm_energy
      logical :: write_xyz
      logical :: write_forces
      logical :: write_charges
      logical :: use_charges_derivatives
      double precision :: pol_H
      double precision :: pol_He
      double precision :: pol_Li
      double precision :: pol_Be
      double precision :: pol_B
      double precision :: pol_C
      double precision :: pol_N
      double precision :: pol_O
      double precision :: pol_F
      double precision :: pol_Ne
      double precision :: pol_Na
      double precision :: pol_Mg
      double precision :: pol_Al
      double precision :: pol_Si
      double precision :: pol_P
      double precision :: pol_S
      double precision :: pol_Cl
      double precision :: pol_Ar

      namelist /ml/ charges_functional,charges_basis,charges_model,mlmm_patch_path,psi4_path,&
                          python_executable,psi4_nthreads,polarize_qm_charges,distort_qm_energy,&
                          write_xyz,write_forces,write_charges,&
                          use_charges_derivatives,pol_H,pol_He,pol_Li,pol_Be,pol_B,pol_C,&
                          pol_N,pol_O,pol_F,pol_Ne,pol_Na,pol_Mg,pol_Al,pol_Si,pol_P,pol_S,&
                          pol_Cl,pol_Ar

      integer :: ierr
      ! Defaults

      charges_functional='wb97x'
      charges_basis='6-31G(d)'
      charges_model='MBIS_CHARGES'
      mlmm_patch_path=''
      psi4_path='psi4'
      python_executable='python'
      psi4_nthreads=1
      polarize_qm_charges = .true.
      distort_qm_energy = .true.
      write_xyz = .false.
      write_forces = .false.
      write_charges = .false.
      use_charges_derivatives = .true.
      pol_H  = -1.d0
      pol_He = -1.d0
      pol_Li = -1.d0
      pol_Be = -1.d0
      pol_B  = -1.d0
      pol_C  = -1.d0
      pol_N  = -1.d0
      pol_O  = -1.d0
      pol_F  = -1.d0
      pol_Ne = -1.d0
      pol_Na = -1.d0
      pol_Mg = -1.d0
      pol_Al = -1.d0
      pol_Si = -1.d0
      pol_P  = -1.d0
      pol_S  = -1.d0
      pol_Cl = -1.d0
      pol_Ar = -1.d0

      ! Read namelist
      rewind 5
      read(5,nml=ml,iostat=ierr)

      if ( ierr > 0 ) then
         call sander_bomb('get_namelist_ml (qm2_extern_mlmm_module)', &
              '&ml namelist read error', &
              'Please check your input.')
      else if ( ierr < 0 ) then
         write(6,*) '&ml amelist read success'
      end if

      mlmm_nml%charges_model=charges_model
      mlmm_nml%charges_functional=charges_functional
      mlmm_nml%charges_basis=charges_basis
      mlmm_nml%mlmm_patch_path=mlmm_patch_path
      mlmm_nml%psi4_path=psi4_path
      mlmm_nml%python_executable=python_executable
      mlmm_nml%psi4_nthreads=psi4_nthreads
      mlmm_nml%polarize_qm_charges=polarize_qm_charges
      mlmm_nml%distort_qm_energy=distort_qm_energy
      mlmm_nml%write_xyz=write_xyz
      mlmm_nml%write_charges=write_charges
      mlmm_nml%write_forces=write_forces
      mlmm_nml%use_charges_derivatives=use_charges_derivatives
      mlmm_nml%pol_H=pol_H
      mlmm_nml%pol_He=pol_He
      mlmm_nml%pol_Li=pol_Li
      mlmm_nml%pol_Be=pol_Be
      mlmm_nml%pol_B=pol_B
      mlmm_nml%pol_C=pol_C
      mlmm_nml%pol_N=pol_N
      mlmm_nml%pol_O=pol_O
      mlmm_nml%pol_F=pol_F
      mlmm_nml%pol_Ne=pol_Ne
      mlmm_nml%pol_Na=pol_Na
      mlmm_nml%pol_Mg=pol_Mg
      mlmm_nml%pol_Al=pol_Al
      mlmm_nml%pol_Si=pol_Si
      mlmm_nml%pol_P=pol_P
      mlmm_nml%pol_S=pol_S
      mlmm_nml%pol_Cl=pol_Cl
      mlmm_nml%pol_Ar=pol_Ar

    end subroutine get_namelist_ml

    subroutine print_namelist(mlmm_nml)
      implicit none
      type(mlmm_nml_type), intent(in) :: mlmm_nml
      write(6,*) 'ML/MM:------------------ML/MM options----------------'
      write(6,*) 'ML/MM:  polarize_qm_charges          ', mlmm_nml%polarize_qm_charges
      write(6,*) 'ML/MM:  distort_qm_energy            ', mlmm_nml%distort_qm_energy
      write(6,*) 'ML/MM:  charges_model                ', mlmm_nml%charges_model
      write(6,*) 'ML/MM:  write_xyz                    ', mlmm_nml%write_xyz
      write(6,*) 'ML/MM:  write_forces                 ', mlmm_nml%write_forces
      write(6,*) 'ML/MM:  write_charges                ', mlmm_nml%write_charges
      write(6,*) 'ML/MM:  use_charges_derivatives      ', mlmm_nml%use_charges_derivatives
      write(6,*) 'ML/MM:  pol_H, pol_C, etc.            specified below'
      write(6,*) 'ML/MM:--------------end ML/MM options----------------'

    end subroutine print_namelist

    ! ----------------------------------------------------------------
    ! Calculates the QMMM forces in the QM and the MM regions analitically.
    ! ----------------------------------------------------------------
    subroutine get_mlmm_ener_forces(nstep,ntpr_internal,nqmatoms,nclatoms,dxyzqm,qmcoords,clcoords,dxyzcl,dxyzqm_qmmm,qmtypes, &
                                    qmcharges,dqmcharges,polarize_qm_charges,distort_qm_energy,alpha, eqmmmtorchani,use_charges_derivatives, &
                                    qm_atomic_symbols,mlmm_patch_path,psi4_path,python_executable, charges_model,&                           
                                    charges_functional, charges_basis, nthreads)

      use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS

      implicit none

      integer, intent(in) :: nstep ! Number of step
      integer, intent(in) :: ntpr_internal ! Printing frequency
      integer, intent(in) :: nqmatoms                       ! Number of QM atoms
      double precision,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
      integer, intent(in) :: nclatoms                       ! Number of MM atoms
      double precision,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
      integer, intent(in) :: qmtypes(nqmatoms)               ! QM atoms atomic numbers
      double precision, intent(in):: qmcharges(nqmatoms)    ! QM atoms charges
      double precision, intent(in) :: dqmcharges(3,nqmatoms,nqmatoms) ! QM atoms charges forces
      ! Note that atomic_charge_derivatives is  an array of shape [3, num_atoms, num_atoms]
      ! where the element [i, j, k] is the derivative of the **charge on
      ! k-th atom** with respect to the **i-th position of the j-th atom**
      double precision, intent(inout) :: dxyzqm(3,nqmatoms) ! QM atoms forces
      double precision, intent(inout) :: dxyzcl(3,nclatoms) ! MM atoms qmmm forces
      double precision, intent(out) :: dxyzqm_qmmm(3,nqmatoms) ! QM atom eqmmm forces
      character(len=2), intent(in) :: qm_atomic_symbols(nqmatoms) ! QM atoms atomic numbers
      character(len=256), intent(in) ::python_executable ! Python executable
      character(len=256), intent(in) ::mlmm_patch_path ! Path to mlmm_patch
      character(len=256), intent(in) ::psi4_path ! Path to psi4
      character(len=256), intent(in) ::charges_model ! Path to psi4
      character(len=256), intent(in) ::charges_functional ! psi4 functional
      character(len=256), intent(in) ::charges_basis ! psi4 basis
      integer, intent(in) :: nthreads ! psi4 nthreads
      logical, intent(in) :: polarize_qm_charges ! QMMM electrostatic contribution
      logical, intent(in) :: distort_qm_energy ! Include distortion correction of QM energy
      double precision distortion_k ! Distortion constant
      logical, intent(in) :: use_charges_derivatives ! consider or not that QM charges depend on QM coordinates
      double precision eqmmmtorchani ! QMMM electrostatic contribution
      double precision distij(nqmatoms,nclatoms) ! QM/MM distance matrix
      double precision, intent(in) :: alpha(nqmatoms) ! QM atoms atomic polarizabilities
      double precision efield(3,nqmatoms) ! Electric field on QM atoms
      double precision epol  !QM/MM polariziton energy (on QM charges)
      double precision edist !QM distortion energy (on QM system)
      double precision :: escf !MAKE IT INOUT <-----------------------------------------------DEBUG
      integer i, j, k

      distortion_k = 0.5d0
      ! Calculates QM/MM distance matrix
      call get_distij_matrix(nqmatoms,nclatoms,qmcoords,clcoords,distij)

      ! Calculates QM/MM coumlomb interaction energy
      eqmmmtorchani=0.d0
      do i=1,nqmatoms
        do j=1,nclatoms
          eqmmmtorchani = eqmmmtorchani + qmcharges(i)*clcoords(4,j)/distij(i,j)
        end do
      end do
      ! Calculates QM/MM polarization energy
      epol=0.d0
      edist=0d0
      ! Only polarization on the QM atoms is accounted
      if (polarize_qm_charges .or. distort_qm_energy) then
        call get_elec_field(nqmatoms,nclatoms,qmcoords,clcoords,distij,efield)
        do i=1,nqmatoms
          epol=epol+alpha(i)*(efield(1,i)**2+efield(2,i)**2+efield(3,i)**2)
        end do
        epol=-0.5d0*epol
        edist=-distortion_k*epol
      end if

      eqmmmtorchani = eqmmmtorchani + epol

      if ( ntpr_internal > 0 .and. mod(nstep+1, ntpr_internal) == 0 ) then
         if (distort_qm_energy) then
             write(6,'(1x,"ML/MM: DISTORTION ENERGY = ",f14.4)') edist*CODATA08_AU_TO_KCAL/CODATA08_A_TO_BOHRS
         end if
         if (polarize_qm_charges) then
             write(6,'(1x,"ML/MM: POLARIZATION ENERGY = ",f14.4)') epol*CODATA08_AU_TO_KCAL/CODATA08_A_TO_BOHRS
             write(6,'(1x,"ML/MM: QM/MM ELECTROSTATIC+POLARIZATION ENERGY=",f14.4)') eqmmmtorchani*CODATA08_AU_TO_KCAL/CODATA08_A_TO_BOHRS
         else
             write(6,'(1x,"ML/MM: QM/MM ELECTROSTATIC ENERGY =",f14.4)') eqmmmtorchani*CODATA08_AU_TO_KCAL/CODATA08_A_TO_BOHRS
         end if
      end if

      eqmmmtorchani = eqmmmtorchani*CODATA08_AU_TO_KCAL/CODATA08_A_TO_BOHRS

      escf = escf + eqmmmtorchani
      ! Calculates QM/MM coumlomb interaction forces
      call get_mlmm_forces_numerical(nqmatoms,nclatoms,dxyzqm,qmcoords,qmtypes,clcoords,dxyzcl,dxyzqm_qmmm, &
                                  qmcharges,dqmcharges,polarize_qm_charges,distort_qm_energy,&
                                  distij,efield,alpha,use_charges_derivatives,&
                                  qm_atomic_symbols,mlmm_patch_path,psi4_path,python_executable, charges_model, charges_functional,&
                                  charges_basis, nthreads)


    end subroutine get_mlmm_ener_forces

    subroutine get_distij_matrix(nqmatoms,nclatoms,qmcoords,clcoords,distij)

      implicit none

      integer, intent(in) :: nqmatoms                       ! Number of QM atoms
      double precision,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
      integer, intent(in) :: nclatoms                       ! Number of MM atoms
      double precision,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
      double precision, intent(out) :: distij(nqmatoms,nclatoms) ! QM atom eqmmm forces
      integer i,j

      do i=1,nqmatoms
        do j=1,nclatoms
          distij(i,j) = (qmcoords(1,i)-clcoords(1,j))**2
          distij(i,j) = distij(i,j) + (qmcoords(2,i)-clcoords(2,j))**2
          distij(i,j) = distij(i,j) + (qmcoords(3,i)-clcoords(3,j))**2
          distij(i,j) = DSQRT(distij(i,j))
        end do
      end do

    end subroutine

    subroutine get_elec_field(nqmatoms,nclatoms,qmcoords,clcoords,distij,efield)

      implicit none

      integer, intent(in) :: nqmatoms                       ! Number of QM atoms
      double precision,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
      integer, intent(in) :: nclatoms                       ! Number of MM atoms
      double precision,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
      double precision, intent(out):: efield(3,nqmatoms)      ! Electric field on QM atoms
      double precision, intent(in) :: distij(nqmatoms,nclatoms) ! QM/MM distance matrix
      double precision efieldvec(3)                        ! atomic electric field vector
      integer i, j

      do i=1,nqmatoms
        efieldvec = 0.d0
        do j=1,nclatoms
          efieldvec(1)=efieldvec(1)+clcoords(4,j)*(clcoords(1,j)-qmcoords(1,i))/(distij(i,j)**3)
          efieldvec(2)=efieldvec(2)+clcoords(4,j)*(clcoords(2,j)-qmcoords(2,i))/(distij(i,j)**3)
          efieldvec(3)=efieldvec(3)+clcoords(4,j)*(clcoords(3,j)-qmcoords(3,i))/(distij(i,j)**3)
        end do
        efield(1,i)=efieldvec(1)
        efield(2,i)=efieldvec(2)
        efield(3,i)=efieldvec(3)
      end do

      end subroutine get_elec_field

      subroutine get_atomic_polarizability(nqmatoms,qmtypes,alpha,pol_H,pol_He,pol_Li,pol_Be,pol_B,&
                 pol_C,pol_N,pol_O,pol_F,pol_Ne,pol_Na,pol_Mg,pol_Al,pol_Si,pol_P,pol_S,pol_Cl,pol_Ar,first_step)

        implicit none

        integer, intent(in) :: nqmatoms  ! Number of QM atoms
        integer, intent(in) :: qmtypes(nqmatoms) ! QM atoms atomic numbers
        double precision, intent(in) :: pol_H
        double precision, intent(in) :: pol_He
        double precision, intent(in) :: pol_Li
        double precision, intent(in) :: pol_Be
        double precision, intent(in) :: pol_B
        double precision, intent(in) :: pol_C
        double precision, intent(in) :: pol_N
        double precision, intent(in) :: pol_O
        double precision, intent(in) :: pol_F
        double precision, intent(in) :: pol_Ne
        double precision, intent(in) :: pol_Na
        double precision, intent(in) :: pol_Mg
        double precision, intent(in) :: pol_Al
        double precision, intent(in) :: pol_Si
        double precision, intent(in) :: pol_P
        double precision, intent(in) :: pol_S
        double precision, intent(in) :: pol_Cl
        double precision, intent(in) :: pol_Ar
        double precision, intent(out) :: alpha(nqmatoms) ! QM atoms atomic polarizabilities
        double precision alphatable(18) ! atomic polarizabilities for atoms from H to Ar
        character(len=2), dimension(18) :: atomic_symbols
        logical, intent(in) :: first_step
        logical element_in_system
        integer i,j
        ! extracted from Litman, J. M., Liu, C., & Ren, P. (2021). Journal of Chemical Information and Modeling, 62(1), 79-87.
        ! corresponds to the average of all hydrogen types
        alphatable(1)=3.08d0
        ! extracted from Schwerdtfeger, P., & Nagle, J. K. (2019). Molecular Physics, 117(9-12), 1200-1225.
        alphatable(2)=1.38375
        alphatable(3)=164.1125
        alphatable(4)=37.74
        alphatable(5)=20.5
        alphatable(6)=11.3
        alphatable(7)=7.4
        alphatable(8)=5.3
        alphatable(9)=3.74
        alphatable(10)=2.6611
        alphatable(11)=162.7
        alphatable(12)=71.2
        alphatable(13)=57.8
        alphatable(14)=37.3
        alphatable(15)=25
        alphatable(16)=19.4
        alphatable(17)=14.6
        alphatable(18)=11.083

        if (pol_H  .ge. 0.d0) alphatable(1)=pol_H
        if (pol_He .gt. 0.d0) alphatable(2)=pol_He
        if (pol_Li .gt. 0.d0) alphatable(3)=pol_Li
        if (pol_Be .gt. 0.d0) alphatable(4)=pol_Be
        if (pol_B  .gt. 0.d0) alphatable(5)=pol_B
        if (pol_C  .gt. 0.d0) alphatable(6)=pol_C
        if (pol_N  .gt. 0.d0) alphatable(7)=pol_N
        if (pol_O  .gt. 0.d0) alphatable(8)=pol_O
        if (pol_F  .gt. 0.d0) alphatable(9)=pol_F
        if (pol_Ne .gt. 0.d0) alphatable(10)=pol_Ne
        if (pol_Na .gt. 0.d0) alphatable(11)=pol_Na
        if (pol_Mg .gt. 0.d0) alphatable(12)=pol_Mg
        if (pol_Al .gt. 0.d0) alphatable(13)=pol_Al
        if (pol_Si .gt. 0.d0) alphatable(14)=pol_Si
        if (pol_P  .gt. 0.d0) alphatable(15)=pol_P
        if (pol_S  .gt. 0.d0) alphatable(16)=pol_S
        if (pol_Cl .gt. 0.d0) alphatable(17)=pol_Cl
        if (pol_Ar .gt. 0.d0) alphatable(18)=pol_Ar

        do i=1,nqmatoms
          alpha(i)=alphatable(qmtypes(i))
        end do

        if(first_step) then
          atomic_symbols(1) = 'H '
          atomic_symbols(2) = 'He'
          atomic_symbols(3) = 'Li'
          atomic_symbols(4) = 'Be'
          atomic_symbols(5) = 'B'
          atomic_symbols(6) = 'C'
          atomic_symbols(7) = 'N'
          atomic_symbols(8) = 'O'
          atomic_symbols(9) = 'F'
          atomic_symbols(10) = 'Ne'
          atomic_symbols(11) = 'Na'
          atomic_symbols(12) = 'Mg'
          atomic_symbols(13) = 'Al'
          atomic_symbols(14) = 'Si'
          atomic_symbols(15) = 'P'
          atomic_symbols(16) = 'S'
          atomic_symbols(17) = 'Cl'
          atomic_symbols(18) = 'Ar'
          write(6,*) "ML/MM: The following atomic polarizabilities (alpha) will be used (A.U.)"
          do i=1,18
            element_in_system=.false.
            do j=1,nqmatoms
              if(qmtypes(j)==i) element_in_system=.true.
            end do
              if(element_in_system) then
                 write(6,'(A, A, A, F10.4)') " ML/MM: alpha(", trim(atomic_symbols(i)), ") = ", alphatable(i)
              end if
          end do

        endif
        !Convertion from A.U. to Angstroms**3
        alpha=alpha*0.1481847d0
      end subroutine get_atomic_polarizability


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Calculates the QMMM forces in the QM and the MM regions numerically
    ! ----------------------------------------------------------------
    subroutine get_mlmm_forces_numerical(nqmatoms,nclatoms,dxyzqm,qmcoords,qmtypes,clcoords,dxyzcl,dxyzqm_qmmm, &
                                   qmcharges,dqmcharges,polarize_qm_charges,distort_qm_energy, &
                                   distij, efield,alpha,use_charges_derivatives,qm_atomic_symbols,mlmm_patch_path,&
                                   psi4_path,python_executable, charges_model, charges_functional,charges_basis, nthreads)

      use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS

      implicit none

      integer, intent(in) :: nqmatoms                       ! Number of QM atoms
      double precision,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
      integer, intent(in) :: nclatoms                       ! Number of MM atoms
      double precision,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
      double precision, intent(in):: qmcharges(nqmatoms)    ! QM atoms charges
      double precision, intent(in) :: dqmcharges(3,nqmatoms) ! QM atoms charges forces
      integer, intent(in) :: qmtypes(nqmatoms)               ! QM atoms atomic numbers
      double precision, intent(inout) :: dxyzqm(3,nqmatoms) ! QM atom forces
      double precision, intent(inout) :: dxyzcl(3,nclatoms) ! MM atom qmmm forces
      double precision, intent(out) :: dxyzqm_qmmm(3,nqmatoms) ! QM atom eqmmm forces
      logical, intent(in) :: polarize_qm_charges ! QMMM electrostatic contribution
      logical, intent(in) :: use_charges_derivatives ! Use derivatives of torchani charges
      logical, intent(in) :: distort_qm_energy ! Include distortion correction of QM energy
      double precision distortion_k ! Distortion constant
      double precision epol, edist ! Polarization and distortion energies
      double precision pert_escf, eqmmmtorchani_plusdelta, eqmmmtorchani_minusdelta
      double precision eqmmmtorchani
      double precision distij(nqmatoms,nclatoms) ! QM/MM distance matrix
      double precision, intent(in) :: alpha(nqmatoms) ! QM atoms atomic polarizabilities
      double precision efield(3,nqmatoms) ! Electric field on QM atoms
      double precision :: pert_qmcoords(3,nqmatoms) ! perturbed QM atom coordinates
      double precision :: pert_clcoords(4,nclatoms) ! perturbed MM atom coordinates and unperturbed charges in au
      double precision :: pert_qmcharges(nqmatoms)    ! perturbed QM atoms charges
      double precision :: pert_dxyzqm(3,nqmatoms) ! perturbed QM atom eqmmm forces
      double precision :: force ! force on a specific coordinate
      double precision :: delta ! coordinates perturbation factor
      character(len=2), intent(in) :: qm_atomic_symbols(nqmatoms) ! QM atoms atomic numbers
      character(len=256), intent(in) ::python_executable ! Python executable
      character(len=256), intent(in) ::mlmm_patch_path ! Path to mlmm_patch
      character(len=256), intent(in) ::psi4_path ! Path to psi4
      character(len=256), intent(in) ::charges_model ! Path to psi4
      character(len=256), intent(in) ::charges_functional ! psi4 functional
      character(len=256), intent(in) ::charges_basis ! psi4 basis
      integer, intent(in) :: nthreads ! psi4 nthreads
      integer k, m, i, j, unit_num

      distortion_k=0.5d0
      dxyzcl=0.d0
      dxyzqm_qmmm=dxyzqm  !saves initial forces in qm region
      pert_qmcharges=qmcharges

      delta=0.005
      do k=1,nqmatoms
        do m=1,3
          ! Summs delta to a temporary (pert_qmcoords)  array, in dimension m
          pert_qmcoords(:,:)=qmcoords(:,:)
          pert_qmcoords(m,k)=pert_qmcoords(m,k)+delta

          ! If necesary, calculates perturbed qm charges (and other things that we don't actually need...)
          if ((trim(charges_model)/='AMBER_CHARGES') .and. use_charges_derivatives) then
             call get_psi4_charges(pert_qmcoords,nqmatoms,pert_qmcharges,qm_atomic_symbols,pert_escf,&
                                 mlmm_patch_path,psi4_path,python_executable, charges_model, charges_functional, charges_basis, nthreads)
          endif

          ! Calculates the QM/MM energy
          ! Calculates QM/MM distance matrix
          call get_distij_matrix(nqmatoms,nclatoms,pert_qmcoords,clcoords,distij)

          ! Calculates QM/MM coumlomb interaction energy
          eqmmmtorchani=0.d0
          do i=1,nqmatoms
            do j=1,nclatoms
               eqmmmtorchani = eqmmmtorchani + pert_qmcharges(i)*clcoords(4,j)/distij(i,j)
            end do
          end do

          ! Calculates QM/MM polarization energy
          epol=0.d0
          ! Only polarization on the QM atoms is accounted
          if (polarize_qm_charges) then
            call get_elec_field(nqmatoms,nclatoms,pert_qmcoords,clcoords,distij,efield)
            do i=1,nqmatoms
              epol=epol+alpha(i)*(efield(1,i)**2+efield(2,i)**2+efield(3,i)**2)
            end do
            epol=-0.5d0*epol
          end if

          ! Calculates QM distortion energy
          edist=0.d0
          if (distort_qm_energy) then
            edist=-distortion_k*epol
          end if

          eqmmmtorchani_plusdelta = eqmmmtorchani + epol + edist
          eqmmmtorchani_plusdelta = eqmmmtorchani_plusdelta*CODATA08_AU_TO_KCAL/CODATA08_A_TO_BOHRS !kcal/mol

          ! Repeats but now substracting delta
          pert_qmcoords(:,:)=qmcoords(:,:)
          pert_qmcoords(m,k)=pert_qmcoords(m,k)-delta

          ! Calculates perturbed qm charges (and other things that we don't actually need...)
          if ((trim(charges_model)/='AMBER_CHARGES') .and. use_charges_derivatives) then
             call get_psi4_charges(pert_qmcoords,nqmatoms,pert_qmcharges,qm_atomic_symbols,pert_escf,&
                                  mlmm_patch_path,psi4_path,python_executable, charges_model, charges_functional, charges_basis, nthreads)
          endif

          ! Calculates QM/MM distance matrix
          call get_distij_matrix(nqmatoms,nclatoms,pert_qmcoords,clcoords,distij)

          ! Calculates QM/MM coumlomb interaction energy
          eqmmmtorchani=0.d0
          do i=1,nqmatoms
            do j=1,nclatoms
              eqmmmtorchani = eqmmmtorchani + pert_qmcharges(i)*clcoords(4,j)/distij(i,j)
            end do
          end do

          ! Calculates QM/MM polarization energy
          epol=0.d0
          ! Only polarization on the QM atoms is accounted
          if (polarize_qm_charges) then
            call get_elec_field(nqmatoms,nclatoms,pert_qmcoords,clcoords,distij,efield)
            do i=1,nqmatoms
              epol=epol+alpha(i)*(efield(1,i)**2+efield(2,i)**2+efield(3,i)**2)
            end do
            epol=-0.5d0*epol
          end if

          ! Calculates QM distortion energy
          edist=0.d0
          if (distort_qm_energy) then
            edist=-distortion_k*epol
          end if

          eqmmmtorchani_minusdelta = eqmmmtorchani + epol + edist
          eqmmmtorchani_minusdelta = eqmmmtorchani_minusdelta*CODATA08_AU_TO_KCAL/CODATA08_A_TO_BOHRS !kcal/mol

          force = (eqmmmtorchani_plusdelta - eqmmmtorchani_minusdelta)/(2.d0*delta) !kcal/(mol Bohr)
          dxyzqm(m,k) = dxyzqm(m,k) + force
        end do
     end do

     ! Does the same but for the MM atoms

     do k=1,nclatoms
       do m=1,3
         ! Summs delta to a temporary (pert_qmcoords)  array, in dimension m
         pert_clcoords(:,:)=clcoords(:,:)
         pert_clcoords(m,k)=pert_clcoords(m,k)+delta

         ! Calculates the QM/MM energy
         ! Calculates QM/MM distance matrix
         call get_distij_matrix(nqmatoms,nclatoms,qmcoords,pert_clcoords,distij)

         ! Calculates QM/MM coumlomb interaction energy
         eqmmmtorchani=0.d0
         do i=1,nqmatoms
           do j=1,nclatoms
              eqmmmtorchani = eqmmmtorchani + qmcharges(i)*pert_clcoords(4,j)/distij(i,j)
           end do
         end do

         ! Calculates QM/MM polarization energy
         epol=0.d0
         ! Only polarization on the QM atoms is accounted
         if (polarize_qm_charges) then
           call get_elec_field(nqmatoms,nclatoms,qmcoords,pert_clcoords,distij,efield)
           do i=1,nqmatoms
             epol=epol+alpha(i)*(efield(1,i)**2+efield(2,i)**2+efield(3,i)**2)
           end do
           epol=-0.5d0*epol
         end if

         ! Calculates QM distortion energy
         edist=0.d0
         if (distort_qm_energy) then
           edist=-distortion_k*epol
         end if

         eqmmmtorchani_plusdelta = eqmmmtorchani + epol + edist
         eqmmmtorchani_plusdelta = eqmmmtorchani_plusdelta*CODATA08_AU_TO_KCAL/CODATA08_A_TO_BOHRS !kcal/mol

         ! Repeats but now substracting delta
         pert_clcoords(:,:)=clcoords(:,:)
         pert_clcoords(m,k)=pert_clcoords(m,k)-delta

         ! Calculates QM/MM distance matrix
         call get_distij_matrix(nqmatoms,nclatoms,qmcoords,pert_clcoords,distij)

         ! Calculates QM/MM coumlomb interaction energy
         eqmmmtorchani=0.d0
         do i=1,nqmatoms
           do j=1,nclatoms
             eqmmmtorchani = eqmmmtorchani + qmcharges(i)*pert_clcoords(4,j)/distij(i,j)
           end do
         end do

         ! Calculates QM/MM polarization energy
         epol=0.d0
         ! Only polarization on the QM atoms is accounted
         if (polarize_qm_charges) then
           call get_elec_field(nqmatoms,nclatoms,qmcoords,pert_clcoords,distij,efield)
           do i=1,nqmatoms
             epol=epol+alpha(i)*(efield(1,i)**2+efield(2,i)**2+efield(3,i)**2)
           end do
           epol=-0.5d0*epol
         end if

         ! Calculates QM distortion energy
         edist=0.d0
         if (distort_qm_energy) then
           edist=-distortion_k*epol
         end if

         eqmmmtorchani_minusdelta = eqmmmtorchani + epol + edist
         eqmmmtorchani_minusdelta = eqmmmtorchani_minusdelta*CODATA08_AU_TO_KCAL/CODATA08_A_TO_BOHRS !kcal/mol

         force = (eqmmmtorchani_plusdelta - eqmmmtorchani_minusdelta)/(2.d0*delta)
         dxyzcl(m,k) = dxyzcl(m,k) + force
       end do
     end do

     dxyzqm_qmmm=dxyzqm-dxyzqm_qmmm  !saves the qmmm forces in the qm region

     end subroutine get_mlmm_forces_numerical

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
     ! Gets the atomic symbols for each element in the QM subsystem
     subroutine get_atomic_symbols(nqmatoms,iqm_atomic_numbers,qm_atomic_symbols)

       implicit none

       integer, intent(in) :: nqmatoms                         ! Number of QM atoms
       integer, intent(in) :: iqm_atomic_numbers(nqmatoms)     ! QM atoms atomic numbers
       character(len=2), intent(out) :: qm_atomic_symbols(nqmatoms) ! QM atoms atomic numbers
       character(len=2), dimension(18) :: atomic_symbols    ! Atomic symbols list
       integer i

       atomic_symbols(1) = 'H '
       atomic_symbols(2) = 'He'
       atomic_symbols(3) = 'Li'
       atomic_symbols(4) = 'Be'
       atomic_symbols(5) = 'B'
       atomic_symbols(6) = 'C'
       atomic_symbols(7) = 'N'
       atomic_symbols(8) = 'O'
       atomic_symbols(9) = 'F'
       atomic_symbols(10) = 'Ne'
       atomic_symbols(11) = 'Na'
       atomic_symbols(12) = 'Mg'
       atomic_symbols(13) = 'Al'
       atomic_symbols(14) = 'Si'
       atomic_symbols(15) = 'P'
       atomic_symbols(16) = 'S'
       atomic_symbols(17) = 'Cl'
       atomic_symbols(18) = 'Ar'

       do i=1,nqmatoms
          qm_atomic_symbols(i)=atomic_symbols(iqm_atomic_numbers(i))
       end do

    end subroutine get_atomic_symbols

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Calculates energy and forces for the QM (ML) region
    subroutine  get_ml_properties(qmcoords,nqmatoms,qm_atomic_symbols,dxyzqm,escf,mlmm_patch_path,python_executable)

       implicit none

       double precision,  intent(in) :: qmcoords(3,nqmatoms)  ! QM atom coordinates
       integer, intent(in) :: nqmatoms                        ! Number of QM atoms
       character(len=2), intent(in) :: qm_atomic_symbols(nqmatoms) ! QM atoms atomic numbers
       character(len=256), intent(in) ::mlmm_patch_path ! Path to mlmm_patch
       character(len=256), intent(in) ::python_executable ! Python executable
       character(len=256) :: read_buffer, call_buffer, line, filename
       double precision, intent(out) :: dxyzqm(3,nqmatoms)! QM atom charges gradient
       double precision, intent(out) :: escf! QM energy
       integer i, ios
       !Write input
       filename='input_ml.tmp'
       open(unit=10, file=trim(filename), status='replace', action='write')
       ! Write the header
       write(10, "('molecule mol {')")
       ! Write coordinates and symbols
       do i = 1, nqmatoms
          write(10, "('',A2,3F15.8)") qm_atomic_symbols(i), qmcoords(1, i), qmcoords(2, i), qmcoords(3, i)
       end do
       ! Write the footer
       write(10, "('}')")
       ! Close the file
       close(unit=10)

       ! Run mlmm_interface (torchani)
       call_buffer = trim(python_executable)//" "//trim(mlmm_patch_path)//'/'//"get_ml_properties.py"
       call system(trim(call_buffer)) ! This generates f_ml.tmp e_ml files

       ! Read outputs
       open (unit=20, file='f_ml.tmp', status='old', action='read')
       i=1
       do while (i .le. nqmatoms)
          read(20,*) dxyzqm(1,i), dxyzqm(2,i), dxyzqm(3,i)
       i = i + 1
       end do
       close(20)
       ! Read outputs
       open (unit=30, file='e_ml.tmp', status='old', action='read')
       read(30,*) escf
       close(30)

       ! Organize the folder
       call_buffer = " mv f_ml.tmp f_ml.tmp.old"
       call system(trim(call_buffer))
       call_buffer = " mv e_ml.tmp e_ml.tmp.old"
       call system(trim(call_buffer))
       call_buffer = " mv input_ml.tmp input_ml.tmp.old"
       call system(trim(call_buffer))

    end subroutine get_ml_properties


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Calculates charges for the QM (ML) region
    subroutine  get_psi4_charges(qmcoords,nqmatoms,qmcharges,qm_atomic_symbols,escf,mlmm_patch_path,psi4_path,&
                                python_executable, charges_model, charges_functional, charges_basis, nthreads)

       implicit none

       double precision,  intent(in) :: qmcoords(3,nqmatoms)  ! QM atom coordinates
       double precision,  intent(inout) :: qmcharges(nqmatoms)  ! QM atom charges
       integer, intent(in) :: nqmatoms                        ! Number of QM atoms
       character(len=2), intent(in) :: qm_atomic_symbols(nqmatoms) ! QM atoms atomic symbols
       character(len=256), intent(in) ::mlmm_patch_path ! Path to mlmm_patch
       character(len=256), intent(in) ::python_executable ! Python executable
       character(len=256), intent(in) ::psi4_path ! Path to psi4
       character(len=256), intent(in) ::charges_model !
       character(len=256), intent(in) ::charges_functional ! psi4 functional
       character(len=256), intent(in) ::charges_basis ! psi4 basis
       integer, intent(in) :: nthreads ! psi4 nthreads
       character(len=256) :: read_buffer, chnthreads, call_buffer,line,filename
       double precision, intent(out) :: escf! QM energy
       integer i, ios

       !Write input
       filename='input_psi4.tmp'
       open(unit=10, file=trim(filename), status='replace', action='write')
       ! Write the header
       write(10, "('molecule mol {')")
       ! Write coordinates and symbols
       do i = 1, nqmatoms
          write(10, "('',A2,3F15.8)") qm_atomic_symbols(i), qmcoords(1, i), qmcoords(2, i), qmcoords(3, i)
       end do
       ! Write the footer
       write(10, "('}')")
       write(10,*)
       write(10, "('set {')")
       line="basis "//trim(charges_basis)
       write(10,"('',A)") trim(line)
       write(10, "('}')")
       write(10,*)
       line="e, wfn = energy('"//trim(charges_functional)//"', return_wfn=True)"
       write(10,"('',A)") trim(line)
       line="oeprop(wfn, '"//trim(charges_model)//"')"
       write(10,"('',A)") trim(line)
       ! Close the file
       close(unit=10)

       if (nthreads .le. 9) write(chnthreads,'(I1)') nthreads
       if (nthreads .gt. 9 .and. nthreads .le. 99) write(chnthreads,'(I2)') nthreads
       if (nthreads .gt. 99) write(chnthreads,'(I3)') nthreads

       ! Run mlmm_interface (torchani)
       call_buffer = trim(python_executable)//" "//trim(mlmm_patch_path)//'/'//"get_psi4_charges.py "//trim(chnthreads)//" "//trim(psi4_path)
       call system(trim(call_buffer)) ! This generates q_psi4.tmp e_psi4.tmp files

      ! Read outputs
      open (unit=20, file='q_psi4.tmp', status='old', action='read')
      i=1
      do while (i .le. nqmatoms)
         read(20,*) qmcharges(i)
      i = i + 1
      end do
      close(20)
      ! Read outputs
      open (unit=30, file='e_psi4.tmp', status='old', action='read')
      read(30,*) escf
      close(30)

      ! Organize the folder
      call_buffer = " mv e_psi4.tmp e_psi4.tmp.old"
      call system(trim(call_buffer))
      call_buffer = " mv q_psi4.tmp q_psi4.tmp.old"
      call system(trim(call_buffer))
      call_buffer = " mv input_psi4.tmp input_psi4.tmp.old"
      call system(trim(call_buffer))
      call_buffer = " mv output_psi4.tmp output_psi4.tmp.old"
      call system(trim(call_buffer))
      call_buffer = " mv output_psi4.log output_psi4.log.old"
      call system(trim(call_buffer))
       
    end subroutine get_psi4_charges

end module qm2_extern_mlmm_module
