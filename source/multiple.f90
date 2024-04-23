subroutine multiple
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Multiple emission
!
! Author    : Arjan Koning and Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
! All global variables
!   numen           ! maximum number of outgoing energies
!   numex           ! maximum number of excitation energies
!   numJ            ! maximum J - value
!   numNchan        ! maximum number of outgoing neutron in individual channel description
!   numpar          ! number of particles
!   numZchan        ! maximum number of outgoing protons in individual channel description
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxNrp          ! maximal number of neutrons away from the initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
!   maxZrp          ! maximal number of protons away from the initial compound nucleus
!   nanglecont      ! number of angles for continuum
!   popeps          ! limit for population cross sections
! Variables for basic reaction
!   flagchannels    ! flag for exclusive channels calculation
!   flagendf        ! flag for information for ENDF - 6 file
!   flagrecoil      ! flag for calculation of recoils
!   flagrpevap      ! flag for evaporation of residual products at high inccident energies
! Variables for input energies
!   flaginitpop     ! flag for initial population distribution
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ntarget         ! neutron number of target nucleus
!   Ztarget         ! charge number of target nucleus
! Variables for output
!   flagcompo       ! flag for output of cross section components
!   flagbinspec     ! flag for output of emission spectrum per excitation bin
!   flagcheck       ! flag for output of numerical checks
!   flagdecay       ! flag for output of decay of each population bin
!   flagddx         ! flag for output of double - differential cross sections
!   flagpop         ! flag for output of population
!   flagspec        ! flag for output of spectra
! Variables for compound reactions
!   adjustTJ        ! logical for energy-dependent TJ adjustment
!   flagcomp        ! flag for compound angular distribution calculation
!   skipCN          ! flag to skip compound nucleus in evaporation chain
! Variables for preequilibrium
!   flag2comp       ! flag for two - component pre - equilibrium model
! Variables for fission
!   flagfisfeed     ! flag for output of fission per excitation bin
!   flagfisout      ! flag for output of fission information
!   flagfission     ! flag for fission
! Variables for discrete levels
!   flaglevels      ! flag for output of discrete level information
! Variables for level density
!   filedensity     ! flag for level densities on separate files
!   flagdensity     ! flag for output of level densities
! Variables for OMP
!   flagompall      ! flag for new optical model calculation for all residual
!   flagomponly     ! flag to execute ONLY an optical model calculation
! Variables for gamma rays
!   fisominit       ! correction factor for isospin forbidden transitions for multiple emission
! Variables for gamma rays
!   fisom           ! correction factor for isospin forbidden transitions for multiple emission
! Variables for energy grid
!   anglecont       ! angle in degrees for continuum
!   deltaE          ! energy bin around outgoing energies
!   ebegin          ! first energy point of energy grid
!   egrid           ! outgoing energy grid
!   Einc            ! incident energy in MeV
! Variables for energies
!   eend            ! last energy point of energy grid
!   eendhigh        ! last energy point for energy grid for any particle
!   mulpreZN        ! logical for multiple pre - equilibrium per nucleus
! Variables for excitation energy grid
!   deltaEx         ! excitation energy bin for population arrays
!   Ex              ! excitation energy
!   Exmax           ! maximum excitation energy for residual nucleus
!   fisfeedJP       ! fission contribution from excitation energy bin per J, P
!   maxex           ! maximum excitation energy bin for residual nucleus
!   maxJ            ! maximal J - value
!   nexmax          ! maximum excitation energy bin for residual nucleus
! Variables for multiple emission
!   Dmulti          ! depletion factor for multiple preequilibrium
!   Eaveragemul     ! average outgoing energy
!   Fcomp           ! compound population fraction per nucleus
!   Fdir            ! direct population fraction per nucleus
!   feedexcl        ! feeding terms from compound excitation ene
!   fisfeedex       ! fission contribution from excitation energy bin
!   Fpreeq          ! preequilibrium population fraction per nucleus
!   mcontrib        ! contribution to emission spectrum
!   mpecontrib      ! contribution to multiple pre - equilibrium emission spectr
!   popexcl         ! population cross section of bin just before decay
!   xsbinspec       ! emission spectrum from compound nucleus per bin
!   xsfeed          ! cross section from compound to residual nucleus
!   xsmpe           ! multiple - preequilibrium cross section per energy bin
!   xsmpeemis       ! multiple - preequilibrium emission spectrum from compound n
!   xsmpetot        ! total multiple - preequilibrium cross section
!   xsmpreeq        ! multiple pre - equilibrium emission spectrum
!   xsmpreeqad      ! multiple preequilibrium angular distribution
!   xsngn           ! total (projectile, gamma - ejectile) cross section
!   xsngnspec       ! total (projectile, gamma - ejectile) spectrum
!   xspartial       ! emitted cross section flux per energy bin
!   xspopcomp       ! compound population cross section per nucleus
!   xspopnuc0       ! population cross section per nucleus
!   xspoppreeq      ! preequilibrium population cross section per nucleus
! Variables for binary reactions
!   feedbinary      ! feeding from first compound nucleus
!   xspopdir        ! direct population cross section per nucleus
! Variables for binary emission spectra
!   xscomp          ! compound elastic cross section
!   xscompad        ! compound emission angular distribution
!   xsemis          ! cross section for emission from compound nucleus
! Variables for compound nucleus from target
!   dExinc          ! excitation energy bin for mother nucleus
!   Exinc           ! excitation energy of entrance bin
!   Fnorm           ! multiplication factor
! Variables for incident channel
!   partdecay       ! total decay per particle
!   popdecay        ! decay from population
!   preeqpopex      ! pre - equilibrium population
!   xsbinary        ! cross section from initial compound to residual nucleus
!   xsngnsum        ! sum over total (projectile, gamma - ejectile) cross section
!   xspop           ! population cross section
!   xspopex         ! population cross section summed over spin and parity
!   xspopexP        ! population cross section per parity
!   xspopnuc        ! population cross section per nucleus
!   xspopnucP       ! population cross section per nucleus per parity
!   xsreacinc       ! reaction cross section for incident channel
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   Nindex          ! neutron number index for residual nucleus
!   NN              ! neutron number of residual nucleus
!   parinclude      ! logical to include outgoing particle
!   parskip         ! logical to skip outgoing particle
!   primary         ! flag to designate primary (binary) reaction
!   strucexist      ! flag to state whether structure info for this nucleus exists
!   strucwrite      ! flag for output of nuclear structure info
!   Zindex          ! charge number index for residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   deg2rad         ! conversion factor for degrees to radians
!   fourpi          ! 4. * pi
!   nuc             ! symbol of nucleus
!   parname         ! name of particle
! Variables for levels
!   tau             ! lifetime of state in seconds
! Variables for fission parameters
!   nfisbar         ! number of fission barrier parameters
! Variables for level density
!   Nlast           ! last discrete level
! Variables for masses
!   S               ! separation energy
! Variables for preequilibrium
!   xspopph         ! population cross section p
!   xspopph2        ! population cross section p
!
! *** Declaration of local data
!
  implicit none
  character(len=13)  :: fisfile               ! fission file
  character(len=13)  :: rpfile                ! file with residual production cross sections
  character(len=60)  :: form1                 ! format
  character(len=60)  :: form2                 ! format
  character(len=132) :: key                   ! keyword
  integer            :: A                     ! mass number of target nucleus
  integer            :: Ares                  ! mass number of residual nucleus
  integer            :: h                     ! help variable
  integer            :: iang                  ! running variable for angle
  integer            :: idensfis              ! identifier to activate possible fission level densities
  integer            :: Ir                    ! spin 
  integer            :: J                     ! spin of level
  integer            :: J2                    ! 2 * J
  integer            :: Jfis                  ! spin of fissioning system
  integer            :: Jres                  ! spin of residual nucleus
  integer            :: N                     ! neutron number of residual nucleus
  integer            :: Ncomp                 ! neutron number index for compound nucleus
  integer            :: nen                   ! energy counter
  integer            :: nex                   ! excitation energy bin of compound nucleus
  integer            :: nexout                ! energy index for outgoing energy
  integer            :: Nix                   ! neutron number index for residual nucleus
  integer            :: NL                    ! last discrete level
  integer            :: Nres                  ! neutron number of residual nucleus
  integer            :: odd                   ! odd (1) or even (0) nucleus
  integer            :: oddres                ! odd (1) or even (0) residual nucleus
  integer            :: p                     ! particle number
  integer            :: Pres                  ! parity of residual nucleus
  integer            :: parity                ! parity
  integer            :: type                  ! particle type
  integer            :: Z                     ! charge number of target nucleus
  integer            :: Zcomp                 ! proton number index for compound nucleus
  integer            :: Zix                   ! charge number index for residual nucleus
  integer            :: Zres                  ! charge number of residual nucleus
  real(sgl)          :: ang                   ! angle
  real(sgl)          :: dEx                   ! excitation energy bin for population arrays
  real(sgl)          :: Eaveragesum           ! help variable
  real(sgl)          :: emissum(0:numpar)     ! integrated binary emission spectrum
  real(sgl)          :: Eout                  ! outgoing energy
  real(sgl)          :: Exm                   ! maximal attainable energy
  real(sgl)          :: Exmin                 ! help variable
  real(sgl)          :: factor                ! help variable
  real(sgl)          :: kalbach               ! Kalbach function
  real(sgl)          :: popepsA               ! limit for population cross sections per energy
  real(sgl)          :: popepsB               ! limit for population cross sections per spin and parity
  real(sgl)          :: Smax                  ! separation energy
  real(sgl)          :: Smin                  ! minimal separation energy
  real(sgl)          :: sumxs                 ! sum over emission channels
  real(sgl)          :: xsdif                 ! difference between in-flux and out-flux per bin
  real(sgl)          :: xsmax                 ! maximum cross sections
  real(sgl)          :: xsp                   ! help variable
  real(sgl)          :: rJ                    ! help variable
  real(sgl)          :: xspopsave(0:numex)    ! help variable for diagnosis
!
! ******************** Loop over nuclei ********************************
!
! Loop over all residual nuclei, starting with the initial compound nucleus (Zcomp=0, Ncomp=0),
! and then according to decreasing Z and N.
!
! excitation : subroutine for excitation energy population
!
  primary = .false.
  if (flagomponly) return
  if (flaginitpop) call excitation
  if (flagpop) write(*, '(/" ########## MULTIPLE EMISSION ##########")')
  do Zcomp = 0, maxZ
    do Ncomp = 0, maxN
!
! Continue for this compound nucleus only when it contains sufficient reaction flux.
! If so, determine the nuclear structure for the residual nuclei.
! Initialize emission spectrum for this compound nucleus.
! Determine excitation energy grid for residual nuclei.
!
! structure : subroutine for nuclear structure parameters
! exgrid    : subroutine to set excitation energy grid
!
      if (skipCN(Zcomp, Ncomp) == 1) cycle
      if (xspopnuc(Zcomp, Ncomp) < popeps) then
        xspopnuc(Zcomp, Ncomp) = 0.
        goto 500
      endif
      Z = ZZ(Zcomp, Ncomp, 0)
      N = NN(Zcomp, Ncomp, 0)
      A = AA(Zcomp, Ncomp, 0)
      if (flagrpevap .and. (Zcomp == maxZrp .or. Ncomp == maxNrp)) then
        xspopnuc0(Z, A) = xspopnuc(Zcomp, Ncomp)
        rpfile = 'rp000000.ex'
        write(rpfile(3:5), '(i3.3)') Z
        write(rpfile(6:8), '(i3.3)') A
        open (unit = 1, file = rpfile, status = 'replace')
        write(1, * ) maxex(Zcomp, Ncomp) + 1, 30, 1, " xs= ", xspopnuc(Zcomp, Ncomp)
        do nex = 0, maxex(Zcomp, Ncomp)
          do parity = - 1, 1, 2
            write(1, '(f10.5, 31es12.5)') Ex(Zcomp, Ncomp, nex), (xspop(Zcomp, Ncomp, nex, J, parity), J = 0, 30)
          enddo
        enddo
        close(1)
        cycle
      endif
      do type = 0, 6
        if (flagspec) then
          do nen = 0, numen
            xsemis(type, nen) = 0.
            xsmpeemis(type, nen) = 0.
            do nex = 0, numex + 1
              xsbinspec(type, nex, nen) = 0.
            enddo
          enddo
        endif
        do nex = 0, numex + 1
          xspartial(type, nex) = 0.
          xsmpe(type, nex) = 0.
          do nexout = 0, numex + 1
            mcontrib(type, nex, nexout) = 0.
            mpecontrib(type, nex, nexout) = 0.
          enddo
        enddo
        xsmpetot(type) = 0.
        if (parskip(type)) cycle
        Zix = Zindex(Zcomp, Ncomp, type)
        Nix = Nindex(Zcomp, Ncomp, type)
        if ( .not. strucexist(Zix, Nix)) then
          call structure(Zix, Nix)
          strucexist(Zix, Nix) = .true.
        endif
      enddo
      do nex = 0, numex
        xspopsave(nex) = 0.
        Dmulti(nex) = 0.
      enddo
      call exgrid(Zcomp, Ncomp)
!
! If a new set of transmission coefficients for each residual nuclide is requested, we calculate inverse cross sections and
! transmission coefficients for this compound nucleus here.
!
! basicxs   : subroutine for basic cross sections and transmission coefficients
!
      if (flagompall) call basicxs(Zcomp, Ncomp)
!
! Write population of compound nucleus before its decay
!
! levelsout    : subroutine for output of discrete levels
! densityout   : subroutine for output of level density parameters
! fissionparout: subroutine for output for fission parameters
!
      dExinc = deltaEx(Zcomp, Ncomp, maxex(Zcomp, Ncomp))
      odd = mod(A, 2)
      if (flagpop) then
        if ( .not. strucwrite(Zcomp, Ncomp)) then
          if (flaglevels) call levelsout(Zcomp, Ncomp)
          if (flagdensity .or. filedensity) call densityout(Zcomp, Ncomp)
          if (flagfisout) call fissionparout(Zcomp, Ncomp)
          strucwrite(Zcomp, Ncomp) = .true.
        endif
        write(*, '(/" Population of Z=", i3, " N=", i3, " (", i3, a2, ") before decay:", es12.5)') Z, N, A, nuc(Z), &
 &        xspopnuc(Zcomp, Ncomp)
        NL = Nlast(Zcomp, Ncomp, 0)
        if (maxex(Zcomp, Ncomp) > NL) then
          write(*, '(" Maximum excitation energy:", f8.3, " Discrete levels:", i3, " Continuum bins:", i3, &
 &          " Continuum bin size:", f8.3)') Exmax(Zcomp, Ncomp), NL, maxex(Zcomp, Ncomp) - NL, dExinc
        else
          write(*, '(" Maximum excitation energy:", f8.3, " Discrete levels:", i3)') Exmax(Zcomp, Ncomp), NL
        endif
        do parity = - 1, 1, 2
          write(*, '(/" Population of Z=", i3, " N=", i3, " (", i3, a2, &
 &          "), Parity=", i2, " before decay:", es12.5)') Z, N, A, nuc(Z), parity, xspopnucP(Zcomp, Ncomp, parity)
          write(*, '(/" bin    Ex     Pop     Pop P=", i2, 9("    J=", f4.1)/)') parity, (J+0.5*odd, J = 0, 8)
          do nex = 0, maxex(Zcomp, Ncomp)
            write(*, '(1x, i3, f8.3, 11es10.3)') nex, Ex(Zcomp, Ncomp, nex), xspopex(Zcomp, Ncomp, nex), &
 &            xspopexP(Zcomp, Ncomp, nex, parity), (xspop(Zcomp, Ncomp, nex, J, parity), J = 0, 8)
          enddo
        enddo
      endif
!
! Optional adjustment factors
!
      call isotrans(Z, N)
      do type = -1, 6
        if (adjustTJ(Zcomp, Ncomp, type)) then
          key = 'tjadjust'
          call adjust(Einc, key, Zcomp, Ncomp, type, 0, factor)
        else
          factor = 1.
        endif
        Fnorm(type) = factor / fisom(type)
      enddo
      if (flagpop) then
        write(*, '(/" Isospin factors to reduce emission for multiple emission for Z=", i3, " N=", i3, " (", i3, a2, ")", /)') &
 &        Z, N, A, nuc(Z)
        do type = 0, 6
          write(*,'(1x, a8, 1x, f8.5)') parname(type), fisom(type)
        enddo
      endif
      do type = 0, 6
        fisom(type) = fisominit(type)
      enddo
!
! Loop: De-excitation of nucleus, starting at the highest excitation energy bin.
!
! Continue for this (Zcomp,Ncomp,nex) only when there is sufficient reaction flux in the excitation energy bin.
! The fission transmission coefficients and level densities only need to be calculated once, at the highest excitation energy.
!
      popepsA = popeps / max(5 * maxex(Zcomp, Ncomp), 1)
      idensfis = 1
      if (flagpop) then
        write(*, '(/" Population of each bin of Z=", i3, " N=", i3, " (", i3, a2, ") before its decay"/)') Z, N, A, nuc(Z)
        write(*, '(" bin    Ex      Pop     P Pop per P", 9("  J=", f4.1, "  ")/)') (J+0.5*odd, J = 0, 8)
      endif
      Smin = S(Zcomp, Ncomp, 1)
      do nex = maxex(Zcomp, Ncomp), 1, - 1
        dExinc = deltaEx(Zcomp, Ncomp, nex)
        if (flagpop) then
          xspopsave(nex) = xspopex(Zcomp, Ncomp, nex)
          do parity = - 1, 1, 2
            write(*, '(1x, i3, f8.3, es10.3, i3, 10es10.3)') nex, Ex(Zcomp, Ncomp, nex), xspopex(Zcomp, Ncomp, nex), &
 &            parity, xspopexP(Zcomp, Ncomp, nex, parity), (xspop(Zcomp, Ncomp, nex, J, parity), J = 0, 8)
          enddo
        endif
!
! For exclusive channel cross section calculations, some variables need to be stored in extra arrays.
!
        if (flagchannels) popexcl(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex)
!
! Discrete levels decay by gamma cascade. Isomers are excluded from gamma cascade.
! Note that we assume that discrete levels cannot particle decay.
! This needs to be added (for light nuclei) in a future version.
!
! cascade: subroutine for gamma-ray cascade
!
        Exinc = Ex(Zcomp, Ncomp, nex)
        if (nex <= Nlast(Zcomp, Ncomp, 0) .and. Exinc <= Smin) then
          if (tau(Zcomp, Ncomp, nex) == 0.) call cascade(Zcomp, Ncomp, nex)
          cycle
        endif
        if (xspopex(Zcomp, Ncomp, nex) < popepsA) cycle
!
! For each mother excitation energy bin, determine the highest possible excitation energy bin nexmax for the residual nuclei.
! As reference, we take the top of the mother bin.
! The maximal residual excitation energy Exm is obtained by subtracting the separation energy from this.
! The bin in which this Exm falls is denoted by nexmax.
!
Loop1:  do type = 1, 6
          if (parskip(type)) cycle
          Exm = Exinc + 0.5 * dExinc - S(Zcomp, Ncomp, type)
          if (type > 1) Exm = Exm - egrid(ebegin(type))
          Zix = Zindex(Zcomp, Ncomp, type)
          Nix = Nindex(Zcomp, Ncomp, type)
          do nexout = maxex(Zix, Nix), 0, - 1
            dEx = deltaEx(Zix, Nix, nexout)
            Exmin = Ex(Zix, Nix, nexout) - 0.5 * dEx
            if (Exmin < Exm) then
              nexmax(type) = nexout
              cycle Loop1
            endif
          enddo
          nexmax(type) = - 1
        enddo Loop1
        nexmax(0) = nex - 1
!
! Prepare transmission coefficients and level densities.
! This can be done outside the loops over spin and parity.
!
! densprepare: subroutine to prepare energy grid and level density information for compound nucleus
!
        call densprepare(Zcomp, Ncomp, idensfis)
        idensfis = 0
!
! Deplete excitation energy bin by multiple pre-equilibrium.
!
! multipreeq2: subroutine for two-component multiple preequilibrium
! multipreeq : subroutine for multiple preequilibrium
!
        if (mulpreZN(Zcomp, Ncomp)) then
          if (flag2comp) then
            call multipreeq2(Zcomp, Ncomp, nex)
          else
            call multipreeq(Zcomp, Ncomp, nex)
          endif
        endif
!
! Compound nucleus decay of mother excitation energy/spin/parity bin.
!
        if (flagcomp) then
          popepsB = popepsA / (5 * maxJ(Zcomp, Ncomp, nex)) * 0.5
          do parity = - 1, 1, 2
            do J = 0, maxJ(Zcomp, Ncomp, nex)
!
! Continue for this (Zcomp,Ncomp,nex,J,P) only when there is sufficient reaction flux in the excitation energy bin.
! The correct value for J is determined.
!
! tfission   : subroutine for fission transmission coefficients
! compound   : subroutine for Hauser-Feshbach model for multiple emission
! tfissionout: subroutine for output of fission transmission coefficients
!
              popdecay = 0.
              partdecay = 0.
              partdecaytot = 0.
              if (xspop(Zcomp, Ncomp, nex, J, parity) < popepsB) cycle
              xsp = xspop(Zcomp, Ncomp, nex, J, parity)
              J2 = 2 * J + odd
              if (flagfission .and. nfisbar(Zcomp, Ncomp) /= 0) call tfission(Zcomp, Ncomp, nex, J2, parity)
              call compound(Zcomp, Ncomp, nex, J2, parity)
              if (flagpop) then
                rJ=0.5*J2
                write(*,'(f4.1,i4,8es10.3)') rJ,parity,xsp,(partdecaytot(type),type=0,6)
                if (flagdecay) then
                  do type = 0, 6
                    if (parskip(type)) cycle
                    Zix = Zindex(Zcomp, Ncomp, type)
                    Nix = Nindex(Zcomp, Ncomp, type)
                    Zres = ZZ(Zcomp, Ncomp, type)
                    Nres = NN(Zcomp, Ncomp, type)
                    Ares = AA(Zcomp, Ncomp, type)
                    oddres = mod(Ares, 2)
                    do Pres = - 1, 1, 2
                      write(*,'(/" Total Pprime=",i2,":",es10.3," via ",a8," emission"/)') Pres,partdecay(type,Pres),parname(type)
                      write(*,'(" bin    Ex",10("    J=",f4.1)/)') (Ir+0.5*oddres,Ir=0,9)
                      do nexout = 0, nexmax(type)
                        write(*, '(1x, i3, f8.3, 10es10.3)') nexout, Ex(Zix, Nix, nexout), &
 &                        (popdecay(type, nexout, Jres, Pres), Jres = 0, 9)
                      enddo
                      write( * , * )
                    enddo
                  enddo
                endif
              endif
            enddo
          enddo
          if (flagfisout) call tfissionout(Zcomp, Ncomp, nex)
        endif
      enddo
!
! Make new population cross section per nucleus
!
      xspopnuc(Zcomp, Ncomp) = xspopex(Zcomp, Ncomp, 0)
      do nex = 1, Nlast(Zcomp, Ncomp, 0)
        if (tau(Zcomp, Ncomp, nex) /= 0.) xspopnuc(Zcomp, Ncomp) = xspopnuc(Zcomp, Ncomp) + xspopex(Zcomp, Ncomp, nex)
      enddo
      if (flagcompo) then
        xsmax = 0.
        Smax = 0.
        do type = 1, 6
          if (xsfeed(Zcomp, Ncomp, type) > xsmax) then
            xsmax = xsfeed(Zcomp, Ncomp, type)
            Smax = S(Zcomp, Ncomp, type)
          endif
        enddo
        do nex = 0, maxex(Zcomp, Ncomp)
          if (Ex(Zcomp, Ncomp, nex) > Smax) exit
          xspoppreeq(Zcomp, Ncomp) = xspoppreeq(Zcomp, Ncomp) + preeqpopex(Zcomp, Ncomp, nex)
        enddo
        xspoppreeq(Zcomp, Ncomp) = min(dble(xspoppreeq(Zcomp, Ncomp)), xspopnuc(Zcomp, Ncomp) - dble(xspopdir(Zcomp, Ncomp)))
        xspopcomp(Zcomp, Ncomp) = max(xspopnuc(Zcomp, Ncomp) - dble(xspoppreeq(Zcomp, Ncomp) - xspopdir(Zcomp, Ncomp)), 0.d0)
        if (xspopnuc(Zcomp, Ncomp) > 0.) then
          Fdir(Zcomp, Ncomp) = xspopdir(Zcomp, Ncomp) / xspopnuc(Zcomp, Ncomp)
          Fpreeq(Zcomp, Ncomp) = xspoppreeq(Zcomp, Ncomp) / xspopnuc(Zcomp, Ncomp)
          Fcomp(Zcomp, Ncomp) = xspopcomp(Zcomp, Ncomp) / xspopnuc(Zcomp, Ncomp)
        endif
      endif
!
! For exclusive channel calculation: Determine feeding terms.
!
      if (flagchannels .and. Zcomp <= numZchan .and. Ncomp <= numNchan) then
        do nex = maxex(Zcomp, Ncomp), 1, - 1
          do type = 0, 6
            if (parskip(type)) cycle
            if (Zcomp == 0 .and. Ncomp == 0 .and. type > 0 .and. .not. flaginitpop) cycle
            Zix = Zindex(Zcomp, Ncomp, type)
            Nix = Nindex(Zcomp, Ncomp, type)
            do nexout = 0, maxex(Zix, Nix)
              feedexcl(Zcomp, Ncomp, type, nex, nexout) = mcontrib(type, nex, nexout)
            enddo
          enddo
        enddo
      endif
      if (Zcomp == 0 .and. Ncomp == 0 .and. flagfission) fisfeedex(0, 0, maxex(0, 0) + 1) = xsbinary( - 1)
!
! Increase emission spectra after decay of mother excitation energy bin.
!
! compemission: subroutine compound emission spectra for continuum
! kalbach     : Kalbach function
!
      if (flagspec) then
        call compemission(Zcomp, Ncomp)
        do type = 0, 6
          if (parskip(type)) cycle
          emissum(type) = 0.
          Eaveragemul(Zcomp, Ncomp, type) = 0.
          Eaveragesum = 0.
          do nen = ebegin(type), eend(type)
            xscomp(type, nen) = xscomp(type, nen) + xsemis(type, nen)
            xsmpreeq(type, nen) = xsmpreeq(type, nen) + xsmpeemis(type, nen)
            if (flagddx .or. flagrecoil) then
              Eout = egrid(nen)
              do iang = 0, nanglecont
                ang = anglecont(iang) * deg2rad
                xscompad(type, nen, iang) = xscompad(type, nen, iang) + xsemis(type, nen) / fourpi
                xsmpreeqad(type, nen, iang) = xsmpreeqad(type, nen, iang) + &
                  xsmpeemis(type, nen) * kalbach(type, Einc, Eout, ang)
              enddo
            endif
            emissum(type) = emissum(type) + (xsemis(type, nen) + xsmpeemis(type, nen)) * deltaE(nen)
            Eaveragesum = Eaveragesum + egrid(nen) * (xsemis(type, nen) + xsmpeemis(type, nen)) * deltaE(nen)
          enddo
          if (emissum(type) > 0.) Eaveragemul(Zcomp, Ncomp, type) = Eaveragesum / emissum(type)
!
! For ENDF-6 files, exclusive (n,gn), (n,gp) ... (n,ga) spectra are required. This is stored here.
!
          if (flagendf .and. parinclude(0) .and. Zcomp == 0 .and. Ncomp == 0) then
              do nen = ebegin(0), eendhigh
                xsngnspec(type, nen) = xsemis(type, nen) + xsmpeemis(type, nen)
              enddo
          endif
        enddo
!
! Output of emission spectrum per bin
!
        if (flagbinspec) then
          do nex = maxex(Zcomp, Ncomp), 1, - 1
            write(*, '(/" Emission spectra from Z=", i3, " N=", i3, " (", i3, a2, "), Ex=", f12.5, " MeV"/)') &
 &            Z, N, A, nuc(Z), Ex(Zcomp, Ncomp, nex)
            write(*, '("  Energy ", 7(2x, a8, 2x)/)') (parname(type), type = 0, 6)
            do nen = ebegin(0), eendhigh
              write(*, '(1x, f8.3, 7es12.5)') egrid(nen), (xsbinspec(type, nex, nen), type = 0, 6)
            enddo
          enddo
        endif
      endif
!
! **** Write partial emission channels and production cross sections ***
!
! Particles
!
      Z = ZZ(Zcomp, Ncomp, 0)
      N = NN(Zcomp, Ncomp, 0)
      A = AA(Zcomp, Ncomp, 0)
      if (flagpop) then
        write(*, '(/" Emitted flux per excitation energy bin of", " Z=", i3, " N=", i3, " (", i3, a2, "):"/)') Z, N, A, nuc(Z)
        write(*, '(" bin    Ex    ", 7(2x, a8, 2x), "  Total     In - out"/)') (parname(type), type = 0, 6)
        do nex = 0, maxex(Zcomp, Ncomp)
          sumxs = 0.
          do type = 0, 6
            sumxs = sumxs + xspartial(type, nex)
          enddo
          xsdif = xspopsave(nex) - sumxs
          write(*, '(1x, i3, f8.3, 9es12.5)') nex, Ex(Zcomp, Ncomp, nex), (xspartial(type, nex), type = 0, 6), sumxs, xsdif
        enddo
!
! Fission
!
        if (flagfission) then
          write(*, '(/" Fission contribution from Z=", i3, " N=", i3, " (", i3, a2, "):"/)') Z, N, A, nuc(Z)
          write(*, '("   Ex    Popul. "/)')
          do nex = 0, maxex(Zcomp, Ncomp)
            write(*, '(1x, f8.3, es12.5)') Ex(Zcomp, Ncomp, nex), fisfeedex(Zcomp, Ncomp, nex)
          enddo
        endif
!
! Multiple pre-equilibrium emission
!
        if (mulpreZN(Zcomp, Ncomp)) then
          write(*, '(/" Multiple preequilibrium emission from ", "Z=", i3, " N=", i3, " (", i3, a2, "):"/)') Z, N, A, nuc(Z)
          write(*, '(61x, "Feeding terms from previous ", "particle-hole configuration"/)')
          if ( .not. flag2comp) then
            write(*, '(" bin    Ex  Mpe ratio  neutron   proton", "  ", 10("  ", i1, "p", i1, "h    "))') &
 &            ((p, h, p = 1, h), h = 1, 4)
            write(*, '("                      emission  emission"/)')
            do nex = Nlast(Zcomp, Ncomp, 0) + 1, maxex(Zcomp, Ncomp)
              write(*, '(1x, i3, f8.3, f8.5, 12es10.3)') nex, Ex(Zcomp, Ncomp, nex), Dmulti(nex), xsmpe(1, nex), &
 &              xsmpe(2, nex), ((xspopph(Zcomp, Ncomp, nex, p, h), p = 1, h), h = 1, 4)
              Dmulti(nex) = 0.
            enddo
          else
            write(*, '(" bin    Ex  Mpe ratio  neutron   proton", 11(2x, 4i2))') 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, &
 &            0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 2, 1, 0, 0, 0, 0, 2, 1, 2, 2, 0, 0, 0, 0, 2, 2, 1, 1, 1, 1
            write(*, '("                      emission  emission"/)')
            do nex = Nlast(Zcomp, Ncomp, 0) + 1, maxex(Zcomp, Ncomp)
              write(*, '(1x, i3, f8.3, f8.5, 13es10.3)') nex, Ex(Zcomp, Ncomp, nex), Dmulti(nex), xsmpe(1, nex), &
 &              xsmpe(2, nex), xspopph2(Zcomp, Ncomp, nex, 1, 1, 0, 0), &
 &              xspopph2(Zcomp, Ncomp, nex, 0, 0, 1, 1), xspopph2(Zcomp, Ncomp, nex, 1, 0, 0, 1), &
                xspopph2(Zcomp, Ncomp, nex, 0, 1, 1, 0), xspopph2(Zcomp, Ncomp, nex, 1, 1, 1, 0), &
                xspopph2(Zcomp, Ncomp, nex, 1, 0, 1, 1), xspopph2(Zcomp, Ncomp, nex, 2, 1, 0, 0), &
                xspopph2(Zcomp, Ncomp, nex, 0, 0, 2, 1), xspopph2(Zcomp, Ncomp, nex, 2, 2, 0, 0), &
                xspopph2(Zcomp, Ncomp, nex, 0, 0, 2, 2), xspopph2(Zcomp, Ncomp, nex, 1, 1, 1, 1)
              Dmulti(nex) = 0.
            enddo
          endif
          write(*, '(/"  Total             ", 2es10.3)') xsmpetot(1), xsmpetot(2)
        endif
!
! Total decay from mother nucleus
!
        write(*, '(/" Emission cross sections to residual ", "nuclei from Z=", i3, " N=", i3, " (", i3, a2, "):"/)') Z, N, A, nuc(Z)
        if (flagfission) write(*, '(" fission  channel", 23x, ":", es12.5)') xsfeed(Zcomp, Ncomp, - 1)
        do type = 0, 6
          if (parskip(type)) cycle
          Z = ZZ(Zcomp, Ncomp, type)
          N = NN(Zcomp, Ncomp, type)
          A = AA(Zcomp, Ncomp, type)
          write(*, '(1x, a8, " channel to Z=", i3, " N=", i3, " (", i3, a2, "):", es12.5)') &
 &          parname(type), Z, N, A, nuc(Z), xsfeed(Zcomp, Ncomp, type)
        enddo
!
! Total emission spectra from mother nucleus.
!
        Z = ZZ(Zcomp, Ncomp, 0)
        N = NN(Zcomp, Ncomp, 0)
        A = AA(Zcomp, Ncomp, 0)
        if (flagspec) then
          write(*, '(/" Emission spectra from Z=", i3, " N=", i3, " (", i3, a2, "):"/)') Z, N, A, nuc(Z)
          write(*, '("  Energy ", 7(2x, a8, 2x)/)') (parname(type), type = 0, 6)
          do nen = ebegin(0), eendhigh
            write(*, '(1x, f8.3, 7es12.5)') egrid(nen), (xsemis(type, nen) + xsmpeemis(type, nen), type = 0, 6)
          enddo
          if (flagcheck) then
            write(*, '(/" ++++++++++ CHECK OF INTEGRATED ", "EMISSION SPECTRA ++++++++++"/)')
            write(*, '(13x, "Cross section   Integrated spectrum", "  Average emission energy"/)')
            do type = 0, 6
              if (parskip(type)) cycle
              write(*, '(1x, a8, 2(4x, es12.5), 10x, f8.3)') parname(type), xsfeed(Zcomp, Ncomp, type), emissum(type), &
 &              Eaveragemul(Zcomp, Ncomp, type)
            enddo
          endif
        endif
!
! Final production cross section for ground state and isomer of nucleus after decay.
!
        write(*, '(/" Final production cross section of Z=", i3, " N=", i3, " (", i3, a2, "):"/)') Z, N, A, nuc(Z)
        write(*, '(" Total       :", es12.5)') xspopnuc(Zcomp, Ncomp)
        write(*, '(" Ground state:", es12.5)') xspopex(Zcomp, Ncomp, 0)
        do nex = 1, Nlast(Zcomp, Ncomp, 0)
          if (tau(Zcomp, Ncomp, nex) /= 0.) write(*, '(" Level", i3, "    :", es12.5)') levnum(Zcomp, Ncomp, nex), &
 &          xspopex(Zcomp, Ncomp, nex)
        enddo
      endif
!
! Fission
!
      if (flagfission .and. flagfisfeed) then
        fisfile = 'fis000000.nex'
        write(fisfile(4:6), '(i3.3)') Z
        write(fisfile(7:9), '(i3.3)') A
        open (unit = 1, file = fisfile, status = 'replace')
        write(1, '("# Reaction: ", g12.4, " MeV ", a8, " on Z=", i3, " N=", i3, " (", i3, a2, ")")') &
 &        einc, parname(k0), Ztarget, Ntarget, Atarget, nuc(Ztarget)
        write(1, '("# Fission contribution from Z=", i3, " N=", i3, " (", i3, a2, ")")') Z, N, A, nuc(Z)
        if (Zcomp == 0 .and. Ncomp == 0) then
          nen = maxex(Zcomp, Ncomp) + 1
        else
          nen = maxex(Zcomp, Ncomp)
        endif
        Jfis = 0
        do nex = 0, nen
          do parity = - 1, 1, 2
            do J = 0, numJ
              if (fisfeedJP(Zcomp, Ncomp, nex, J, parity) > 0.) Jfis = max(Jfis, J)
            enddo
          enddo
        enddo
        write(1, '("# # energies =", i6)') nen+1
        write(1, '("# # spins    =", i4)') Jfis+1
        form1='("#    Ex   Population",xx(5x,i2,"+",9x,i2,"-",4x))'
        write(form1(25:26), '(i2.2)') Jfis+1
        form2='(1x,f8.3,es12.5,xx(2es12.5))'
        write(form2(17:18), '(i2.2)') Jfis+1
        write(1, fmt = form1) (J, J, J = 0, Jfis)
        do nex = 0, nen
          write(1, fmt = form2) Ex(Zcomp, Ncomp, nex), fisfeedex(Zcomp, Ncomp, nex), &
 &          ((fisfeedJP(Zcomp, Ncomp, nex, J, parity), parity = - 1, 1, 2), J = 0, Jfis)
        enddo
        close (unit = 1)
      endif
!
! ******* Add binary cross sections to initial compound nucleus ********
!
! This is done to keep proper track of exclusive cross sections and channels such as (n,gn).
!
  500 if (Zcomp == 0 .and. Ncomp == 0 .and. .not. flaginitpop) then
        xsngnsum = 0.
        do type = - 1, 6
          if (parinclude(type)) then
            xsngn(type) = xsfeed(0, 0, type)
            if (type /= 0) xsngnsum = xsngnsum + xsngn(type)
            xsfeed(0, 0, type) = xsfeed(0, 0, type) + xsbinary(type)
          endif
        enddo
        if (flagchannels) then
          popexcl(0, 0, maxex(0, 0) + 1) = xsreacinc
          if (flagfission) fisfeedex(0, 0, maxex(0, 0) + 1) = xsbinary( - 1)
          do type = 0, 6
            if (parskip(type)) cycle
            Zix = Zindex(Zcomp, Ncomp, type)
            Nix = Nindex(Zcomp, Ncomp, type)
            do nexout = 0, maxex(Zix, Nix)
              feedexcl(0, 0, type, maxex(0, 0) + 1, nexout) = feedbinary(type, nexout)
            enddo
          enddo
        endif
      endif
    enddo
  enddo
  return
end subroutine multiple
! Copyright A.J. Koning 2021
