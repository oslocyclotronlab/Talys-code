      subroutine gammaout(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 15, 2016
c | Task  : Output of gamma-ray strength functions, transmission
c |         coefficients and cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*25 model
      integer      Zcomp,Ncomp,Z,N,A,l,nen
      real         e,fstrength
c
c ***** Gamma-ray strength functions and transmission coefficients *****
c
c Zcomp     : charge number index for compound nucleus
c Ncomp     : neutron number index for compound nucleus
c ZZ,Z      : charge number of residual nucleus
c NN,N      : neutron number of residual nucleus
c AA,A      : mass number of residual nucleus
c nuc       : symbol of nucleus
c gamgam    : total radiative width
c dgamgam   : uncertainty in gamgam
c gamgamth  : theoretical total radiative width
c D0        : experimental s-wave resonance spacing in eV
c dD0       : uncertainty in D0
c D0theo    : theoretical s-wave resonance spacing
c swaveth   : theoretical strength function for s-wave
c gnorm     : gamma normalization factor
c strength  : model for E1 gamma-ray strength function
c model     : string for gamma-ray strength function
c strengthM1: model for M1 gamma-ray strength function
c etable    : constant to adjust tabulated strength functions
c ftable    : constant to adjust tabulated strength functions
c upbend    : properties of the low-energy upbend of given multipolarity
c nTqrpa    : number of temperatures for QRPA
c gammax    : number of l-values for gamma multipolarity
c sgr       : strength of GR
c ngr       : number of GR
c egr       : energy of GR
c ggr       : width of GR
c epr       : energy of PR
c gpr       : width of PR
c tpr       : strength of PR
c kgr       : constant for gamma-ray strength function
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c egrid,e   : outgoing energy grid
c fstrength : gamma-ray strength function
c Tjl       : transmission coefficients as a function of particle type,
c             energy, spin and l-value
c
      Z=ZZ(Zcomp,Ncomp,0)
      N=NN(Zcomp,Ncomp,0)
      A=AA(Zcomp,Ncomp,0)
      write(*,'(/" ########## GAMMA STRENGTH FUNCTIONS, TRANSMISSION",
     +  " COEFFICIENTS AND CROSS SECTIONS ##########")')
      write(*,'(/" Gamma-ray information for Z=",i3," N=",i3,
     +  " (",i3,a2,") "/)') Z,N,A,nuc(Z)
      write(*,'(" S-wave strength function parameters:"/)')
      write(*,'(" Exp. total radiative width=",f10.5," eV +/-",f8.5,
     +  " Theor. total radiative width for l=0:",f15.5," eV")')
     +  gamgam(Zcomp,Ncomp),dgamgam(Zcomp,Ncomp),gamgamth(Zcomp,Ncomp,0)
      write(*,'(53x," Theor. total radiative width for l=1:",f15.5,
     +  " eV")') gamgamth(Zcomp,Ncomp,1)
      write(*,'(" Exp. D0                   =",f10.2," eV +/-",f8.2,
     +  " Theor. D0                   =",f15.5," eV")')
     +  D0(Zcomp,Ncomp),dD0(Zcomp,Ncomp),D0theo(Zcomp,Ncomp)
      write(*,'(53x," Theor. D1                   =",f15.5," eV")')
     +  D1theo(Zcomp,Ncomp)
      write(*,'(" Theor. S-wave strength f. =",f10.5,"E-4")')
     +  1.e4*swaveth(Zcomp,Ncomp)
      write(*,'(" Normalization factor      =",f10.5)') gnorm
      if (strength.eq.1) model="Kopecky-Uhl              "
      if (strength.eq.2) model="Brink-Axel               "
      if (strength.eq.3) model="Goriely HFbcs tables     "
      if (strength.eq.4) model="Goriely HFB tables       "
      if (strength.eq.5) model="Goriely Hybrid model     "
      if (strength.eq.6) model="Goriely T-dep. HFB Tables"
      if (strength.eq.7) model="Goriely T-dep. RMF Tables"
      if (strength.eq.8) model="Gogny D1M HFB+QRPA Tables"
      write(*,'(/" Gamma-ray strength function model for E1: ",a25)')
     +  model
      if (strengthM1.eq.1) model="RIPL-1                   "
      if (strengthM1.eq.2) model="RIPL-2                   "
      write(*,'(/" Gamma-ray strength function model for M1: ",a25)')
     +  model
      if (strengthM1.eq.8) model="Gogny D1M HFB+QRPA Tables"
      write(*,'(/" Gamma-ray strength function model for M1: ",a25)')
     +  model
      if (Z.eq.N)
     +  write(*,'(/" Isospin effect on g-strength for incident ",a1,
     +  ": fiso=",f5.2)') parsym(k0),fiso(Zcomp,Ncomp,k0)
      if (strength.eq.3.or.strength.eq.4.or.strength.ge.6) then
        write(*,'(/" Adjustable parameters for E1 : etable=",f10.5,
     +    " ftable=",f10.5," number of T=",i3)')
     +    etable(Zcomp,Ncomp,1,1),ftable(Zcomp,Ncomp,1,1),nTqrpa
      endif
      if (strengthM1.eq.8) then
        write(*,'(/" Adjustable parameters for M1 : etable=",f10.5,
     +    " ftable=",f10.5)')
     +    etable(Zcomp,Ncomp,0,1),ftable(Zcomp,Ncomp,0,1)
      endif
      if (upbend(Zcomp,Ncomp,1,1,1).ne.0.) then
        write(*,'(/" Inclusion of an E1 upbend C exp(-eta*E) with C=",
     +    es10.2," eta=",f8.2)')
     +    upbend(Zcomp,Ncomp,1,1,1),upbend(Zcomp,Ncomp,1,1,2)
      endif
      if (upbend(Zcomp,Ncomp,0,1,1).ne.0.) then
        write(*,'(/" Inclusion of an M1 upbend C exp(-eta*E) with C=",
     +    es10.2," eta=",f8.2)')
     +    upbend(Zcomp,Ncomp,0,1,1),upbend(Zcomp,Ncomp,0,1,2)
      endif
      do 10 l=1,gammax
        write(*,'(/" Normalized gamma-ray strength functions and ",
     +    "transmission coefficients for l=",i2,/)') l
        write(*,'(" Giant resonance parameters :"/)')
        if (ngr(Zcomp,Ncomp,1,l).eq.2) then
          write(*,'(" sigma0(M",i1,") =",f8.3,"       sigma0(E",i1,
     +      ") =",f8.3," and ",f8.3,
     +      "    PR: sigma0(M",i1,") =",f8.3,"       sigma0(E",i1,
     +      ") =",f8.3)')
     +      l,sgr(Zcomp,Ncomp,0,l,1),l,
     +      sgr(Zcomp,Ncomp,1,l,1),sgr(Zcomp,Ncomp,1,l,2),
     +      l,tpr(Zcomp,Ncomp,0,l,1),l,tpr(Zcomp,Ncomp,1,l,1)
        else
          write(*,'(" sigma0(M",i1,") =",f8.3,"       sigma0(E",i1,
     +      ") =",f8.3,
     +      "    PR: sigma0(M",i1,") =",f8.3,"       sigma0(E",i1,
     +      ") =",f8.3)')
     +      l,sgr(Zcomp,Ncomp,0,l,1),l,sgr(Zcomp,Ncomp,1,l,1),
     +      l,tpr(Zcomp,Ncomp,0,l,1),l,tpr(Zcomp,Ncomp,1,l,1)
        endif
        if (ngr(Zcomp,Ncomp,1,l).eq.2) then
          write(*,'("      E(M",i1,") =",f8.3,"            E(E",i1,
     +      ") =",f8.3," and ",f8.3,
     +    "    PR:      E(M",i1,") =",f8.3,"            E(E",i1,
     +      ") =",f8.3)')
     +      l,egr(Zcomp,Ncomp,0,l,1),l,
     +      egr(Zcomp,Ncomp,1,l,1),egr(Zcomp,Ncomp,1,l,2),
     +      l,epr(Zcomp,Ncomp,0,l,1),l,epr(Zcomp,Ncomp,1,l,1)
        else
          write(*,'("      E(M",i1,") =",f8.3,"            E(E",i1,
     +      ") =",f8.3,
     +    "    PR:      E(M",i1,") =",f8.3,"            E(E",i1,
     +      ") =",f8.3)')
     +      l,egr(Zcomp,Ncomp,0,l,1),l,egr(Zcomp,Ncomp,1,l,1),
     +      l,epr(Zcomp,Ncomp,0,l,1),l,epr(Zcomp,Ncomp,1,l,1)
        endif
        if (ngr(Zcomp,Ncomp,1,l).eq.2) then
          write(*,'("  gamma(M",i1,") =",f8.3,"        gamma(E",i1,
     +      ") =",f8.3," and ",f8.3,
     +      "    PR:  gamma(M",i1,") =",f8.3,"        gamma(E",i1,") =",
     +      f8.3)')
     +      l,ggr(Zcomp,Ncomp,0,l,1),l,
     +      ggr(Zcomp,Ncomp,1,l,1),ggr(Zcomp,Ncomp,1,l,2),
     +      l,gpr(Zcomp,Ncomp,0,l,1),l,gpr(Zcomp,Ncomp,1,l,1)
        else
          write(*,'("  gamma(M",i1,") =",f8.3,"        gamma(E",i1,
     +      ") =",f8.3,
     +      "    PR:  gamma(M",i1,") =",f8.3,"        gamma(E",i1,") =",
     +      f8.3)')
     +      l,ggr(Zcomp,Ncomp,0,l,1),l,ggr(Zcomp,Ncomp,1,l,1),
     +      l,gpr(Zcomp,Ncomp,0,l,1),l,gpr(Zcomp,Ncomp,1,l,1)
        endif
        write(*,'("      k(M",i1,") =",es14.5,"      k(E",i1,") =",
     +    es14.5/)') l,kgr(Zcomp,Ncomp,0,l),l,kgr(Zcomp,Ncomp,1,l)
        write(*,'("      E       f(M",i1,")        f(E",i1,")",
     +    "        T(M",i1,")        T(E",i1,")"/)')  l,l,l,l
        do 20 nen=ebegin(0),eend(0)
          e=egrid(nen)
          write(*,'(1x,f8.3,4es13.5)') e,
     +      fstrength(Zcomp,Ncomp,0.,e,0,l)*gnorm,
     +      fstrength(Zcomp,Ncomp,0.,e,1,l)*gnorm,Tjl(0,nen,0,l),
     +      Tjl(0,nen,1,l)
   20   continue
   10 continue
c
c **************** Cross sections for inverse channels *****************
c
c xsreac: reaction cross section
c
      write(*,'(/" Photoabsorption cross sections"/)')
      write(*,'("     E      reaction"/)')
      do 110 nen=ebegin(0),eend(0)
        write(*,'(1x,f8.3,es12.4)') egrid(nen),xsreac(0,nen)
  110 continue
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely
