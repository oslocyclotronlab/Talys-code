      subroutine breakupAVR
c
c +---------------------------------------------------------------------
c | Author: Marilena Avrigeanu
c | Date  : October 10 2015
c | Task  : Deuteron breakup fractions, Avrigeanu model
c +---------------------------------------------------------------------
c
c This model is based on  Eqs. (1-3), PRC 88,014612(2013)
c            and cross sections, Eqs. (2-5), PRC 89,044613 (2014)
c This module includes also subroutine checkBU and function GAUSS
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,nen
      real    fracBUp,fracBUE,fracBUBF,fracBUT,EnormEB,RnormEB,TnormBU
      real    ebreakCM(2),ebreakLS(2)
      real    spec03LS(2,0:numen),spec03CM(2,0:numen),
     +        speccorCM(2),speccorLS(2),sumtest(2),BCout(2)
      real    Sigr,EmaxnLS,EmaxnCM,E0n03,xsaBU,xsBUT,
     +        Gauss,E0n03CM,width,w03,Eout,En,EnCM,
     +        BCin,Bdeut,AT13,arg,arg2,argcm
c
c ************************** Avrigeanu model ***************************
c
      write(8,*)"    "
      write(8,*)'  DEUTERON break-up Parameterization, Avrigeanu',
     1   'et al., PRC88,014612, PRC89,044613, and Refs. therein'
c
c DEUTERON Break-up model, Avrigeanu et al., PRC88,014612 & PRC89,044613
c
c Calculation of terms independent of emission energy.
c
c Einc            : incident energy in MeV
c Atarget         : mass number of target nucleus
c BCin,BCout      : effective incident/outgoing Coulomb barrier
c Ztarget         : charge number of target nucleus
c xsreacinc, Sigr : reaction cross section for incident channel
c Bdeut           : deuteon binding energy
c AT13            : A**1/3
c EmaxnCM         : maximum breakup nucleons energy, C.M. System
c EmaxnLS         : maximum breakup nucleons energy, Lab. System
c fracBUE         : elastic breakup fraction
c fracBUBF        : nucleon inelastic breakup fraction
c fracBUp         : nucleon breakup fraction
c xsEB            : elastic breakup cross section
c xsBF            : nucleon inelastic breakup cross section
c xsBUnuc         : nucleon breakup cross section
c xsBUT           : TOTAL breakup cross section
c Gauss           : function for Gaussian
c xsaBU           : breakup cross section
c
c       xsBUnuc = xsEB + xsBF,  Eq. 4  PRC89, 044613.
c
c Avoiding to count xsEB twice, TOTAL breakup cross section is:
c       xsBUT = xsEB + 2*xsBF,  Eq. 5  PRC89, 044613.
c
      Sigr=xsreacinc
      Bdeut=2.225
      BCin=Ztarget/9.5
      BCout(1)=0.
      BCout(2)=BCin
      AT13=Atarget**onethird
c
c maximum breakup nucleon energy:
c Laboratory System
c
      EmaxnLS=Einc*(Atarget+1.)/(Atarget+2.)-Bdeut*(Atarget+1.)/Atarget
c
c Center of Mass System
c
      EmaxnCM  = Einc*Atarget/(Atarget+2.)-Bdeut
      if(EmaxnLS.le.0.) then
         EmaxnLS=0.
         write(8,*)" deuteron energy lower than Bd, no breakup"
         write(*,*)" deuteron energy lower than Bd, no breakup"
         go to 399
      endif
c
      write(8,*)"    "
      write(8,*)" deuteron energy & maximum outgoing fragments energy"
      write(8,*)"    Einc   EmaxnLS   EmaxnCM  "
      write(8,889)Einc,EmaxnLS,EmaxnCM
c
      RnormEB=1.
      TnormBU=1.
      xsBUT=0.
      xsaBU=0.
c
      do 103 type =1,2
        ebreakLS(type)=0.
        ebreakCM(type)=0.
        sumtest(type)=0.
        speccorCM(type)=0.
        speccorLS(type)=0.
        do 103 nen=0,eend(type)
          spec03LS(type,nen)=0.
          spec03CM(type,nen)=0.
          xspreeqbu(type,nen)=0.
103   continue
c
      call checkBU(EnormEB,RnormEB,TnormBU)
c
      write(8,*)" Ei_EB_norm   EB_norm ","  total_BU_norm "
      write(8,889)EnormEB,RnormEB,TnormBU
c
c arg2: squre of energy
c argcm: C.M. energy
c Tnorm: normalization factor
c
      arg=Einc
      arg2=arg*arg
      argcm=arg*Atarget/(Atarget+2)
c
c      ********   breakup fractions        ********
c
c       NUCLEON TOTAL breakup fraction, PRC88,014612, Eq. 1
       fracBUp=0.087-0.0066*Ztarget+0.00163*Ztarget*AT13
     +        +0.0017*AT13*arg-0.000002*Ztarget*arg2
      fracBUp = TnormBU*fracBUp
c
      if (fracBUp.le.0.) then
         fracBUp = 0.
         go to 499
      endif
c
c       ELASTIC breakup fraction, PRC88,014612, Eq. 2
        fracBUE=0.031-0.0028*Ztarget+0.00051*Ztarget*AT13
     1        +0.0005*AT13*arg-0.000001*Ztarget*arg2
c
      if(Einc.lt.EnormEB) then
         fracBUE=TnormBU*fracBUE
      else
         fracBUE=TnormBU*RnormEB*fracBUp
      endif
c
c       nucleon INELASTIC breakup fraction, PRC88,014612, Eq. 3
      fracBUBF = fracBUp-fracBUE
c
c       TOTAL breakup fraction,  PRC89,044613, Eq. 5
c fracBUT: total breakup fraction
c
      fracBUT = 2.D+00*fracBUBF+fracBUE
c
c      ********     end    breakup fractions       ********
c
c   ******************************
c    TOTAL BREAKUP cross section
      xsBUT=Sigr*fracBUT
c    crossafterbreakup  xsaBU
      xsaBU=Sigr-xsBUT
c   ******************************
c
c
c-------------           deuteron fragments   loop         --------
c
      do 301, type=1,2
c
c   ******************************
c
c    ELASTIC breakup cross section
cccc      sigEB(type)=Sigr*fracBUE
        xsEB(type)=Sigr*fracBUE
c    INELASTIC breakup cross section (BF) equal for neutron and proton
cccc      sigBF(type)=Sigr*fracBUBF
        xsBF(type)=Sigr*fracBUBF
c    TOTAL NUCLEON breakup cross section = xsEB+xsBF
cccc      sigBUp(type)=Sigr*fracBUp
        xsBUnuc(type)=Sigr*fracBUp
c
c   ******************************
c
      if(type.eq.1) then
         write(8,*)" breakup fractions, equal for neutron and",
     1      "  protons, and the corresponding cross sections"
            write(8,*)"   Einc    fracBUp   fracBUE",
     1      "   fracBUBF xsBUnuc       xsEB      xsBF"
            write(8,889)Einc,fracBUp,fracBUE,fracBUBF,
     1            xsBUnuc(1),xsEB(1),xsBF(1)
          write(8,*)" Total nucleon(n/p) emission break-up",
     1      " xsBUnucc, sum of inelastic (BF) & elastic (EB)",
     2      " break-up components:"
          write(8,*)"    xsBUnuc = xsEB + xsBF = ",xsBUnuc(Type),
     1       ",   PRC89, 044613, Eq. 4. "
          write(8,*)" Avoiding to count xsEB twice, TOTAL  ",
     1      " breakup cross section is:"
          write(8,*)"    xsBUT = xsEB + 2*xsBF = ", xsBUT,
     1      ",   PRC89, 044613, Eq. 5."
      endif
c
c  ebreakLS, ebreakCM  : Centroid energy of the breakup nucleon energy
c                        distributions, Kalbach 2003 in Laboratory, and
c                        respectively Center of Mass Systems
c   width              : Full width at half maximum of the breakup nucleon
c                        energy distribution, Kalbach 2003
c
      ebreakCM(type) = (0.5*(argcm - Bdeut - BCin) + BCout(type))
      ebreakLS(type) = 0.5*(Atarget+1.)/(Atarget+2.)*Einc+
     1   0.5*(Atarget+1.)/Atarget*(-Bdeut-BCin+2.*BCout(type))
      width = 1.15 + 0.12*Einc - Atarget/140.
c
      if (width.le.0.) then
           width=0.
           go to 299
      endif
c
      if (ebreakLS(type).le.0.) ebreakLS(type)=0.01
      if (ebreakCM(type).le.0.) ebreakCM(type)=0.01
c
c ebreakLS: LAB break-up energy
c ebreakCM: C.M. break-up energy
c E0n03: break-up energy
c E0n03CM: C.M. break-up energy
c w03    : width
c
      E0n03=ebreakLS(type)
      E0n03CM=ebreakCM(type)
      w03=width
c
      speccorCM(type)=0.
      speccorLS(type)=0.
c
c------------------    breakup nucleon energy   loop   -----------------
c
c ebegin               : first energy point of energy grid
c eend                 : last energy point of energy grid
c Eout                 : outgoing energy
c egrid                : outgoing energy
c EnCM                 : C.M. energy
c xspreeqbu            : nucleon breakup spectrum in Center of Mass System
c spec03CM             : nucleon breakup spectrum in Center of Mass System
c spec03LS             : nucleon breakup spectrum in Laboratory System
c speccorCM            : correction factors for breakup spectra
c speccorLS            : correction factors for breakup spectra
c sumtest              : check of breakup spectrum
c
c
      DO 501, nen=ebegin(type),eend(type)
c
c
        Eout=egrid(nen)   !ine breakup
        En=Eout
        EnCM=Eout
c
c      breakup nucleon spectra in Laboratory System
c
        if (Eout.le.EmaxnLS) then
          spec03LS(type,nen) = GAUSS(En,E0n03,w03)
c
c   BREAKUP THRESHOLD:  En>Emax BU
c
        else
          spec03LS(type,nen) = 0.D+00
        endif
c
c      breakup nucleon spectra in Center of Mass
c
        if  (Eout.le.EmaxnCM) then
          spec03CM(type,nen) = GAUSS(EnCM,E0n03CM,w03)
c
c   BREAKUP THRESHOLD: EnCM>Emax BU
c
        else
          spec03CM(type,nen) = 0.D+00
        endif
c
        speccorCM(type)=speccorCM(type)+spec03CM(type,nen)*deltaE(nen)
        speccorLS(type)=speccorLS(type)+spec03LS(type,nen)*deltaE(nen)
c
501   continue
c
c-------------------  END   breakup nucleon energy  loop ---------------
c
299   continue
c
c
301   continue
c
c----------------------     END   deuteron fragments      --------------
c
ccc            write(8,*)"#  breakup nucleon spectra:"
ccc            write(8,*)"#   Eout    spec_CM "
      do 505 type=1,2
        if(type.eq.1) then
ccc       write(8,*)"# breakup NEUTRON"
        else
ccc           write(8,*)"   "
ccc    write(8,*)"# breakup PROTON"
        endif
        if (speccorCM(type).eq.0..or.speccorLS(type).eq.0.) goto 505
        sumtest(type)=0.D+00
        do 504 nen=ebegin(type),eend(type)
          Eout=egrid(nen)
          spec03LS(type,nen) =
     1         xsBUnuc(type)*spec03LS(type,nen)/speccorLS(type)
          spec03CM(type,nen )=
     1         xsBUnuc(type)*spec03CM(type,nen)/speccorCM(type)
            sumtest(type) =
     1   sumtest(type)+spec03CM(type,nen)*deltaE(nen)
c
c   ******************************
      xspreeqbu(type,nen) = spec03CM(type,nen)
c   ******************************
ccc    if(Eout.le.EmaxnCM)then
ccc    write(8,889)Eout,xspreeqbu(type,nen)
ccc    endif
504   continue
505   continue
c
c
c
399   write(8,*)" incident deuteron energy too low"
      goto 699
c
499   write(8,*)" nucleon total breakup fraction zero"
      goto 699
c
699   continue
      write(8,*)'  end ..ine breakupAVR'
c
889     FORMAT(3x,f6.3,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,3x,f8.2,
     1         3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.2)
c
c
      return
      end
      subroutine checkBU(EnormEB,RnormEB,TnormBU)
c
c +---------------------------------------------------------------------
c | Author: Marilena Avrigeanu
c | Date  : August 10, 2015
c | Task  : check the elastic and nucleon breakup fraction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer IRDIM
      parameter (IRDIM=6001)
      integer Inorm,I
c
      real EnormEB,RnormEB,TnormBU,bind,arg,arg2,AT13
      real Edi(IRDIM),Tnorm(0:IRDIM),fBUT(IRDIM),fBUE(0:IRDIM)
      real FEDrat(IRDIM),fBUBF(IRDIM),fBU(IRDIM)
c
c   ***  Deuteron breakup fractions, Eqs.(1-3), PRC 88,014612(2013)  ***
c
c fBUE               : elastic breakup fraction
c IRDIM              : dimension fo breakup subroutine
c fBUBF              : nucleon inelastic breakup fraction
c fBUT               : nucleon breakup fraction
c fBU                : total breakup fraction
c RnormEB            : normalization factor for elastic breakup fraction
c EnormEB            : energy for elastic breakup fraction
c TnormBU            : normalization factor for total breakup fraction
c
c  According to the CDCC predictions for the elastic breakup cross
c  section behavior(PRC 82, 037601 (2010)), for higher energy than the
c  energetic domain (~30 MeV) where the parametrization was obtained,
c  the elastic breakup fraction is mentained constant, through the
c  RnormEB factor,    fBUE=RnormEB*fBUT.
c
c  Mainly for heavy nuclei, A~200, and  for Einc~Coulomb barrier, the
c  normalization constant TnormBu prevents the total breakup cross
c  section to exceed the reaction cross section.
c
c bind: energy bin
c I : counter
c Inorm: help variable
c
      AT13=Atarget**onethird
      RnormEB=1.
      TnormBU=1.
      bind=0.1
      Inorm=80./0.1 + 1.1
c
      do I=1,IRDIM
         Edi(I)=0.
         FEDrat(I)=0.
         Tnorm(I)=0.
         fBUT(I)=0.
         fBUE(I)=0.
         fBUBF(I)=0.
         fBU(I)=0.
      enddo
      Tnorm(0)=0.
      fBUE(0)=0.
c
c----------------------     deuteron incident energy   -----------------
c
      do 31, i=1,Inorm
c
c Edi: energy
c
      Edi(I)=i*bind
      arg=Edi(I)
      arg2=arg*arg
c
c      ********   PRC 89,044613 (2014)... parametrization        *****
c
        fBUT(I)=0.087-0.0066*Ztarget+0.00163*Ztarget*AT13+
     1          0.0017*AT13*arg-0.000002*Ztarget*arg2
         if(fBUT(I).le.0.d+00) then
            fBUT(I) = 0.
         go to 29
         endif
c
        fBUE(I)=0.031-0.0028*Ztarget+0.00051*Ztarget*AT13+
     1          0.0005*AT13*arg-0.000001*Ztarget*arg2
c
c FEDrat: ratio for breakup calculation
c
      FEDrat(I) = fBUE(I)/fBUT(I)
c
c      elastic breakup normalization for Ed>25 MeV
c
       if(fBUE(I).ge.fBUE(I-1))go to 25
         if((fBUE(I)/fBUE(I-1)).ge.0.9999) then
            INORM=I-1
            EnormEB=INORM*bind
            RnormEB=FEDrat(INORM)
         endif
c
       fBUE(I)=RnormEB*fBUT(I)
c
c     end elastic breakup normalization for Ed>25 MeV
c
25    continue
c
          fBUBF(i) = fBUT(i)-fBUE(i)
      fBU(i) = 2*fBUT(i)-fBUE(i)
c
c    total normalization for very heavy nuclei
c
      if(fBU(I).le.0.9D+00) go to 27
      Tnorm(I)=0.9D+00/fBU(I)
c
      if (Tnorm(I).lt.Tnorm(I-1)) TnormBU=Tnorm(I)
c
c      end total normalization for very heavy nuclei
c
  27  continue
  29  continue
  31  continue
c
c     ********     end PRC 89,044613 (2014)... parametrization    *******
c
c w0: width
c XG : term
c E0n : help variable
c
      return
      end
      FUNCTION GAUSS(En,E0n,w0)
      include "talys.cmb"
      real En,E0n,w0,term1,arg,XG,Gauss
      PI=ACOS(-1.)
      term1=1./(w0*sqrt(2.D+00*pi))
      arg=-(En-E0n)*(En-E0n)/(2.D+00*w0*w0)
      XG=term1*exp(-arg)
      GAUSS=term1*exp(arg)
      return
      end
