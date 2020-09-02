      subroutine buenhance(Z00,N00,BFsum)
c +---------------------------------------------------------------------
c
c | Author: Marilena Avrigeanu
c | Date  : February 2015
c | Task  : Inelastic breakup enhancement brought by breakup neutrons
c           and protons to any (Z,N) residual nucleus from d+Atarget
c           interaction, PRC 89,044613 (2014) and References therein
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,Z,N,nen,Zrez,Nrez,Z00,N00
      real BFenhance(2),enhance03(2),conv03(2),Ebreak(2),BCout(2)
      real sumtest0(2,numen),sumtest(2),sumGauss(2)
      real Emaxn,E0n03,term03,Gauss,width,w03,En,BCin,Bdeut,BFsum
c
c ************************** Avrigeanu model ***************************
c
c DEUTERON Break-up model by M.Avrigeanu et al., PRC89,044613(2014),
c                                                NDS118,301(2014)
c
c Einc        : incident energy in MeV
c Atarget     : mass number of target nucleus
c BCin,BCout  : effective Coulomb barrier
c Ztarget     : charge number of target nucleus
c Bdeut       : deuteon binding energy
c Emaxn       : maximum breakup nucleons energy in Laboratory System
c Ebreak      : Centroid energy of the breakup nucleon energy
c               distributions, Kalbach 2003, in Laboratory System
c width       : Full width at half maximum of the breakup nucleon
c               energy distribution, Kalbach 2003
c Z00         : charge number for residual nucleus
c N00         : neutron number for residual nucleus
c Z           : charge number for residual nucleus
c N           : neutron number for residual nucleus
c ENHratio    : breakup nucleons enhancing reaction cross
c               sections, PRC 89,044613 Eq. (7),
c               n + Atarget  sig(n,Z,A,Eout)/sig_Total(n,Enout);
c               p + Atarget  sig(p,Z,A,Eout)/sig_Reaction(p,Epout)
c BFenhance   : inelastic breakup enhancement brought by breakup
c               nucleons to any (Z,N) residual nucleus from d+Atarget
c               interaction
c BFsum       : total inelastic breakup enhancement to (Z,N) residual
c               nucleus population
c
      BFsum = 0.
c
      do type=1,2
        Ebreak(type)=0.
        conv03(type)=0.
        enhance03(type)=0.
        sumtest(type)=0.
        sumGauss(type)=0.
        BFenhance(type)=0.
      enddo
c
      do 103 type=1,2
      do 103 nen=1,eend(type)
            sumtest0(type,nen)=0.
103      continue
c
      Z=Z00
      N=N00
      Bdeut = 2.225
      BCin=Ztarget/9.5
      BCout(1)=0.
      BCout(2)=BCin
       Emaxn = Einc*(Atarget+1.)/(Atarget+2.)-Bdeut*(Atarget+1.)/Atarget
c
c      Breakup threshold
c
      if(Emaxn.le.0.d+00) then
            Emaxn=0.
            write(8,*)" deuteron energy lower than Bd, no breakup"
            go to 499
      endif
c
c
cBUenhance
cBUenhance      In the following loops, 109 and 115 , in order to check
CBUenhance      the implementation of the inelastic breakup enhancement
cBUenhance     in TALYS, arbitrary values are given to ENHratios, until
cBUenhance      the corresponding external file/files will be created
cBUenhance         from TENDL, and then properly read in subroutine
cBUenhance      BUiniial.f (still commented its call in preeqcomplex.f,
cBUenhance      line 41))
cBUenhance
c
c Zrez: counter
c Nrez: counter
c
      do 109 Zrez=0,numZ
            do 109 Nrez=0,numN
                  do 107 nen=0,numen
                  ENHratio(1,Zrez,Nrez,nen)=0.02
107        continue
109      continue
c
      do 115 Zrez=0,numZ
            do 115 Nrez=0,numN
                  do 113 nen=0,numen
                        ENHratio(2,Zrez,Nrez,nen)=0.01
113        continue
115      continue
c
c
c----------------------     fragnent type  -----------------------------
c
      do 401, type=1,2
c
c
      Ebreak(type) = 0.5*(Atarget+1.)/(Atarget+2.)*Einc+
     1      0.5*(Atarget+1.)/Atarget*(-Bdeut-BCin+2.*BCout(type))
      width = 1.15 + 0.12*Einc - Atarget/140.
C
      if(width.le.0.d+00) then
            width=0.
            go to 499
      endif
c
      if(Ebreak(type).le.0.d+00) then
      Ebreak(type)=0.01
      endif
c
      E0n03=Ebreak(type)
      w03=width
c
      conv03(type)     = 0.
      sumGauss(type)   = 0.
      BFenhance(type)  = 0.
c
c----------------------     breakup nucleon energy   -------------------
c
c xsBF                         : nucleon inelastic breakup cross section
c                                calculated in subroutine breakupAVR
c ebegin                       : first energy point of energy grid
c eend                         : last energy point of energy grid
c En,egrid                     : outgoing energy
c term03                       : help variable for  breakup enhancement
c                                calculations
c conv03               : help variable for  breakup enhancement
c                                calculations
c enhance03    : help variable for  breakup enhancement
c                                calculations
c sumtest0                     : help variables for ckhecking the
c                                calculation procedures
c sumtest                      : help variables for ckhecking the
c                                calculation procedures
c sumGauss                     : help variables for ckhecking the
c                                calculation procedures
c
      do 301, nen=ebegin(type),eend(type)
            En=egrid(nen)
         if(En.le.Emaxn) then
            sumtest0(type,nen) = GAUSS(En,E0n03,w03)
         term03 = ENHratio(type,Z,N,nen)*GAUSS(En,E0n03,w03)*deltaE(nen)
c
c      BREAKUP THRESHOLD:  En>Emax BU
c
         else
            sumtest0(type,nen) = 0.
            term03 = 0.
          endif
c
            sumGauss(type)=sumGauss(type)+sumtest0(type,nen)*deltaE(nen)
            conv03(type) = conv03(type) + term03
301      continue
c
c--------------------   END   breakup nucleon energy   -----------------
c
            enhance03(type) = xsBF(type)*conv03(type)
            BFenhance(type) = enhance03(type)/sumGauss(type)
            BFsum= BFsum + BFenhance(TYPE)
c
401      continue
c
c----------------------     END   fragnent type  -----------------------
c
c
499      continue
c
c
             write(8,*)'   '
        write(8,*)"   Ed    type    BFenhance    Emaxn    ",
     1  " E0_SL     w03   sumtest     xsBUnuc"
c
c
            do 506 type=1,2
            if(sumGauss(type).eq.0.d+00)cycle
              do 505 nen=ebegin(type),eend(type)
                  sumtest0(type,nen)=xsBUnuc(type)*
     1              sumtest0(type,nen)/sumGauss(type)
              sumtest(type)=sumtest(type)+sumtest0(type,nen)*deltaE(nen)
505              continue
c
        write(8,999)Einc,type,BFenhance(type),Emaxn,Ebreak(type),
     1                    width,sumtest(type),xsBUnuc(type)

c
c
506      continue
c
          write(8,*)' end .. ine BUenhance:  Z =',Z,' N=',N,' BFsum=',
     1                                    BFsum
c
999      Format(2x,f6.2,1x,I3,5x,f8.3,4x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,
     1       4x,f8.3,3x,f8.3,4x,f8.3,4x,f8.3,4x,f8.3,4x,f8.3)
c
c
      return
      end
