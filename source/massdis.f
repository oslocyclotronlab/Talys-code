      subroutine massdis
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 18, 2019
c | Task  : Fission fragment yields
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      include "gef.cmb"
      integer numZff,numNff
      parameter (numZff=80,numNff=150)
      character*13 fffile
      character*90 gefpath
      real         Exfis(1000),xsfis(1000),Jfis,fisepsA,fisepsB,
     +             partfisxs,sumpre,sumpost,Fmulti,beldm1(136,203),
     +             ushell1(136,203),xstabtot(numZff,numNff),
     +             xstabcomp(0:numZ,0:numN,numZff,numNff),
     +             popfpEx(numZff,numNff,0:numpop),
     +             popfpJ(numZff,numNff,0:numJ),Ebin(0:numpop),
     +             term,Etabtot(numZff,numNff,1000),partfisJ(0:numJ),
     +             Jtabtot(numZff,numNff,100),sumJ,sumxs,sumE,sum,
     +             xsfisFF,Ytabtot(numZff,numNff)
      integer      iz,ia,in,i,j,gefwrite,Zcomp,Ncomp,Z,A,odd,Zix,Nix,
     +             nexend,iskip,istep,nex,nen,type,iza,
     +             nexgef,Jgef,nb,parity
c
c ************************** Mass yields *******************************
c
c xsfistot   : total fission cross section
c nummass    : number of masses
c xsfpApre   : pre-neutron emission cross section
c xsfpApost  : post-neutron emission corrected cross section
c yieldApre  : pre-neutron emission fission yield
c yieldApost : post-neutron emission corrected fission yield
c numelem    : number of elements
c numZff: number of Z of fission fragments
c numNff: number of N of fission fragments
c fffile     : fission fragment file
c xsfpZApre  : pre-neutron emission isotopic cross section
c xsfpZApost : post-neutron emission corrected isotopic cross section
c yieldZApre : pre-neutron emission isotopic yield
c yieldZApost: post-neutron emission corrected isotopic yield
c fymodel    : fission yield model, 1: Brosa 2: GEF
c path       : directory containing structure files to be read
c gefpath    : path for GEF files
c popfpEx    : energy population of FF
c popfpJ     : spin population of FF
c flagoutfy  : flag for output detailed fission yield calculation
c gefwrite   : integer for output detailed fission yield calculation
c Rfiseps    : ratio for limit for fission cross section per nucleus
c fiseps     : limit for fission cross section per excitation energy bin
c
c Initialization
c
      do ia=1,nummass
        xsfpApre(ia)=0.
        xsfpApost(ia)=0.
        yieldApre(ia)=0.
        yieldApost(ia)=0.
        do type=0,numpar
          nuA(type,ia)=0.
        enddo
      enddo
      do in=1,numneu
        yieldNpre(in)=0.
        yieldNpost(in)=0.
        do iz=1,numelem
          xsfpZApre(iz,in)=0.
          xsfpZApost(iz,in)=0.
          yieldZApre(iz,in)=0.
          yieldZApost(iz,in)=0.
          do nex=0,1
            xsfpex(iz,in,nex)=0.
            yieldfpex(iz,in,nex)=0.
            fpratio(iz,in,nex)=0.
          enddo
        enddo
      enddo
      yieldtotpre=0.
      xsfptotpre=0.
      yieldtotpost=0.
      xsfptotpost=0.
      do type=0,6
        do i=1,numnu
          Pdisnu(type,i)=0.
        enddo
        nubar(type)=0.
      enddo
      do iz=1,numZff
        do in=1,numNff
          xstabtot(iz,in)=0.
          do nex=1,1000
            Etabtot(iz,in,nex)=0.
          enddo
          do J=1,100
            Jtabtot(iz,in,J)=0.
          enddo
          do nex=0,numpop
            popfpEx(iz,in,nex)=0.
          enddo
          do J=0,numJ
            popfpJ(iz,in,J)=0.
          enddo
          do Zcomp=0,maxZ
            do Ncomp=0,maxN
              xstabcomp(Zcomp,Ncomp,iz,in)=0.
            enddo
          enddo
        enddo
      enddo
      fpeps=Rfiseps*xsfistot
      if (fpeps.eq.0.) return
      if (fymodel.ge.2) then
        gefpath=trim(path)//'fission/gef/'
        if (flagoutfy) then
          gefwrite=1
        else
          gefwrite=0
        endif
c
c Read nuclear structure information for GEF
c
c beldm1: binding energy from liquid drop model
c fisepsA : fission tolerance
c fisepsB : fission tolerance
c ushell1: shell correction
c
        open (unit=4,file=trim(gefpath)//'beldm.dat',status='old')
        read(4,*) beldm1
        close(4)
        do i=1,203
          do j=1,136
           beldm(i,j)=beldm1(j,i)
          end do
        end do
        open (unit=4,file=trim(gefpath)//'ushell.dat',status='old')
        read(4,*) ushell1
        close(4)
        do i=1,203
          do j=1,136
           ushel(i,j)=ushell1(j,i)
          end do
        end do
        open (unit=4,file=trim(gefpath)//'nucprop.dat',status='old')
        do i=1,3885
          read(4,*) (RNucTab(i,j),j=1,8)
        end do
        close(4)
      endif
c
c Loop over nuclides
c
c Zcomp      : charge number index for compound nucleus
c maxZ       : maximal number of protons away from the initial
c              compound nucleus
c Ncomp      : neutron number index for compound nucleus
c maxN       : maximal number of neutrons away from the initial
c              compound nucleus
c ZZ,Z       : charge number of residual nucleus
c AA,A       : mass number of residual nucleus
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c maxex      : maximum excitation energy bin for compound nucleus
c iskip,istep: help variables
c xsbinary   : cross section from initial compound to residual nucleus
c Ex         : excitation energy
c partfisxs  : partial fission cross section
c fisfeedex  : fission contribution from excitation energy bin
c brosafy    : subroutine for fission fragment yields based on Brosa
c              model
c disa       : normalised fission fragment mass yield per excitation
c              energy bin
c disacor    : normalised fission product mass yield per excitation
c              energy bin
c disaz      : normalised fission fragment isotope yield
c              per excitation energy bin
c disazcor   : normalised fission product isotope yield
c              per excitation energy bin
c gefran     : number of random events for GEF calculation
c Exfis: excitation energy for fission
c
      do 20 Zcomp=0,maxZ
        do 25 Ncomp=0,maxN
          Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          odd=mod(A,2)
          Zix=Zindex(Zcomp,Ncomp,0)
          Nix=Nindex(Zcomp,Ncomp,0)
          if (xsfeed(Zcomp,Ncomp,-1).le.fpeps) goto 25
          if (Zcomp.eq.0.and.Ncomp.eq.0) then
            nexend=maxex(Zcomp,Ncomp)+1
          else
            nexend=maxex(Zcomp,Ncomp)
          endif
          fisepsA=fpeps/max(3*maxex(Zcomp,Ncomp),1)
          iskip=0
          istep=4
          if (fymodel.eq.2) then
            do 30 nex=1,1000
              Exfis(nex)=0.
              xsfis(nex)=0.
   30       continue
            nen=0
          endif
          do 40 nex=nexend,0,-1
            if (Zcomp.eq.0.and.Ncomp.eq.0.and.
     +        nex.eq.maxex(Zcomp,Ncomp)+1) then
              excfis=Etotal
              partfisxs=xsbinary(-1)
              do 42 J=0,numJ
                partfisJ(J)=0.
                do 44 parity=-1,1,2
                  partfisJ(J)=partfisJ(J)+
     +              fisfeedJP(Zcomp,Ncomp,nex,J,parity)
   44           continue
   42         continue
            else
              if (mod(iskip,istep).ne.0) then
                iskip=iskip+1
                goto 40
              endif
              if (nex-istep+1.lt.0) goto 40
              if (Ex(Zcomp,Ncomp,nex-istep+1).ge.30.) then
                partfisxs=0.
                do 45 J=0,numJ
                  partfisJ(J)=0.
   45           continue
                do 50 i=0,istep-1
                  partfisxs=partfisxs+fisfeedex(Zcomp,Ncomp,nex-i)
                  do 55 J=0,numJ
                    do 57 parity=-1,1,2
                      partfisJ(J)=partfisJ(J)+
     +                  fisfeedJP(Zcomp,Ncomp,nex-i,J,parity)
   57               continue
   55             continue
   50           continue
                if (partfisxs.ne.0) then
                  excfis=0.
                  do 60 i=0,istep-1
                    excfis=excfis+fisfeedex(Zcomp,Ncomp,nex-i)*
     +                Ex(Zcomp,Ncomp,nex-i)
   60             continue
                  excfis=excfis/partfisxs
                endif
                iskip=1
              else
                excfis=Ex(Zcomp,Ncomp,nex)
                partfisxs=fisfeedex(Zcomp,Ncomp,nex)
                do 62 J=0,numJ
                  partfisJ(J)=0.
                  do 64 parity=-1,1,2
                    partfisJ(J)=partfisJ(J)+
     +                fisfeedJP(Zcomp,Ncomp,nex,J,parity)
   64             continue
   62           continue
              endif
            endif
            if (partfisxs.gt.fisepsA) then
c
c Brosa
c Normalization: sum over disa = 2.
c
              if (fymodel.eq.1) then
                call brosafy(Zix,Nix)
                do 70 ia=1,A
                  xsfpApre(ia)=xsfpApre(ia)+0.5*disa(ia)*partfisxs
                  xsfpApost(ia)=xsfpApost(ia)+0.5*disacor(ia)*partfisxs
                  do 72 iz=1,Z
                    in=ia-iz
                    if (in.lt.1.or.in.gt.numneu) goto 72
                    xsfpZApre(iz,in)=xsfpZApre(iz,in)+
     +                0.5*disaz(ia,iz)*partfisxs
                    xsfpZApost(iz,in)=xsfpZApost(iz,in)+
     +                0.5*disazcor(ia,iz)*partfisxs
   72             continue
   70           continue
              endif
c
c GEF
c
c xsfis: fission cross section
c xsfisFF: fission cross section per FF
c
              if (fymodel.eq.2) then
                nen=nen+1
                Exfis(nen)=excfis
                xsfis(nen)=partfisxs
              endif
c
c GEF + TALYS evaporation
c Normalization: sum over Ytab, Etab, Jtab = 1
c
c flagfisout   : flag for output of fission information
c Jfis    : spin of fissioning system
c partfisJ : partial fission spin distribution
c xstabtot: total cross section from GEF
c Jtabtot: total spin from GEF
c xstabcomp: Z, N cross section from GEF
c Ytabtot: yield from GEF
c
              if (fymodel.eq.3.and.A.le.350) then
                fisepsB=fisepsA/(5*maxJ(Zcomp,Ncomp,nex))*0.5
                do 75 J=0,maxJ(Zcomp,Ncomp,nex)
                  if (partfisJ(J).lt.fisepsB) goto 75
                  Jfis=real(J)+0.5*odd
                  call gefsub(Z,A,excfis,Jfis)
                  write(*,*) " FF excitation for Z= ",
     +              Z," A= ",A," Ex= ",excfis," J= ",Jfis," xs= ",
     +              partfisJ(J)," N_cases:",N_cases
                  do 80 iza=1,N_cases
                    iz=NZMkey(iza,3)
                    in=NZMkey(iza,2)
                    if (iz.gt.numZff.or.in.gt.numNff) goto 80
                    term=Ytab(iza)*partfisJ(J)
                    xstabtot(iz,in)=xstabtot(iz,in)+term
                    xstabcomp(Zcomp,Ncomp,iz,in)=
     +                xstabcomp(Zcomp,Ncomp,iz,in)+term
                    do 85 nexgef=1,1000
                      Etabtot(iz,in,nexgef)=Etabtot(iz,in,nexgef)+
     +                  term*Etab(iza,nexgef)
   85               continue
                    do 87 Jgef=1,100
                      Jtabtot(iz,in,Jgef)=Jtabtot(iz,in,Jgef)+
     +                  term*Jtab(iza,Jgef)
   87               continue
   80             continue
   75           continue
              endif
            endif
   40     continue
c
c GEF
c Normalization: sum over ysum,yAz= 2. * sigma_fission
c
c flagffspin: flag to use spin distribution in initial population
c Fmulti: factor for multi-chance fission
c
          if (fymodel.eq.2.and.A.le.350) then
            call geftalys(real(Z),real(A),nen,Exfis,xsfis,gefwrite,
     +        gefran)
            do ia=1,A
              xsfpApre(ia)=xsfpApre(ia)+0.5*ysum(ia)
              xsfpApost(ia)=xsfpApost(ia)+0.5*ysump(ia)
              if (ia.le.200) then
                do iz=1,Z
                  in=ia-iz
                  if (in.ge.1.and.in.le.numneu) then
                    xsfpZApre(iz,in)=xsfpZApre(iz,in)+0.5*yAZ(ia,iz)
                    xsfpZApost(iz,in)=xsfpZApost(iz,in)+0.5*yAZp(ia,iz)
                  endif
                enddo
              endif
            enddo
            if (xsfistot.gt.0.) then
              Fmulti=Ncomp*xsfeed(Zcomp,Ncomp,-1)
              do i=1,numnu
                if (ann_sum(i).gt.0.)
     +            Pdisnu(1,i)=Pdisnu(1,i)+(Fmulti+ann_sum(i))/xsfistot
              enddo
              do ia=1,A
                if (anpre_sum(ia).gt.0.)
     +            nuA(1,ia)=nuA(1,ia)+(Fmulti+anpre_sum(ia))/xsfistot
              enddo
              nubar(1)=nubar(1)+(Fmulti+anMean_sum)/xsfistot
            endif
          endif
   25   continue
   20 continue
      do type=0,6
        sum=0.
        do 90 i=1,numin
          sum=sum+Pdisnu(type,i)
   90   continue
        if (sum.gt.0.) then
          do 95 i=1,numin
            Pdisnu(type,i)=Pdisnu(type,i)/sum
   95     continue
        endif
      enddo
c
c GEF + TALYS evaporation
c
c Ebin: energy of bin
c Etabtot: tabulated energy
c sumJ: sum over spin distribution
c Jgef: counter
c nexbeg: first energy index
c nexend: last energy index
c nexgef: counter
c
      if (fymodel.eq.3) then
        Ebin(0)=0.
        do i=1,numpop
          Ebin(i)=0.1*i
        enddo
        sumxs=0.
        do 131 iz=1,numZff
          do 132 in=1,numNff
            if (xstabtot(iz,in).eq.0.) goto 132
            sumxs=sumxs+xstabtot(iz,in)
            sumE=0.
            do nexgef=1,1000
              sumE=sumE+Etabtot(iz,in,nexgef)
            enddo
            if (sumE.gt.0.) then
              do nexgef=1,1000
                Etabtot(iz,in,nexgef)=Etabtot(iz,in,nexgef)/sumE
              enddo
            endif
            sumJ=0.
            do Jgef=1,100
              sumJ=sumJ+Jtabtot(iz,in,Jgef)
            enddo
            if (sumJ.gt.0.) then
              do Jgef=1,100
                Jtabtot(iz,in,Jgef)=Jtabtot(iz,in,Jgef)/sumJ
              enddo
            endif
  132     continue
  131   continue
        if (sumxs.gt.0.) then
          do iz=1,numZff
            do in=1,numNff
              Ytabtot(iz,in)=xstabtot(iz,in)/sumxs
            enddo
          enddo
        endif
        do iz=1,numZff
          do in=1,numNff
            ia=iz+in
            xsfisFF=xsfistot*Ytabtot(iz,in)
            xsfpApre(ia)=xsfpApre(ia)+xsfisFF
            xsfpZApre(iz,in)=xsfpZApre(iz,in)+xsfisFF
            do nexgef=1,1000
              popfpEx(iz,in,nexgef)=popfpEx(iz,in,nexgef)+
     +          xsfisFF*Etabtot(iz,in,nexgef)
            enddo
            do J=1,30
              popfpJ(iz,in,J)=Jtabtot(iz,in,J)
            enddo
          enddo
        enddo
        do 160 iz=1,Z
          do 170 ia=1,A
            in=ia-iz
            if (in.lt.1.or.in.gt.numneu) goto 170
            if (xsfpZApre(iz,in).eq.0.) goto 170
            sumJ=0.
            do J=1,30
              sumJ=sumJ+popfpJ(iz,in,J)
            enddo
            if (sumJ.gt.0.) then
              do J=1,30
                popfpJ(iz,in,J)=popfpJ(iz,in,J)/sumJ
              enddo
            endif
            nb=numpop
            do 165 nex=100,numpop
              if (popfpEx(iz,in,nex).lt.1.e-9) then
                nb=nex
                goto 167
              endif
  165       continue
  167       nb=max(nb,1)
            fffile='ff000000.ex'
            write(fffile(3:5),'(i3.3)') iz
            write(fffile(6:8),'(i3.3)') ia
            open (unit=1,file=fffile,status='replace')
            if (flagffspin) then
              write(1,*) nb+1,30,1," xs= ",xsfpZApre(iz,in)
              do nex=0,nb
                write(1,'(f10.5,30es12.5)') Ebin(nex),
     +            (0.5*popfpEx(iz,in,nex)*popfpJ(iz,in,J),
     +            J=1,30)
              enddo
            else
              write(1,*) nb+1,0,1," xs= ",xsfpZApre(iz,in)
              do nex=0,nb
                write(1,'(f10.5,es12.5)') Ebin(nex),popfpEx(iz,in,nex)
              enddo
            endif
            close(1)
  170     continue
  160   continue
      endif
c
c Normalization to fission yields (sum = 2)
c
c sumpre: sum over pre-neutron FP's
c
      sumpre=0.
      sumpost=0.
      do ia=1,Atarget
        sumpre=sumpre+xsfpApre(ia)
        sumpost=sumpost+xsfpApost(ia)
      enddo
      do 220 iz=1,Ztarget
        do 230 ia=iz+1,Atarget
          in=ia-iz
          if (in.gt.numneu) goto 230
          if (xsfpZApre(iz,in).eq.0.) goto 230
          if (sumpre.gt.0.) yieldZApre(iz,in)=2.*xsfpZApre(iz,in)/sumpre
          yieldApre(ia)=yieldApre(ia)+yieldZApre(iz,in)
          yieldNpre(in)=yieldNpre(in)+yieldZApre(iz,in)
          yieldtotpre=yieldtotpre+yieldZApre(iz,in)
          xsfptotpre=xsfptotpre+xsfpZApre(iz,in)
          if (fymodel.le.2) then
            if (sumpost.gt.0.) yieldZApost(iz,in)=
     +        2.*xsfpZApost(iz,in)/sumpost
            yieldApost(ia)=yieldApost(ia)+yieldZApost(iz,in)
            yieldNpost(in)=yieldNpost(in)+yieldZApost(iz,in)
            yieldtotpost=yieldtotpost+yieldZApost(iz,in)
            xsfptotpost=xsfptotpost+xsfpZApost(iz,in)
          endif
  230   continue
  220 continue
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely
