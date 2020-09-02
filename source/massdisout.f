      subroutine massdisout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 29, 2016
c | Task  : Output of fission fragment yields
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  ftype(2)
      character*12 isostring(0:1)
      character*90 yieldfile,fpfile
      integer      i,iz,ia,in,nen,nex
c
c ****************** Output of fission yields **************************
c
      isostring(0)='ground state'
      isostring(1)='isomer      '
      write(*,'(/" ++++++++++ FISSION YIELDS ++++++++++"/)')
      ftype(1)='A'
      ftype(2)='N'
      do i=1,2
c
c Write results to separate files
c
c ftype     : type of yield distribution
c yieldfile : file with fission yields
c natstring : string extension for file names
c iso       : counter for isotope
c Einc0     : incident energy in MeV
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c Starget   : symbol of target nucleus
c Ztarget   : charge number of target nucleus
c yieldApre : pre-neutron emission mass yield
c yieldApost: post-neutron emission corrected mass yield
c
        yieldfile='yield'//ftype(i)//'0000.000.fis'//natstring(iso)
        if (Einc0.lt.0.001) then
          write(yieldfile(7:13),'(es7.1)') Einc0
        else
          write(yieldfile(7:14),'(f8.3)') Einc0
          write(yieldfile(7:10),'(i4.4)') int(Einc0)
        endif
        open (unit=1,file=yieldfile,status='replace')
        write(1,'("# ",a1," + ",i3,a2,": fission yields")')
     +    parsym(k0),Atarget,Starget
        write(1,'("# E-incident = ",es12.5)') Einc0
        write(1,'("# ")')
        write(1,'("# ")')
        if (i.eq.1) then
          write(*,'(" Fission yields as function of A"/)')
          write(*,'("  ",a1,"    Pre-n yield   Post-n yield ",
     +      "      Pre-n xs      Post-n xs")') ftype(i)
          write(1,'("# ",a1,"    Pre-n yield   Post-n yield ",
     +      "      Pre-n xs      Post-n xs")') ftype(i)
          do 10 ia=1,Atarget
            write(*,'(i3,4es15.4)') ia,yieldApre(ia),yieldApost(ia),
     +        xsfpApre(ia),xsfpApost(ia)
            write(1,'(i3,4es15.4)') ia,yieldApre(ia),yieldApost(ia),
     +        xsfpApre(ia),xsfpApost(ia)
  10      continue
          write(*,'(/"Tot",4es15.4)') yieldtotpre,yieldtotpost,
     +      xsfptotpre,xsfptotpost
        else
          write(*,'(/" Fission yields as function of N")')
         write(*,'(/"  ",a1,"    Pre-n yield   Post-n yield")') ftype(i)
         write(1,'("# ",a1,"    Pre-n yield   Post-n yield")') ftype(i)
          do 20 in=1,Ntarget
            write(*,'(i3,2es15.4)') in,yieldNpre(in),yieldNpost(in)
            write(1,'(i3,2es15.4)') in,yieldNpre(in),yieldNpost(in)
  20      continue
          write(*,'(/"Tot",2es15.4)') yieldtotpre,yieldtotpost
        endif
        close (unit=1)
      enddo
c
c Write ff/fp production
c
c fpexist    : flag for existence of fission product
c fpfile     : file with fission product
c yieldZApre : pre-neutron emission isotopic yield
c yieldZApost: post-neutron emission corrected isotopic yield
c
      write(*,'(/" FF/FP production"/)')
      write(*,'("    Z    A iso Pre-n yield   Post-n yield",
     +  "      Pre-n xs      Post-n xs      Ratio")')
      do 210 ia=1,Atarget
        do 220 iz=1,Ztarget
          in=ia-iz
          if (in.lt.1.or.in.gt.Ninit) goto 220
          if (xsfpZApre(iz,in).lt.fpeps.and.
     +      xsfpZApost(iz,in).lt.fpeps.and..not.fpexist(iz,in))
     +      goto 220
          fpfile='fp000000.tot'//natstring(iso)
          write(fpfile(3:8),'(2i3.3)') iz,ia
          if (.not.fpexist(iz,in)) then
            open (unit=1,file=fpfile,status='replace')
            write(1,'("# ",a1," + ",i3,a2,": Fission product yield of ",
     +        i3,a2)') parsym(k0),Atarget,Starget,ia,nuc(iz)
            write(1,'("# ")')
            write(1,'("# # energies =",i6)') numinc
            write(1,'("# ")')
            write(1,'("# E-incident Pre-n yield Post-n yield ",
     +        "  Pre-n xs   Post-n xs")')
            do 230 nen=1,nin0-1
              write(1,'(5es12.5)') eninc(nen),0.,0.,0.,0.
  230       continue
          else
            open (unit=1,file=fpfile,status='old')
            do 240 nen=1,nin0+4
              read(1,*,end=250,err=250)
  240       continue
          endif
          if (xsfpZApre(iz,in).ge.fpeps.and.xsfpZApost(iz,in).ge.fpeps)
     +      write(*,'(2i5,4es15.4)') iz,ia,yieldZApre(iz,in),
     +      yieldZApost(iz,in),xsfpZApre(iz,in),xsfpZApost(iz,in)
          write(1,'(5es12.5)') Einc0,yieldZApre(iz,in),
     +      yieldZApost(iz,in),xsfpZApre(iz,in),xsfpZApost(iz,in)
  250     close (unit=1)
          if (xsfpex(iz,in,1).gt.0.) then
            do 260 nex=0,1
              write(fpfile(10:12),'("L",i2.2)') nex
              if (.not.fpexist(iz,in)) then
                open (unit=1,file=fpfile,status='replace')
                write(1,'("# ",a1," + ",i3,a2,": Fission product yield",
     +            " of ",i3,a2,1x,a12)') parsym(k0),Atarget,Starget,ia,
     +            nuc(iz),isostring(nex)
                write(1,'("# ")')
                write(1,'("# # energies =",i6)') numinc
                write(1,'("# ")')
                write(1,'("# E-incident Post-n yield  Post-n xs ",
     +            "  Ratio ")')
                do 270 nen=1,nin0-1
                  write(1,'(5es12.5)') eninc(nen),0.,0.,0.,0.
  270           continue
              else
                open (unit=1,file=fpfile,status='old')
                do 280 nen=1,nin0+4
                  read(1,*,end=290,err=290)
  280           continue
              endif
              if (xsfpZApre(iz,in).ge.fpeps.and.
     +          xsfpZApost(iz,in).ge.fpeps)
     +          write(*,'(2i5,i3,12x,es15.4,15x,2es15.4)') iz,ia,nex,
     +          yieldfpex(iz,in,nex),xsfpex(iz,in,nex),
     +          fpratio(iz,in,nex)
              write(1,'(4es12.5)') Einc0,yieldfpex(iz,in,nex),
     +          xsfpex(iz,in,nex),fpratio(iz,in,nex)
  290         close (unit=1)
  260       continue
          endif
          if (.not.fpexist(iz,in)) fpexist(iz,in)=.true.
  220   continue
c
c Write cumulative ff/fp production
c
        if (xsfpApre(ia).lt.fpeps.and.
     +    xsfpApost(ia).lt.fpeps.and..not.fpaexist(ia)) goto 210
        fpfile='fp000000.tot'//natstring(iso)
        write(fpfile(6:8),'(i3.3)') ia
        if (.not.fpaexist(ia)) then
          fpaexist(ia)=.true.
          open (unit=1,file=fpfile,status='replace')
          write(1,'("# ",a1," + ",i3,a2,": Fission product yield of A=",
     +      i3)') parsym(k0),Atarget,Starget,ia
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("# ")')
          write(1,'("# E-incident Pre-n yield Post-n yield",
     +      "   Pre-n xs   Post-n xs")')
          do 310 nen=1,nin0-1
            write(1,'(5es12.5)') eninc(nen),0.,0.,0.,0.
  310     continue
        else
          open (unit=1,file=fpfile,status='old')
          do 320 nen=1,nin0+4
            read(1,*,end=330,err=330)
  320     continue
        endif
        write(1,'(5es12.5)') Einc0,yieldApre(ia),yieldApost(ia),
     +    xsfpApre(ia),xsfpApost(ia)
  330   close (unit=1)
  210 continue
      write(*,'(/"Total     ",4es15.4)') yieldtotpre,yieldtotpost,
     +  xsfptotpre,xsfptotpost
      write(*,'(/" Cumulative FF/FP production"/)')
      write(*,'("    A    Pre-n yield   Post-n yield ",
     +  "      Pre-n xs      Post-n xs")')
      do 410 ia=1,Atarget
        write(*,'(i5,4es15.4)') ia,yieldApre(ia),yieldApost(ia),
     +    xsfpApre(ia),xsfpApost(ia)
  410 continue
      write(*,'(/"Total",4es15.4)') yieldtotpre,yieldtotpost,
     +  xsfptotpre,xsfptotpost
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
