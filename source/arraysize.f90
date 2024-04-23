subroutine arraysize
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of size of largest arrays
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2023-03-10: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! *** Declaration of local data
!
! *************************** Array size *******************************
!
  numZ1=numZ+1
  numN1=numN+1
  numZph1=numZph+1
  numNph1=numNph+1
  numZchan1=numZchan+1
  numNchan1=numNchan+1
  numZastro1=numZastro+1
  numNastro1=numNastro+1
  numlev1=numlev+1
  numpar1=numpar+1
  numex1=numex+1
  numen1=numen+1
  numin1=numin+1
  numip1=numip+1
  numid1=numid+1
  numit1=numit+1
  numih1=numih+1
  numia1=numia+1
  numchantot1=numchantot+1
  numgamqrpa1=numgamqrpa+1
  numJ1=numJ+1
  numdens1=numdens+1
  numbar1=numbar+1
  numexc1=numexc+1
  nummt1=nummt+1
  numisom1=numisom+1
  numen61=numen6+1
  numen21=numen2+1
  numparx1=numexc/2+1
  numenlow1=numenlow+1
  numenin1=numenin+1
  numenrp1=numenrp+1
  numangcont1=numangcont+1
  numl1=numl+1
  numJph1=numJph+1
  numenrec1=numenrec+1
  numangrec1=numangrec+1
  write(*, '(/"  Size of some large arrays"/)')
  write(*, '(" gamexist",t20,i9)') numZ1*numN1*numlev1*numlev1
  write(*, '(" chanisoexist",t20,i9)') numin1*numip1*numid1*numit1*numih1*numia1*numlev1
  write(*, '(" specexcl",t20,i9)') numchantot1*numpar1*numex1*numen1
  write(*, '(" feedexcl",t20,i9)') numZchan1*numNchan1*numpar1*numex1*numex1
  write(*, '(" fqrpa",t20,i9)') numZ1*numN1*numgamqrpa1*2*numgam
  write(*, '(" bassign",t20,i9)') numZ1*numN1*numlev1*numlev1
  write(*, '(" ldtable",t20,i9)') numZ1*numN1*numdens1*numJ1*3*numbar
  write(*, '(" phtable2",t20,i9)') 2*2*numexc1*numexc1*numexc1*numexc1*numdens1
  write(*, '(" xspopph2",t20,i9)') numZph1*numNph1*numex1*numparx1*numparx1*numparx1*numparx1
  write(*, '(" frescue",t20,i9)') nummt1*numisom1*numen61
  write(*, '(" preeqpop",t20,i9)') numZ1*numN1*numex1*numJ1*3
  write(*, '(" xspop",t20,i9)') numZ1*numN1*numex1*numJ1*3
  write(*, '(" wemission2",t20,i9)') numpar1*numparx1*numparx1*numparx1*numparx1*numen1
  write(*, '(" xspopph2",t20,i9)') numZph1*numNph1*numparx1*numparx1*numparx1*numparx1*numex1
  write(*, '(" xsrp",t20,i9)') numZ1*numN1*numisom1*numenrp1
  write(*, '(" Tjlnex",t20,i9)') numpar1*numex1*3*numl1
  write(*, '(" fxsgamdischan",t20,i9)') numenlow1*numchantot1*numlev*numlev
  write(*, '(" areaejlab",t20,i9)') numpar1*numen21*2*numangcont1
  write(*, '(" areareclab",t20,i9)') numZ1*numN1*numenrec1*2*numangrec1
  write(*, '(" xsastroex",t20,i9)') numZastro1*numNastro1*numenin1*numlev1
  write(*, '(" ddxrec",t20,i9)') numZ1*numN1*numex1*numenrec1*2*numangrec1
  write(*, '(" sfactor",t20,i9)') numZ1*numN1*numex1*numJ1*3
  write(*, '(" Tjlnex",t20,i9)') numpar1*numex1*3*numl1
  write(*, '(" phdensjp",t20,i9)') numZ1*numN1*numdens1*numJph1*3
  write(*, '()')
  return
end subroutine arraysize
! Copyright A.J. Koning 2023
