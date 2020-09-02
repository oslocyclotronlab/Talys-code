      program talys
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning, Stephane Hilaire and Stephane Goriely
c | Date  : December 28, 2019
c | Task  : Main program
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c                            TALYS-1.95
c
c                   (Version: December 28, 2019)
c
c             Nuclear reactions from 1 keV to 200 MeV.
c
c Consult the manual for a detailed outline of TALYS.
c Please send any comments, corrections or improvements to one of
c the three principal authors:
c
c       |-------------------------------------------------------|
c       |                 Arjan Koning                          |
c       |                                                       |
c       | IAEA - International Atomic Energy Agency             |
c       | Wagrammerstrasse 5                                    |
c       | P.O. Box 100, 1400 Vienna Austria                     |
c       | Email: A.koning@iaea.org                              |
c       |-------------------------------------------------------|
c
c       |-------------------------------------------------------|
c       |                 Stephane Hilaire                      |
c       |                                                       |
c       | Commissariat a l'Energie Atomique, Bruyeres-le-Chatel |
c       | Service de Physique et Techniques Nucleaires          |
c       | B.P. 12, F-91680 Bruyeres-le-Chatel, France           |
c       | Phone: (+33) 1 6926 6028 FAX: (+33) 1 6926 7063       |
c       | Email: stephane.hilaire@cea.fr                        |
c       |-------------------------------------------------------|
c
c       |-------------------------------------------------------|
c       |                 Stephane Goriely                      |
c       |                                                       |
c       | Institut d'Astronomie et d'Astrophysique              |
c       | Universite Libre de Bruxelles                         |
c       | Campus de la Plaine, CP-226, 1050 Brussels, Belgium   |
c       | Phone: (+32) 2 650 2843                               |
c       | Email: sgoriely@astro.ulb.ac.be                       |
c       |-------------------------------------------------------|
c
c ************* Input, initialization and reaction models **************
c
c machine      : subroutine for machine dependent statements
c constants    : subroutine for constants and initialization
c talysinput   : subroutine for user input and defaults
c talysinitial : subroutine for initialization of nuclear structure
c talysreaction: subroutine with reaction models
c flagnatural  : flag for calculation of natural element
c natural      : subroutine for calculation of natural element
c
      call machine
      call constants
      call talysinput
      call talysinitial
      call talysreaction
      if (flagnatural) call natural
      end
Copyright (C)  2019 A.J. Koning, S. Hilaire and S. Goriely
