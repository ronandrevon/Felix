#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# felixsim
#
# Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
#
# (C) 2013/14, all rights reserved
#
# Version: :VERSION:
# Date:    :DATE:
# Time:    :TIME:
# Status:  :RLSTATUS:
# Build:   :BUILD:
# Author:  :AUTHOR:
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#  This file is part of felixsim.
#
#  felixsim is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  felixsim is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# $Id: makefile-CSC.GF,v 1.2 2014/03/27 18:13:30 phsht Exp $
# cygwin needs: gcc, gfortran, make, git, openmpi, lapack, fftw3
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#To run in windows command line (cmd), set up environmental variables for mpi, mkl and ifort:
#mpi
#in: $C:\Program Files (x86)\IntelSWTools\mpi\5.1.1.110\intel64\bin
#'CALL mpivars.bat intel64'
#
#ifort
#in: C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016\windows\bin
#'CALL ifortvars.bat intel64'
#
#mkl:
#in: C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016\windows\mkl\bin\intel64
# 'CALL mklvars_intel64'

F90 = mpiifort
FC = ifort
#F90FLAGS =
#F77FLAGS =

#links to BLAS, LAPACK and fftw3 through mkl
INCLUDEDIR = "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016\windows\mkl\include"

#no idea what these flags do - got them through the mkl_link_tool:
#C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016\windows\mkl\tools
#Qmkl seemed to be the automatic way of setting up the mkl libraries
LIBFLAGS= mkl_intel_lp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib /Qopenmp /Qmkl


SIMOBJECTFILES=gmodules.obj smodules.obj inputmodules.obj errorcodes.obj message.obj felixsim.obj util.obj in.obj out.obj readinput.obj \
readciffile.obj experimentalsetup.obj imagesetup.obj montagesetup.obj refineutils.obj\
structurefactorsetup.obj ciftbx.obj hash_funcs.obj eigen.obj invert.obj microscopy.obj \
diffractionpatterndefinitions.obj diffractionpatterninitialisation.obj crystallographyinitialisation.obj \
crystallography.obj bloch.obj Ug.obj image.obj errorchecks.obj writeoutput.obj dqng.obj xerror.obj d1mach.obj

REFINEOBJECTFILES=gmodules.obj smodules.obj inputmodules.obj message.obj errorcodes.obj felixrefine.obj util.obj in.obj out.obj readciffile.obj \
experimentalsetup.obj imagesetup.obj montagesetup.obj structurefactorsetup.obj readinput.obj\
eigen.obj invert.obj microscopy.obj diffractionpatterndefinitions.obj diffractionpatterninitialisation.obj \
crystallographyinitialisation.obj crystallography.obj bloch.obj Ug.obj image.obj ciftbx.obj simplex.obj\
hash_funcs.obj RefineWriteOut.obj felixfunction.obj refineutils.obj writeoutput.obj dqng.obj xerror.obj d1mach.obj errorchecks.obj\
symmetry.obj

%.obj: %.f90
	$(F90) $(F90FLAGS) -c $< -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -c $< $(INCLUDEDIR) $(LIBFLAGS)

%.obj: %.F90
	$(F90) $(F90FLAGS) -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.obj: %.f
	$(FC) $(F77FLAGS) -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.obj: %.F
	$(FC) $(F77FLAGS) -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

all:	felixsim felixrefine

felixsim: $(SIMOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(SIMOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(SIMOBJECTFILES) $(LIBFLAGS)

felixrefine: $(REFINEOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(REFINEOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(REFINEOBJECTFILES) $(LIBFLAGS)


print:	
	a2ps -E -o print.ps gmodules.f90 smodules.f90 message.f90 errorcodes.f90 felixsim.f90 util.f90 \
in.f90 out.f90 readinput.f90 readciffile.f90 experimentalsetup.f90 imagesetup.f90 \
montagesetup.f90 structurefactorsetup.f90 eigen.f90 invert.f90 microscopy.f90 \
diffraction.f90 diffractionpatterninitialisation.f90 crystallographyinitialisation.f90 \
crystallography.f90 bloch.f90 Ug.f90 image.f90 ciftbx.f hash_funcs.f \
writeoutput.f90 dqng.f xerror.f d1mach.f felixdraw.f90 felixrefine.f90 refineutils.f90 \
makefile-cygwin.GF; convert -density 150 print.ps print.pdf

#awful windows command line syntax, same comand as 'rm -f'
clean:	
	DEL /F /S /Q /A "./*.mod"
	DEL /F /S /Q /A "*.obj"
#clean all:	
#	rm -f core *.mod *.obj *.exe
