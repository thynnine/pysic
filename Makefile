#
# Makefile for Pysic
#

# fortran compiler
FORTRAN       = gfortran
# mpi fortran compiler
MPIFORTRAN    = mpif90

# build directory
BUILDDIR      = build
TMPDIR        = fortran_tmp

# compiler flags
DEBUG_FLAGS   = "-g -fbounds-check"
F90_FLAGS     = 

# module name suffix
DEBUG_SUFFIX  = _debug
NOOPT_SUFFIX  = _noopt
SERIAL_SUFFIX = _serial

# source files, do not edit
FORTRAN_FILES = pysic_fortran.pyf Mersenne.F90 MPI.F90 Quaternions.F90 Utility.F90 Geometry.F90 Potentials.F90 Core.F90 PyInterface.F90

.PHONY: help clean test debug noopt serial parallel

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  help      print this message"
	@echo "  clean     trash previous build"
	@echo "  debug     compile debug version"
	@echo "  noopt     compile unoptimized version"
	@echo "  serial    compile serial version"
	@echo "  parallel  compile parallel version"

test:
	mkdir -p $(TMPDIR)	
	cp ./fortran/MPI.f90 $(TMPDIR)/MPI.F90
	cp ./fortran/Quaternions.f90 $(TMPDIR)/Quaternions.F90
	cp ./fortran/Utility.f90 $(TMPDIR)/Utility.F90
	cp ./fortran/Potentials.f90 $(TMPDIR)/Potentials.F90
	cp ./fortran/Core.f90 $(TMPDIR)/Core.F90
	cp ./fortran/Geometry.f90 $(TMPDIR)/Geometry.F90
	cp ./fortran/PyInterface.f90 $(TMPDIR)/PyInterface.F90
	cp ./fortran/Mersenne.F90 $(TMPDIR)/Mersenne.F90

	cd $(TMPDIR); f2py -m pysic_fortran -h pysic_fortran.pyf PyInterface.F90; \
	f2py -c --noopt --fcompiler=$(FORTRAN) --f90exec=$(FORTRAN) $(FORTRAN_FILES)
	rm $(TMPDIR)/*
	rmdir $(TMPDIR)

	@echo
	@echo "Test build finished."


clean:
	-rm -rf $(TMPDIR)
	-rm -rf $(BUILDDIR)

debug:
	mkdir -p $(BUILDDIR)
	mkdir $(BUILDDIR)/pysic$(DEBUG_SUFFIX)
	cp -r ./pysic/* $(BUILDDIR)/pysic$(DEBUG_SUFFIX)

	mkdir $(TMPDIR)	
	cp ./fortran/MPI.f90 $(TMPDIR)/MPI.F90
	cp ./fortran/Quaternions.f90 $(TMPDIR)/Quaternions.F90
	cp ./fortran/Utility.f90 $(TMPDIR)/Utility.F90
	cp ./fortran/Potentials.f90 $(TMPDIR)/Potentials.F90
	cp ./fortran/Core.f90 $(TMPDIR)/Core.F90
	cp ./fortran/Geometry.f90 $(TMPDIR)/Geometry.F90
	cp ./fortran/PyInterface.f90 $(TMPDIR)/PyInterface.F90
	cp ./fortran/Mersenne.F90 $(TMPDIR)/Mersenne.F90

	cd $(TMPDIR); f2py --debug-capi -m pysic_fortran -h pysic_fortran.pyf PyInterface.F90; \
	f2py --debug-capi -c --debug --noopt --fcompiler=$(FORTRAN) --f90exec=$(FORTRAN) --f90flags=$(DEBUG_FLAGS) $(FORTRAN_FILES)
	mv $(TMPDIR)/pysic_fortran.so $(BUILDDIR)/pysic$(DEBUG_SUFFIX)
	rm $(TMPDIR)/*
	rmdir $(TMPDIR)

	@echo
	@echo "Build finished. The module pysic$(DEBUG_SUFFIX) is in the directory $(BUILDDIR)."

noopt:
	mkdir -p $(BUILDDIR)
	mkdir $(BUILDDIR)/pysic$(NOOPT_SUFFIX)
	cp -r ./pysic/* $(BUILDDIR)/pysic$(NOOPT_SUFFIX)

	mkdir -p $(TMPDIR)	
	cp ./fortran/MPI.f90 $(TMPDIR)/MPI.F90
	cp ./fortran/Quaternions.f90 $(TMPDIR)/Quaternions.F90
	cp ./fortran/Utility.f90 $(TMPDIR)/Utility.F90
	cp ./fortran/Potentials.f90 $(TMPDIR)/Potentials.F90
	cp ./fortran/Core.f90 $(TMPDIR)/Core.F90
	cp ./fortran/Geometry.f90 $(TMPDIR)/Geometry.F90
	cp ./fortran/PyInterface.f90 $(TMPDIR)/PyInterface.F90
	cp ./fortran/Mersenne.F90 $(TMPDIR)/Mersenne.F90

	cd $(TMPDIR); f2py -m pysic_fortran -h pysic_fortran.pyf PyInterface.F90; \
	f2py -c --noopt --fcompiler=$(FORTRAN) --f90exec=$(FORTRAN) $(FORTRAN_FILES)
	mv $(TMPDIR)/pysic_fortran.so $(BUILDDIR)/pysic$(NOOPT_SUFFIX)
	rm $(TMPDIR)/*
	rmdir $(TMPDIR)

	@echo
	@echo "Build finished. The module pysic$(NOOPT_SUFFIX) is in the directory $(BUILDDIR)."

serial:
	mkdir -p $(BUILDDIR)
	mkdir $(BUILDDIR)/pysic$(SERIAL_SUFFIX)
	cp -r ./pysic/* $(BUILDDIR)/pysic$(SERIAL_SUFFIX)

	mkdir -p $(TMPDIR)	
	cp ./fortran/MPI.f90 $(TMPDIR)/MPI.F90
	cp ./fortran/Quaternions.f90 $(TMPDIR)/Quaternions.F90
	cp ./fortran/Utility.f90 $(TMPDIR)/Utility.F90
	cp ./fortran/Potentials.f90 $(TMPDIR)/Potentials.F90
	cp ./fortran/Core.f90 $(TMPDIR)/Core.F90
	cp ./fortran/Geometry.f90 $(TMPDIR)/Geometry.F90
	cp ./fortran/PyInterface.f90 $(TMPDIR)/PyInterface.F90
	cp ./fortran/Mersenne.F90 $(TMPDIR)/Mersenne.F90

	cd $(TMPDIR); f2py -m pysic_fortran -h pysic_fortran.pyf PyInterface.F90; \
	f2py -c --noopt --fcompiler=$(FORTRAN) --f90exec=$(FORTRAN) $(FORTRAN_FILES)
	mv $(TMPDIR)/pysic_fortran.so $(BUILDDIR)/pysic$(SERIAL_SUFFIX)
	rm $(TMPDIR)/*
	rmdir $(TMPDIR)

	@echo
	@echo "Build finished. The module pysic$(SERIAL_SUFFIX) is in the directory $(BUILDDIR)."


parallel:
	mkdir -p $(BUILDDIR)
	mkdir $(BUILDDIR)/pysic
	cp -r ./pysic/* $(BUILDDIR)/pysic

	mkdir -p $(TMPDIR)	
	cp ./fortran/MPI.f90 $(TMPDIR)/MPI.F90
	cp ./fortran/Quaternions.f90 $(TMPDIR)/Quaternions.F90
	cp ./fortran/Utility.f90 $(TMPDIR)/Utility.F90
	cp ./fortran/Potentials.f90 $(TMPDIR)/Potentials.F90
	cp ./fortran/Core.f90 $(TMPDIR)/Core.F90
	cp ./fortran/Geometry.f90 $(TMPDIR)/Geometry.F90
	cp ./fortran/PyInterface.f90 $(TMPDIR)/PyInterface.F90
	cp ./fortran/Mersenne.F90 $(TMPDIR)/Mersenne.F90

	cd $(TMPDIR); f2py -m pysic_fortran -h pysic_fortran.pyf PyInterface.F90; \
	f2py -c --fcompiler=$(FORTRAN) --f90exec=$(MPIFORTRAN) --f90flags="-D MPI $(F90_FLAGS)" $(FORTRAN_FILES)
	mv $(TMPDIR)/pysic_fortran.so $(BUILDDIR)/pysic
	rm $(TMPDIR)/*
	rmdir $(TMPDIR)

	@echo
	@echo "Build finished. The module pysic is in the directory $(BUILDDIR)."
