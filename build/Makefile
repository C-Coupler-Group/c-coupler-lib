# Common Makefile: a framework for building all libraries
#
# This Makefile was originally from CCSM, and then modified by Dr. Cheng Zhang
# and Dr. Li Liu (liuli-cess@tsinghua.edu.cn)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# parse cmd-line and establish values for EXEC, VPATH, SRCS, OBJS, etc
#-------------------------------------------------------------------------------

EXEC     := libc_coupler.a
VPFILE   := srcfiles
VPATH    := .
SRCFILE  := NONE
SRCS     := NONE
DEPGEN   := ./makdep  # an externally provided dependency generator
MOD_SUFFIX := mod

ifneq ($(VPATH),.)
  # this variable was specified on cmd line or in an env var
else
  ifneq ($(VPFILE),NONE)
    # explicit list of VPATH dirs is provided
    VPATH := $(wildcard . $(shell cat $(VPFILE) ) )
  endif
endif

ifneq ($(SRCS),NONE)
  # this variable was specified on cmd line or in an env var
else
  ifneq ($(SRCFILE),NONE)
    # explicit list of src files is provided
    SRCS := $(shell cat $(SRCFILE) )
  else
    # list of src files is all .F90 .F .c files in VPATH
    SRCS := $(wildcard $(addsuffix /*.cxx , $(VPATH)) \
                     $(addsuffix /*.F90, $(VPATH)) \
		       $(addsuffix /*.[F], $(VPATH)) )
  endif
endif

OBJS  := $(addsuffix .o, $(sort $(basename $(notdir $(SRCS)))))
DEPS  := $(addsuffix .d, $(sort $(basename $(notdir $(SRCS)))))
INCS  := $(patsubst %,-I%, $(VPATH) )
RM    := rm

.SUFFIXES:
.SUFFIXES: .cxx .F90 .F .o

all: $(EXEC)

$(DEPGEN) : makdep.c
	$(CC) -o $(DEPGEN) $<
#-------------------------------------------------------------------------------
# include the file that provides macro definitions required by build rules
# note: the MACFILE may not be needed for certain goals
#-------------------------------------------------------------------------------

ifneq ($(MAKECMDGOALS), db_files)
  -include $(MACFILE)
  ifeq ($(USE_OMP),TRUE)
        CFLAGS      += -openmp
        CXXFLAGS    += -openmp
        FFLAGS      += -openmp
        LDFLAGS     += -openmp
  endif
endif

#-------------------------------------------------------------------------------
# echo file names, paths, compile flags, etc. used during build
#-------------------------------------------------------------------------------

db_files:
	@echo " "
	@echo "* EXEC    := $(EXEC)"
	@echo "* MACFILE := $(MACFILE)"
	@echo "* VPFILE  := $(VPFILE)"
	@echo "* VPATH   := $(VPATH)"
	@echo "* SRCFILE := $(SRCFILE)"
	@echo "* INCS    := $(INCS)"
	@echo "* SRCS    := $(SRCS)"
	@echo "* OBJS    := $(OBJS)"
	@echo "* DEPS    := $(DEPS)"
db_flags:
	@echo " "
	@echo "* cpp     := $(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR)"
	@echo "* cc      := $(CC)  $(CFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F.o    := $(FC)  $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F90.o  := $(FC)  $(FFLAGS) $(FREEFLAGS) $(INCS) $(INCLDIR)"

#-------------------------------------------------------------------------------
# build rules: MACFILE, cmd-line, or env vars must provide the needed macros
#-------------------------------------------------------------------------------

$(EXEC): $(OBJS)
	$(AR) $(ARFLAGS) $(EXEC) $(OBJS)

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCS) $(CPPDEFS) $(INCLDIR) $<

#.c.o:
#	$(CC) -c $(INCLDIR) $(INCS) $(CPPDEFS) $(CFLAGS)  $<

ifeq ($(CPP),NONE)
.f.o:
	$(FC) $(INCLDIR) $(INCS) $(FFLAGS) $(FIXEDFLAGS) $(CPPFLAGS) $(CPPDEFS) $<

.f90.o:
	$(FC) $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS)  $(CPPFLAGS) $(CPPDEFS) $<

.F.o:
	$(FC) $(INCLDIR) $(INCS) $(FFLAGS) $(FIXEDFLAGS) $(CPPFLAGS) $(CPPDEFS) $<

.F90.o:
	$(FC) $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS)  $(CPPFLAGS) $(CPPDEFS) $<
else

.f.o:
	$(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR) $< > $*.F 
	$(FC) $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR) $*.F  

.f90.o:
	$(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR) $< > $*.F90
	$(FC) $(FFLAGS) $(FREEFLAGS) $(INCS) $(INCLDIR) $*.F90  

.F.o:
	$(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR) $< > $*.f 
	$(FC) $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR) $*.f  

.F90.o:
	$(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR) $< > $*.f90
	$(FC) $(FFLAGS) $(FREEFLAGS) $(INCS) $(INCLDIR) $*.f90  
endif

mostlyclean:
	$(RM) -f *.f *.f90 

clean:
	$(RM) -f *.f *.f90 *.d *.a *.$(MOD_SUFFIX) $(OBJS) $(DEPGEN)

realclean:
	$(RM) -f *.f *.f90 *.d *.a *.$(MOD_SUFFIX) $(OBJS) $(DEPGEN)  $(EXEC)

#-------------------------------------------------------------------------------
# Build & include dependency files
#-------------------------------------------------------------------------------
# ASSUMPTIONS:
# o an externally provided dependency generator, $(DEPGEN), is available, 
#   its cmd line syntax is compatible with the build rules below.   Eg, for 
#   each .o file, there is a corresponding .d (dependency) file, and both
#   will be dependent on the same src file, eg.    foo.o foo.d : foo.F90
#   Also, the dependancy genorator's capabilities, limitations, and assumptions
#   are understood & accepted.
#-------------------------------------------------------------------------------

#%.d : %.c $(DEPGEN)
#	@ echo "Building dependency for $@"
#	@ $(DEPGEN) -f $(INCS) $< | head -3  > $@
%.d : %.F $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@
%.d : %.F90 $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@
%.d : %.H $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@

# the if-tests prevent DEPS files from being created when they're not needed
ifneq ($(MAKECMDGOALS), db_files)
ifneq ($(MAKECMDGOALS), db_flags)
ifneq ($(MAKECMDGOALS), mostlyclean)
ifneq ($(MAKECMDGOALS), clean)
ifneq ($(MAKECMDGOALS), realclean)
    -include $(DEPS)
endif
endif
endif
endif
endif
