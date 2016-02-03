#=============================================================================
# Makefile for vertical transport model
#=============================================================================
SHELL := /bin/bash

include Makefile_header.mk
#=============================================================================
#  Compile Flags
#=============================================================================
# IFORT compilation options
#FFLAGS = -cpp -w -O3 -auto -noalign -convert big_endian 

# IFORT compilation options (for debugging: no optimization)
#FFLAGS = -cpp -w -O0 -auto -noalign -convert big_endian -g -traceback -CB -vec-report0

# Add options for medium memory model.  This is to prevent G-C from 
# running out of memory at hi-res, especially when using netCDF I/O. 
#FFLAGS += -mcmodel=medium -i-dynamic
#FFLAGS += -shared-intel -mcmodel=large -i-dynamic

#=============================================================================
#  Libraries
#=============================================================================
# NetCDF
FFLAGS += -I$(GC_INCLUDE)
LINK_NC := $(shell $(GC_BIN)/nf-config --flibs)
LINK_NC += $(shell $(GC_BIN)/nc-config --libs)
LINK_NC += $(filter -l%,$(NCL))

ifeq ($(LINK_NC),)
LINK_NC   := -lnetcdf -lhdf5_hl -lhdf5 -lm -z
endif

# Prepend the library directory path to the linking sequence
LINK_NC   := -L$(LIB_NETCDF) $(LINK_NC)

#=============================================================================
#  Compile Command
#=============================================================================
# Compile command
#F90  = ifort $(FFLAGS)         # Single processor
#F90 += -openmp -Dmultitask -mp # Multiprocessor


#=============================================================================
#  Objects and Modules
#=============================================================================
MODS =                        \
ifort_errmsg.o		      \
error_mod.o		      \
inquireMod.o		      \
gigc_errcode_mod.o	      \
charpak_mod.o		      \
CMN_SIZE_mod.o		      \
CMN_GCTM_mod.o		      \
julday_mod.o 		      \
m_do_err_out.o		      \
m_netcdf_io_checks.o	      \
m_netcdf_io_close.o	      \
m_netcdf_io_create.o	      \
m_netcdf_io_define.o	      \
m_netcdf_io_get_dimlen.o      \
m_netcdf_io_handle_err.o      \
m_netcdf_io_open.o	      \
m_netcdf_io_read.o	      \
m_netcdf_io_readattr.o	      \
m_netcdf_io_write.o	      \
ncdf_mod.o		      \
time_mod.o		      \
pressure_mod.o		      \
grid_mod.o		      \
pjc_pfix_mod.o		      \
tpcore_fvdas_mod.o	      \
geosfp_read_mod.o	      \
cleanup.o		      \
initialize_mod.o	      \
advection_mod.o		      \
chem_mod.o				  \
diagnostic_mod.o	      \
main.o

#=============================================================================
#  Executable
#=============================================================================
fzppm:  $(MODS)
	$(F90)  $(MODS) $(LINK_NC) -o fzppm

#=============================================================================
#  Dependencies Listing
#=============================================================================
cleanup.o		      : cleanup.F
CMN_GCTM_mod.o		      : CMN_GCTM_mod.F
CMN_SIZE_mod.o		      : CMN_SIZE_mod.F	  gigc_errcode_mod.o
charpak_mod.o		      : charpak_mod.F
diagnostic_mod.o	      : diagnostic_mod.F90			\
				julday_mod.o
error_mod.o		      : error_mod.F
geosfp_read_mod.o	      : geosfp_read_mod.F90  CMN_SIZE_mod.o     
gigc_errcode_mod.o	      : gigc_errcode_mod.F90
grid_mod.o		      : grid_mod.F90	     gigc_errcode_mod.o	\
				CMN_GCTM_mod.o	     error_mod.o
ifort_errmsg.o		      : ifort_errmsg.F
inquireMod.o		      : inquireMod.F90
initialize_mod.o	      : initialize_mod.F90   grid_mod.o		\
				CMN_SIZE_mod.o	     CMN_GCTM_mod.o	\
				error_mod.o
julday_mod.o		      : julday_mod.F
m_do_err_out.o		      : m_do_err_out.F90
m_netcdf_io_checks.o	      : m_netcdf_io_checks.F90
m_netcdf_io_close.o	      : m_netcdf_io_close.F90 	m_do_err_out.o
m_netcdf_io_create.o	      : m_netcdf_io_create.F90	m_do_err_out.o
m_netcdf_io_define.o	      : m_netcdf_io_define.F90
m_netcdf_io_io_get_dimlen.o   : m_netcdf_io_get_dimlen.F90		\
				m_do_err_out.o
m_netcdf_io_handle_err.o      : m_netcdf_io_handle_err.F90		\
				m_do_err_out.o
m_netcdf_io_open.o	      : m_netcdf_io_open.F90	m_do_err_out.o
m_netcdf_io_read.o	      : m_netcdf_io_read.F90	m_do_err_out.o
m_netcdf_io_readattr.o	      : m_netcdf_io_readattr.F90		\
				m_do_err_out.o
m_netcdf_io_write.o	      : m_netcdf_io_write.F90	m_do_err_out.o
ncdf_mod.o		      : ncdf_mod.F90	    m_netcdf_io_open.o	\
				m_netcdf_io_get_dimlen.o		\
				m_netcdf_io_read.o  m_netcdf_io_close.o \
				m_netcdf_io_readattr.o			\
				m_netcdf_io_create.o			\
				m_netcdf_io_define.o			\
				m_netcdf_io_write.o m_netcdf_io_checks.o 	
main.o                        : main.F90          CMN_SIZE_mod.o 	\
				CMN_GCTM_mod.o 	  error_mod.o		\
				geosfp_read_mod.o time_mod.o 		\
				tpcore_fvdas_mod.o			\
				pjc_pfix_mod.o	  grid_mod.o		\
				initialize_mod.o  advection_mod.o	\
				diagnostic_mod.o 
pjc_pfix_mod.o		      : pjc_pfix_mod.F	  CMN_SIZE_mod.o	\
				CMN_GCTM_mod.o	  error_mod.o		\
				pressure_mod.o	  grid_mod.o
pressure_mod.o		      : pressure_mod.F	  CMN_SIZE_mod.o   	\
				error_mod.o
time_mod.o		      : time_mod.F	  julday_mod.o
tpcore_fvdas_mod.o	      : tpcore_fvdas_mod.F90
advection_mod.o		      : advection_mod.F90 	pressure_mod.o	\
				CMN_SIZE_mod.o		CMN_GCTM_mod.o	\
				error_mod.o		gigc_errcode_mod.o \
				grid_mod.o		tpcore_fvdas_mod.o \
				pjc_pfix_mod.o
chem_mod.o				  : chem_mod.F90			grid_mod.o	\
				CMN_SIZE_mod.o		error_mod.o

#=============================================================================
#  Other Makefile Commands
#=============================================================================
clean:
	rm -rf *.o *.mod ifc* fzppm rii_files

.PHONY: clean

.SUFFIXES: .f .F .f90 .F90
.f.o:			      ; $(F90) -c $*.f
.F.o:			      ; $(F90) -c $*.F
.f90.o:                       ; $(F90) -c -free $*.f90 
.F90.o:                       ; $(F90) -c -free $*.F90 

#=============================================================================
#  End
#=============================================================================
