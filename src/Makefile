#### Settings..

compiler       = gfortran # ifort or gfortran
compiler_flags = debug    # normal, debug, or timing
linker_style   = dynamic  # dynamic or static
include_netcdf = false    # true or false

###################################################

#### Declaration of source files..

general_SOURCES = allocate_carbon.f90   canopy.f90       config_tools.f90     global_variables.f90  \
                  leaf.f90              light.f90        linked_lists.f90     log_tools.f90         \
                  main.f90              math_tools.f90   soil_air.f90         soil_functions.f90    \
                  spa_cmd_line.f90      spa_config.f90   spa_initialise.f90   spa_io.f90            \
                  spa_io_csv.f90

netcdf_SOURCES = netcdf_tools.f90  spa_io_netcdf.f90  spa_restart.f90

###################################################

#### Declaration of general compiler/linker flags..
# ( which one actually gets used depends upon what )
# (  the user specifies at the top of this file.   )

timing = -pg

# For the Intel fortran compiler..(v11 or newer!)
ifort_normal  = -fpp -vec-report=0
ifort_debug   = $(ifort_normal) -check all -debug-parameters all -ftrapuv -g -fpe0 \
                -implicitnone -O0 -p -stand f03 -traceback 
ifort_timing  = $(ifort_normal) $(timing) 
ifort_dynamic = -i-dynamic                   # use the dynamic intel libs
ifort_static  = -static -i-static -I/usr/lib # only use the static intel libraries

# For the GNU fortran compiler..(v4.4.5 or newer!)
gfortran_normal  = -cpp -ffree-line-length-none 
gfortran_debug   = $(gfortran_normal) -fbounds-check -fimplicit-none -frange-check \
                   -ftree-vectorizer-verbose=0 -ggdb -O0 -pedantic -std=f2003 -Wall
gfortran_timing  = $(gfortran_normal) $(timing)
gfortran_dynamic =                           #
gfortran_static  = -static # -static-libgcc  # only use the static gfortran libraries

#############################################
# PLEASE DO NOT TOUCH BELOW HERE UNLESS YOU #
# ARE HAPPY YOU KNOW HOW MAKEFILES WORK.    #
#############################################

#### Trim any white-space off the user's entries..

FC=$(strip $(compiler))
COMP=$(strip $(compiler_flags))
LINK=$(strip $(linker_style))
INC_NC=$(strip $(include_netcdf))

#### If user selected gfortran, make sure it's a new enough version for us to use...

#
# THERE IS _NO_ WAY TO EXIT A MAKEFILE, AS IT IS NOT A SCRIPT, RUN IN LINEAR ORDER!
#
#ifeq ($(FC),gfortran)
#	FC_VERSION=$(shell gfortran --version|grep -i "GNU Fortran (GCC)"| cut -d " " -f4)
#endif
#ifeq ($(FC_VERSION),4.1.2)
#	exit
#endif

#### Select the compiler/linker flags the user asked for...

ifeq ($(INC_NC),true)	# combine all linker flags
	NC_PREPROCESSOR_FLAG = -D USE_NETCDF
	NC_INC_DIR = $(shell nc-config --prefix)/include
	NC_LIB_DIR = $(shell nc-config --prefix)/lib
	NC_LINK_FLAGS = -L$(NC_LIB_DIR) -lnetcdf -lnetcdff
endif
ifeq ($(FC),gfortran)
	NC_INC_FLAGS = -I$(ncpath)/include -I/usr/lib64/gfortran/modules/
else
	NC_INC_FLAGS = -I$(ncpath)/include
endif

FC_FLAGS   = $($(FC)_$(COMP)) $(NC_PREPROCESSOR_FLAG)
INC_FLAGS  = $(NC_INC_FLAGS)
LINK_FLAGS = $($(FC)_$(LINK)) $(NC_LINK_FLAGS)

#### Select the sources files the user asked for..

ifeq ($(INC_NC),true)
	spa_SOURCES = $(general_SOURCES) $(netcdf_SOURCES)
else
	spa_SOURCES = $(general_SOURCES)
endif

#### Instructions for building SPA...

spa:	$(spa_SOURCES:%.f90=%.o)
	$(FC) -o $@ $^ $($(COMP)) $(LINK_FLAGS)
	mv spa ../

clean:
	rm -f *.mod *.o make_dep ../spa ../gmon.out

%.o:	%.f90
	$(FC) $(FC_FLAGS) $(INC_FLAGS) -c $< -o $@

Makefile:	make_dep

make_dep:	$(spa_SOURCES)
		./f90_dependencies.py $^ > $@

-include make_dep
