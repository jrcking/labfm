FC := gfortran
LD := gfortran
FFLAGS :=-fopenmp -fbounds-check -ffpe-trap=zero -O3 -Wall -g -J./obj -I./obj -m64
LDFLAGS := -fopenmp -lstdc++ -m64 -lopenblas

## libopenblas is required for this code.

ifeq ($(ob),1)
FFLAGS += -Dob
endif

#ale
ifeq ($(ale),1)
FFLAGS += -Dale=1
endif

# a cylinder obstacle at the origin
ifeq ($(cyl),1)
FFLAGS += -Dcyl=1
endif

# A line of high-resolution
ifeq ($(line),1)
FFLAGS += -Dline=1
endif

## specify order between 2 and 12
ifeq ($(order),2)
FFLAGS += -Dorder=2
endif
ifeq ($(order),3)
FFLAGS += -Dorder=3
endif
ifeq ($(order),4)
FFLAGS += -Dorder=4
endif
ifeq ($(order),5)
FFLAGS += -Dorder=5
endif
ifeq ($(order),6)
FFLAGS += -Dorder=6
endif
ifeq ($(order),7)
FFLAGS += -Dorder=7
endif
ifeq ($(order),8)
FFLAGS += -Dorder=8
endif
ifeq ($(order),9)
FFLAGS += -Dorder=9
endif
ifeq ($(order),10)
FFLAGS += -Dorder=10
endif
ifeq ($(order),11)
FFLAGS += -Dorder=11
endif
ifeq ($(order),12)
FFLAGS += -Dorder=12
endif

SUB_DIRS := para base
SRC_DIR  := $(addprefix source/,$(SUB_DIRS))

#common_parameter needs to be first as mod file is depended upon by nearly everything.
OBJ_FILES := obj/kind_parameters.o obj/common_parameter.o obj/common_2d.o
OBJ_FILES += obj/sphtools.o obj/analytic_functions.o obj/neighbours.o
OBJ_FILES += obj/svd_lib.o
OBJ_FILES += obj/nodes.o obj/moments.o
OBJ_FILES += obj/derivatives.o
OBJ_FILES += obj/basic_convergence_studies.o  
OBJ_FILES += obj/filtering.o
OBJ_FILES += obj/burgers_equation.o obj/td_swarm.o
OBJ_FILES += obj/ns_equations.o obj/ns_equations_cylinder.o obj/ns_equations_tg.o
OBJ_FILES += obj/swarm.o
OBJ_FILES += $(foreach sdir,$(SRC_DIR),$(patsubst $(sdir)/%.F90,obj/%.o,$(wildcard $(sdir)/*.F90)))

HDEPS := $(OBJ_FILES:.o=.d)

vpath %.F90 $(SRC_DIR)

#-------
default: labfm
labfm: $(OBJ_FILES)
	$(LD) -o $@ $^ $(LDFLAGS)

obj/%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -vf ./obj/*.o
	rm -vf ./obj/*.mod
	rm -vf ./labfm
	rm -vf ./fort.*
	rm -vf ./data_out/uv/uv*
	rm -vf ./paraview_files/PART*
	rm -vf ./plots/*.png
	rm -vf ./plots/*.eps

