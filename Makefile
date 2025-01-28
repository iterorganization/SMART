# SET COMPILER OPTIONS
ifeq ($(FC), ifx)
# INTEL - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  FCFLAGS=-fPIC -fpp -extend-source # FPIC AND PREPROCESSING OPTIONS
  #F90FLAGS+=-g -debug -O0 -fpe-all=0 -no-ftz -traceback -check bounds
else ifeq ($(FC), ifort)
# INTEL - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  FCFLAGS=-fPIC -fpp -extend-source# FPIC AND PREPROCESSING OPTIONS
  #F90FLAGS+=-g -O0 # For debug
  FCFLAGS+=-check all -warn all -gen_interfaces -fpe0 -ftrapuv -traceback -g # For debug
else ifeq ($(FC), gfortran)
# GFORTRAN - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  FCFLAGS=-fPIC -cpp -ffixed-line-length-none # FPIC AND PREPROCESSING OPTIONS
  FCFLAGS+=-Wall -g -fcheck=bounds -O0 -ffpe-trap=invalid,zero,overflow -Wuninitialized
else
  $(error Unsupported Fortran compiler $(FC); exit 1)
endif
FCINC=-I. `pkg-config al-fortran --cflags`
FCLIB=`pkg-config al-fortran --libs`

# XMLLIB LIBRARIES
FCLIB+=`pkg-config xmllib --libs`
FCINC+=`pkg-config xmllib --cflags`

# Constants
FCINC+=`pkg-config fundamental-constants --cflags`

# LIST OF FORTRAN FILES
OBJS_ACTOR = codeparam_smart.o smart.o pelIMAS1.o SMTH.o stepup.o runnt.o
OBJS_STDA  = codeparam_standalone.o standalone.o

all: lib exe

lib: libsmart.a

exe: smart

# STANDALONE COMPILATION
smart: libsmart.a ${OBJS_STDA}
	$(FC) ${FCFLAGS} -o smart ${OBJS_STDA} $(FCLIB) -L. -lsmart

# LIBRARY COMPILATION
libsmart.a: ${OBJS_ACTOR}
	ar -cr $@ $^

# TEST THE CONSISTENCY OF THE XML AND XSD FILES FOR INPUT PARAMETERS
validate:
	xmllint --noout input/standalone.xsd input/standalone.xml
	xmllint --noout --schema input/standalone.xsd input/standalone.xml
	xmllint --noout input/smart.xsd input/smart.xml
	xmllint --noout --schema input/smart.xsd input/smart.xml

# COMPILE THE FC2K ACTOR
actor: libsmart.a 
	sed 's/__COMPILER__/${FC}/g' smart_template.yaml > smart.yaml
	iwrap -f smart.yaml -i $(ACTOR_FOLDER)

# RULES FOR ALL OBJECT FILES
%.o:%.f90
	$(FC) ${FCFLAGS} -c  $< $(FCINC)

%.o:%.f
	$(FC) ${FCFLAGS} -c  $< $(FCINC)

# CLEAN DIRECTORY
clean:
	rm -f smart *.a *.mod *.o smart.xml smart.yaml *.dat
