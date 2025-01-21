include ./compiler.mk
# SET COMPILER OPTIONS
ifeq ($(F90), ifx)
# INTEL - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  F90FLAGS=-fPIC -fpp -extend-source # FPIC AND PREPROCESSING OPTIONS
  #F90FLAGS+=-g -debug -O0 -fpe-all=0 -no-ftz -traceback -check bounds
else ifeq ($(F90), ifort)
# INTEL - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  F90FLAGS=-fPIC -fpp -extend-source# FPIC AND PREPROCESSING OPTIONS
  #F90FLAGS+=-g -O0
else ifeq ($(F90), gfortran)
# GFORTRAN - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  F90FLAGS=-fPIC -cpp -ffixed-line-length-none # FPIC AND PREPROCESSING OPTIONS
  F90FLAGS+=-Wall -g -fcheck=bounds -O0 -ffpe-trap=invalid,zero,overflow -Wuninitialized
else
  $(error Unsupported Fortran compiler $(F90); exit 1)
endif
F90INC=-I. `pkg-config al-fortran --cflags`
F90LIB=`pkg-config al-fortran --libs`

# XMLLIB LIBRARIES
F90LIB+=`pkg-config xmllib --libs`
F90INC+=`pkg-config xmllib --cflags`

# Constants
F90INC+=`pkg-config fundamental-constants --cflags`

# LIST OF FORTRAN FILES
OBJS_ACTOR = codeparam_smart.o smart.o pelIMAS1.o SMTH.o stepup.o runnt.o
OBJS_STDA  = codeparam_standalone.o standalone.o

all: lib exe

lib: libsmart.a

exe: smart

# STANDALONE COMPILATION
smart: libsmart.a ${OBJS_STDA}
	$(F90) ${F90FLAGS} -o smart ${OBJS_STDA} $(F90LIB) -L. -lsmart

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
	sed 's/__COMPILER__/${F90}/g' smart_template.yaml > smart.yaml
	iwrap -f smart.yaml -i $(ACTOR_FOLDER)

# RULES FOR ALL OBJECT FILES
%.o:%.f90
	$(F90) ${F90FLAGS} -c  $< $(F90INC)

%.o:%.f
	$(F90) ${F90FLAGS} -c  $< $(F90INC)

# CLEAN DIRECTORY
clean:
	rm -f smart *.a *.mod *.o smart.xml smart.yaml *.dat
