# DEFINE F90
COMPILER = ifort
include ./compiler.mk

# SET COMPILER OPTIONS
ifeq ($(F90), ifort)
# INTEL - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  F90FLAGS=-fPIC -fpp -g # FPIC AND PREPROCESSING OPTIONS
else ifeq ($(F90), gfortran)
# GFORTRAN - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  F90FLAGS=-fPIC -cpp # FPIC AND PREPROCESSING OPTIONS
else
  $(warning Unsupported Fortran compiler $(F90), proceed with care...)
  F90FLAGS=-fPIC -cpp # FPIC AND PREPROCESSING OPTIONS
endif
F90INC=-I. `pkg-config al-fortran --cflags`
F90LIB=`pkg-config al-fortran --libs`

# XMLLIB LIBRARIES
F90LIB+=`pkg-config xmllib --libs`
F90INC+=`pkg-config xmllib --cflags`

# LIST OF FORTRAN FILES
OBJS_ACTOR = codeparam_smart.o smart.o pelIMAS.o SMTH.o
OBJS_STDA  = codeparam_standalone.o smart.o standalone.o pelIMAS.o SMTH.o

all: lib actor

lib: libsmart.a

exe: standalone.exe

# STANDALONE COMPILATION
standalone.exe: libsmart.a ${OBJS_STDA}
	$(F90) ${F90FLAGS} -o standalone.exe ${OBJS_STDA} $(F90LIB) -L. -lsmart

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
	rm -f *.exe *.a *.mod *.o smart.xml smart.yaml

