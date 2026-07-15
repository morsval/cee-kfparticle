# Compiler: g++
CXX = $(shell root-config --cxx)

# Root Compiler flags
ROOTFLAGS = $(shell root-config --cflags)

VMCWORKDIR ?= $(error VMCWORKDIR is not set!)
FAIRROOTPATH ?= $(error FAIRROOTPATH is not set!)

CEEROOT_INCLUDES = -I${VMCWORKDIR} \
                   -I${VMCWORKDIR}/CeeData \
                   -I${VMCWORKDIR}/CEETPC_Tracker \
                   -I${VMCWORKDIR}/t0 \
                   -I${VMCWORKDIR}/t0/edm \
                   -I${VMCWORKDIR}/tpc \
                   -I${VMCWORKDIR}/tpc/edm \
                   -I${VMCWORKDIR}/itof \
                   -I${VMCWORKDIR}/itof/edm \
                   -I${VMCWORKDIR}/etof \
                   -I${VMCWORKDIR}/etof/edm \
                   -I${VMCWORKDIR}/zdc \
                   -I${VMCWORKDIR}/zdc/edm \
                   -I${VMCWORKDIR}/mwdc \
                   -I${VMCWORKDIR}/mwdc/edm 

#CFLAGS = ${ROOTFLAGS} -g -I${FAIRROOTPATH}/include ${CEEROOT_INCLUDES}
CFLAGS = ${ROOTFLAGS} -g -I/software/apps/FairSoft_install_jun19p2/include -I${FAIRROOTPATH}/include ${CEEROOT_INCLUDES}
# Root Library
ROOTLIBS = $(shell root-config --libs)

FAIRLIBDIR = ${FAIRROOTPATH}/lib
CEELIBDIR = ${VMCWORKDIR}/build/lib

LDFLAGS = ${ROOTLIBS} \
          -L${FAIRLIBDIR} -lBase \
					-L${CEELIBDIR} -lCeeData -lCeeT0 -lCeeTpc -lCeeeTOF -lCeeiTOF -lCeeZDC -lCeeMWDC -lCEETPC_Tracker

OBJDIR = obj
TARGET = Analyser.bin

OBJS = $(OBJDIR)/Analyser.o \
       $(OBJDIR)/AnaMaker.o \
       $(OBJDIR)/MyTrack.o \
       $(OBJDIR)/MyLambda.o

all: $(TARGET)

$(shell mkdir -p $(OBJDIR))

$(TARGET): $(OBJS)
	$(CXX) $(CFLAGS) $^ -o $@ $(LDFLAGS)
	strip --strip-all $@

$(OBJDIR)/Analyser.o: analysis.cxx
	$(CXX) -c $(CFLAGS) $< -o $@

$(OBJDIR)/AnaMaker.o: AnaMaker/AnaMaker.cxx AnaMaker/AnaMaker.h
	$(CXX) -c $(CFLAGS) $< -o $@

$(OBJDIR)/MyTrack.o: AnaMaker/MyTrack.cxx AnaMaker/MyTrack.h
	$(CXX) -c $(CFLAGS) $< -o $@

$(OBJDIR)/MyLambda.o: AnaMaker/MyLambda.cxx AnaMaker/MyLambda.h
	$(CXX) -c $(CFLAGS) $< -o $@

# Object files
clean:
	rm -vf $(OBJECTS) $(EXECUTABLE)
	rm -rf $(OBJDIR)
	rm -f $(TARGET)

.PHONY: clean all
