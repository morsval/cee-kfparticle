TARGET := CeeV0Analysis.bin
OBJDIR := obj
LITEDIR := ../CeeV0Lite

# Default make, `make all`, and combined goals such as `make clean all` need
# the experiment environment. Pure clean/unit-test goals remain standalone.
NEEDS_ROOT := $(if $(MAKECMDGOALS),$(filter all shared libCeeKFParticleAnalysis.so $(TARGET),$(MAKECMDGOALS)),all)
ifneq ($(NEEDS_ROOT),)
CXX := $(shell root-config --cxx)
VMCWORKDIR ?= $(error VMCWORKDIR is not set)
FAIRROOTPATH ?= $(error FAIRROOTPATH is not set)
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)

CEE_INCLUDES := -I$(VMCWORKDIR) \
  -I$(VMCWORKDIR)/CeeData \
  -I$(VMCWORKDIR)/CEETPC_Tracker \
  -I$(VMCWORKDIR)/t0 -I$(VMCWORKDIR)/t0/edm \
  -I$(VMCWORKDIR)/tpc -I$(VMCWORKDIR)/tpc/edm \
  -I$(VMCWORKDIR)/itof -I$(VMCWORKDIR)/itof/edm \
  -I$(VMCWORKDIR)/etof -I$(VMCWORKDIR)/etof/edm

EXTRA_INCLUDES ?=
EXTRA_LDLIBS ?=
CPPFLAGS := -Iinclude -I$(LITEDIR)/include -I$(FAIRROOTPATH)/include \
  $(CEE_INCLUDES) $(EXTRA_INCLUDES)
CXXFLAGS := $(ROOT_CFLAGS) $(if $(filter -std=%,$(ROOT_CFLAGS)),,-std=c++11) -O2 -g -Wall -Wextra -fPIC -MMD -MP
LDLIBS := $(ROOT_LIBS) -L$(FAIRROOTPATH)/lib -lBase \
  -L$(VMCWORKDIR)/build/lib -lCeeData -lCeeT0 -lCeeTpc \
  -lCeeeTOF -lCeeiTOF -lCEETPC_Tracker $(EXTRA_LDLIBS)
else
CXX ?= g++
endif

ANALYSIS_SOURCES := src/CeeV0AnalysisAlgorithms.cxx \
  src/CeeV0AnalysisConfig.cxx src/CeeV0Analyzer.cxx src/main.cxx \
  src/CeeKFParticleInterface.cxx src/CeeKFParticleAnalysisMaker.cxx \
  src/CeeKFParticleObservables.cxx src/CeeKFParticleRootSource.cxx src/CeeKFParticleRootWriter.cxx
LITE_SOURCES := $(LITEDIR)/src/CeeV0LiteProductionVertex.cxx \
  $(LITEDIR)/src/CeeV0LiteMath.cxx \
  $(LITEDIR)/src/CeeV0LiteTypes.cxx \
  $(LITEDIR)/src/CeeV0LiteTransport.cxx \
  $(LITEDIR)/src/CeeV0LiteVertexFitter.cxx \
  $(LITEDIR)/src/CeeV0LiteFinder.cxx
SOURCES := $(ANALYSIS_SOURCES) $(LITE_SOURCES)
OBJECTS := $(patsubst src/%.cxx,$(OBJDIR)/analysis_%.o,$(ANALYSIS_SOURCES)) \
  $(patsubst $(LITEDIR)/src/%.cxx,$(OBJDIR)/lite_%.o,$(LITE_SOURCES))
DEPFILES := $(OBJECTS:.o=.d)

.PHONY: all clean unit-test config-test syntax-test field-cli-test check
all: $(TARGET) libCeeKFParticleAnalysis.so

shared: libCeeKFParticleAnalysis.so

libCeeKFParticleAnalysis.so: $(filter-out $(OBJDIR)/analysis_main.o,$(OBJECTS))
	$(CXX) -shared $^ -o $@ $(LDLIBS)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(OBJDIR)/analysis_%.o: src/%.cxx | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/lite_%.o: $(LITEDIR)/src/%.cxx | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $@

# This target intentionally needs no ROOT, FairRoot or CEE installation.
unit-test:
	$(CXX) -std=c++11 -Wall -Wextra -pedantic -Iinclude \
	  src/CeeV0AnalysisAlgorithms.cxx tests/test_algorithms.cxx \
	  -o /tmp/test_cee_v0_analysis_algorithms
	/tmp/test_cee_v0_analysis_algorithms

config-test:
	$(CXX) -std=c++11 -Wall -Wextra -pedantic \
	  -Iinclude -I$(LITEDIR)/include \
	  src/CeeV0AnalysisAlgorithms.cxx src/CeeV0AnalysisConfig.cxx \
	  tests/test_config.cxx $(LITE_SOURCES) \
	  -o /tmp/test_cee_v0_analysis_config
	/tmp/test_cee_v0_analysis_config

# Parse the complete analysis source with API-shaped stubs. This catches C++
# regressions on machines without ROOT; it is not a substitute for a real CEE
# link and ROOT-I/O test.
syntax-test:
	$(CXX) -std=c++11 -Wall -Wextra -Wpedantic -fsyntax-only \
	  -Itests/stubs -Iinclude -I$(LITEDIR)/include \
	  src/CeeV0AnalysisConfig.cxx src/CeeV0Analyzer.cxx src/main.cxx \
	  src/CeeKFParticleRootSource.cxx src/CeeKFParticleRootWriter.cxx

field-cli-test:
	python3 tests/test_field_cli.py

check: unit-test config-test syntax-test field-cli-test maker-test

clean:
	rm -rf $(OBJDIR) $(TARGET) libCeeKFParticleAnalysis.so

-include $(DEPFILES)

.PHONY: maker-test root-writer-test shared
MAKER_SOURCES := src/CeeKFParticleInterface.cxx src/CeeKFParticleAnalysisMaker.cxx src/CeeKFParticleObservables.cxx
maker-test:
	$(CXX) -std=c++11 -O2 -Wall -Wextra -Wpedantic -Iinclude -I$(LITEDIR)/include \
	  $(MAKER_SOURCES) $(LITE_SOURCES) tests/test_maker.cxx -o /tmp/test_cee_kfp_maker
	/tmp/test_cee_kfp_maker
# Requires ROOT only; no CEE/FairRoot and no compile-shape stubs.
root-writer-test:
	$(CXX) $$(root-config --cflags) -std=c++17 -O2 -Iinclude -I$(LITEDIR)/include \
	  $(MAKER_SOURCES) src/CeeKFParticleRootWriter.cxx $(LITE_SOURCES) \
	  tests/test_root_writer.cxx $$(root-config --libs) -o /tmp/test_cee_kfp_writer
	/tmp/test_cee_kfp_writer

# Real ROOT I/O with synthetic CEE-shaped dictionaries, NOT the experiment EDM.
.PHONY: root-source-test
root-source-test:
	mkdir -p /tmp/cee_kfp_source_test
	rootcling -f /tmp/cee_kfp_source_test/dictionary.cxx -Itests/root_fixtures \
	  tests/root_fixtures/CeeRootTestTypes.h tests/root_fixtures/LinkDef.h
	$(CXX) $$(root-config --cflags) -std=c++17 -O2 -I. -Itests/root_fixtures -Iinclude -I$(LITEDIR)/include \
	  $(MAKER_SOURCES) src/CeeKFParticleRootSource.cxx src/CeeV0AnalysisAlgorithms.cxx \
	  src/CeeV0AnalysisConfig.cxx $(LITE_SOURCES) /tmp/cee_kfp_source_test/dictionary.cxx \
	  tests/test_root_source.cxx $$(root-config --libs) -o /tmp/cee_kfp_source_test/test
	/tmp/cee_kfp_source_test/test
