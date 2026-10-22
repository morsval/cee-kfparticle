CXX ?= g++
AR ?= ar

CPPFLAGS += -Iinclude
CXXFLAGS += -O2 -g -std=c++11 -Wall -Wextra -Wpedantic

OBJDIR := obj
LIB := libCeeV0Lite.a

SOURCES := src/CeeV0LiteProductionVertex.cxx \
           src/CeeV0LiteMath.cxx \
           src/CeeV0LiteTypes.cxx \
           src/CeeV0LiteTransport.cxx \
           src/CeeV0LiteVertexFitter.cxx \
           src/CeeV0LiteFinder.cxx

HEADERS := $(wildcard include/*.h)

OBJECTS := $(patsubst src/%.cxx,$(OBJDIR)/%.o,$(SOURCES))
DEPFILES := $(OBJECTS:.o=.d)

.PHONY: all clean test

all: $(LIB)

$(LIB): $(OBJECTS)
	$(AR) rcs $@ $^

$(OBJDIR)/%.o: src/%.cxx
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

test: $(OBJDIR)/test_math_transport $(OBJDIR)/test_end_to_end \
      $(OBJDIR)/test_cee_adapter $(OBJDIR)/test_mother_covariance \
      $(OBJDIR)/test_field_by $(OBJDIR)/test_production_state $(OBJDIR)/test_charged_production \
      $(OBJDIR)/test_transport_jacobian $(OBJDIR)/test_pv_dca_reference
	./$(OBJDIR)/test_math_transport
	./$(OBJDIR)/test_end_to_end
	./$(OBJDIR)/test_cee_adapter
	./$(OBJDIR)/test_mother_covariance
	./$(OBJDIR)/test_field_by
	./$(OBJDIR)/test_production_state
	./$(OBJDIR)/test_charged_production
	./$(OBJDIR)/test_transport_jacobian
	./$(OBJDIR)/test_pv_dca_reference

$(OBJDIR)/test_transport_jacobian: tests/test_transport_jacobian.cxx $(SOURCES) $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter %.cxx,$^) -o $@

$(OBJDIR)/test_charged_production: tests/test_charged_production.cxx $(SOURCES) $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter %.cxx,$^) -o $@

$(OBJDIR)/test_production_state: tests/test_production_state.cxx $(SOURCES) $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter %.cxx,$^) -o $@

.PHONY: reference-test
reference-test: $(OBJDIR)/test_production_state $(OBJDIR)/test_charged_production
	python3 tests/verify_production_reference.py
	python3 tests/verify_charged_reference.py

$(OBJDIR)/test_field_by: tests/test_field_by.cxx $(SOURCES) $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter %.cxx,$^) -o $@

$(OBJDIR)/test_math_transport: tests/test_math_transport.cxx $(SOURCES) $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter %.cxx,$^) -o $@

$(OBJDIR)/test_end_to_end: tests/test_end_to_end.cxx $(SOURCES) $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter %.cxx,$^) -o $@

$(OBJDIR)/test_cee_adapter: tests/test_cee_adapter.cxx \
                           $(SOURCES) $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter %.cxx,$^) -o $@

$(OBJDIR)/test_mother_covariance: tests/test_mother_covariance.cxx $(SOURCES) $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter %.cxx,$^) -o $@

-include $(DEPFILES)

clean:
	rm -rf $(OBJDIR) $(LIB)

$(OBJDIR)/test_pv_dca_reference: tests/test_pv_dca_reference.cxx $(SOURCES) $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(filter %.cxx,$^) -o $@
