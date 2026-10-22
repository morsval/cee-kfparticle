# CEE build entry point. Load the CeeRoot environment before `make`.
.DEFAULT_GOAL := all
.PHONY: all lite analysis legacy check clean help
all: lite analysis
lite:
	+$(MAKE) -C CeeV0Lite all
analysis:
	+$(MAKE) -C CeeV0Analysis all
legacy:
	+$(MAKE) -f Makefile.legacy all
check:
	+$(MAKE) -C CeeV0Lite test
	+$(MAKE) -C CeeV0Analysis check
clean:
	+$(MAKE) -C CeeV0Lite clean
	+$(MAKE) -C CeeV0Analysis clean
help:
	@echo 'make [-jN]: CeeV0Lite + analysis executable + Maker shared library'
	@echo 'make check: standalone regression tests; make legacy: old Analyser.bin'
