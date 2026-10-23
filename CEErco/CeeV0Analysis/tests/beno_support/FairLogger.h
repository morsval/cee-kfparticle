// Test-only logging shim. Physics EDM and ROOT I/O use actual CeeRoot sources.
#ifndef CEE_BENO_TEST_FAIR_LOGGER_H
#define CEE_BENO_TEST_FAIR_LOGGER_H
#include <iostream>
#define LOG(level) std::cerr
#endif
