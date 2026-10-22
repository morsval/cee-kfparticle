#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
g++ -O3 -std=c++11 -Wall -Wextra -fPIC -shared -fopenmp -I../../include \
  core_bridge.cxx ../../src/CeeV0Lite{Math,Types,Transport,VertexFitter,ProductionVertex,Finder}.cxx \
  -o libtoy_core.so
