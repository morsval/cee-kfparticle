#!/usr/bin/env python3
"""Exercise the real CLI parser; stub only ROOT event processing."""
import pathlib
import subprocess
import tempfile

root = pathlib.Path(__file__).resolve().parents[1]
lite = root.parent / 'CeeV0Lite'
with tempfile.TemporaryDirectory(prefix='cee_by_cli_') as tmp:
    tmp = pathlib.Path(tmp)
    stub = tmp / 'analyzer_stub.cxx'
    stub.write_text('''#include "CeeV0Analyzer.h"
#include <iostream>
namespace CeeV0Analysis {
struct Analyzer::Implementation {
  AnalysisConfig config;
  explicit Implementation(const AnalysisConfig& c): config(c) {}
};
Analyzer::Analyzer(const AnalysisConfig& c): fImplementation(new Implementation(c)) {}
Analyzer::~Analyzer() { delete fImplementation; }
int Analyzer::Run() {
  std::cout << "ByTesla=" << fImplementation->config.finder.magneticFieldByTesla << "\\n";
  std::cout << "ProductionState=" << fImplementation->config.finder.buildProductionState << "\\n";
  std::cout << "Reference=" << fImplementation->config.finder.usePrimaryVertexDcaReference << "\\n";
  std::cout << "ReferenceLimit=" << fImplementation->config.finder.referenceTransportConfig.maximumAbsPathCm << "\\n";
  return 0;
}
}
''')
    binary = tmp / 'cli'
    command = ['c++', '-std=c++11', '-Wall', '-Wextra', '-Wpedantic',
               '-I'+str(root/'include'), '-I'+str(lite/'include'),
               str(root/'src/main.cxx'), str(root/'src/CeeV0AnalysisConfig.cxx'),
               str(root/'src/CeeV0AnalysisAlgorithms.cxx'), str(stub)]
    command += [str(p) for p in sorted((lite/'src').glob('*.cxx'))]
    subprocess.run(command + ['-o', str(binary)], check=True)

    def run(arguments, expected=0):
        result = subprocess.run([str(binary)] + arguments, text=True, capture_output=True)
        assert result.returncode == expected, (arguments, result.stdout, result.stderr)
        return result

    assert '--by-tesla' in run(['--help']).stdout
    assert '--no-production-state' in run(['--help']).stdout
    base = ['--input', 'unused.root', '--mode', 'simulation']
    assert 'Reference=1' in run(base).stdout
    assert 'Reference=0' in run(base+['--magnetic-track-reference','stored']).stdout
    assert 'Reference=1' in run(base+['--magnetic-track-reference','pv-dca']).stdout
    assert 'ReferenceLimit=200' in run(base+['--reference-max-path-cm','200']).stdout
    run(base+['--magnetic-track-reference','bad'],1)
    run(base+['--reference-max-path-cm','0'],1)
    run(base+['--reference-max-path-cm','-1'],1)
    assert 'ProductionState=1' in run(base).stdout
    assert 'ProductionState=0' in run(base+['--no-production-state']).stdout
    for field in ('0.5', '-0.5'):
        result = run(base + ['--momentum', 'tracking', '--by-tesla', field])
        assert 'ByTesla='+field in result.stdout
    assert 'ByTesla=0' in run(base + ['--momentum', 'tof', '--by-tesla', '0']).stdout
    result = run(base + ['--momentum', 'tof', '--bz-tesla', '0'])
    assert 'Deprecated' in result.stderr and 'ByTesla=0' in result.stdout
    for field in ('0.5', '-0.5'):
        result = run(base + ['--momentum', 'tracking', '--bz-tesla', field], 1)
        assert 'CEE uses B=(0,By,0)' in result.stderr
    run(base + ['--momentum', 'tof', '--by-tesla', '0.5'], 1)
    run(base + ['--momentum', 'tracking', '--by-tesla', '0'], 1)
    run(base + ['--momentum', 'tracking', '--by-tesla', 'nan'], 1)
print('CEE signed By CLI and legacy Bz rejection tests passed')
