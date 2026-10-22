# v0.7.0 validation

- Core: 9 executable test groups passed, including zero/finite field, production states, covariance, transport derivatives and magnetic PV-DCA reference.
- Legacy analysis: algorithms, timing/configuration, C++ syntax and real CLI parsing passed.
- Maker: identical mass and 7D covariance to direct Finder on the same selected input in zero-field TOF and By=0.5 T. PID masks, 2-of-202 subset, post-fit cuts, role aliases, duplicate IDs/indices, clearing between events and error cleanup tested.
- ROOT writer: actual ROOT 6.40 file write/read round trip, selected branches, custom column, int64 event identity, covariance arrays, invalid-state NaN reset and CREATE overwrite protection passed. The expected existing-file diagnostic is part of that test.
- ROOT source: actual ROOT with synthetic CEE-shaped streamable dictionaries; two-file input list, missing branch in the second file, original track IDs, preselection, track-limit error and Finder integration passed.
- Direct macro: Cling loaded a shared library built with the same synthetic dictionaries, ran runCeeKFParticle.C on the magnetic fixture and wrote one mother with DcaProton. This is a macro/API integration test, not an experiment-dictionary link test.
- Manual: 129 topics on 129 intended pages; new chapter, variable tables, cover, contents and version page rendered and visually checked.

Real CEE/FairRoot dictionaries and real detector data are not available here. The supplied synthetic dictionaries must never be used as production libraries. Build `make shared` in the experiment environment before running on actual files.

Core physics source files are byte-identical to the delivered v0.6.5 baseline. The GitHub PR also imports those previously delivered modules, since remote main did not contain them.
