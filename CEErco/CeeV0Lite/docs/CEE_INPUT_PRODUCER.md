# CEE producer requirements

The analysis adapter can consume the new `CeeTpcTrack` and `CeeTpcVertex`
fields, but an accessor alone does not populate a ROOT file. In the current
repository no tracking task calls `CeeTpcTrack::SetVertexSeed`, and no uploaded
vertex producer calls `AddXErr/AddYErr/AddZErr`.

## Present zero-field files

Zero-field TOF reconstruction remains runnable. If a complete valid seed is
absent, the adapter uses `GetX/Y/Z` and `GetPx/Py/Pz` together with the
configured fallback covariance. The analysis output records whether the seed
or the legacy fallback was used, so this is visible rather than silent.

## Required before magnetic-field production

The intended producer reference is the measured track DCA point to the reconstructed primary vertex (not the beam-axis DCA). `SetVertexSeed` only stores its arguments. From v0.6.5, the formal Finder also reanchors valid magnetic inputs to the current event PV DCA in memory, including full state/covariance transport. This accepts older valid seeds at other references without editing the ROOT file. It does not constrain daughter tracks to pass through the PV.

At the same reference surface used for the stored Cartesian state, the track
producer must save

```text
(x,y,z,px,py,pz) and its symmetric 6x6 covariance
```

in units cm and GeV/c, then call

```cpp
track.SetVertexSeed(positionCm, momentumGeV, covariance6x6);
```

For GenFit this should be obtained atomically from the fitted state on one
`MeasuredStateOnPlane` (for example through its 6D-state/covariance conversion),
not by combining a position from one extrapolation with a covariance from
another. The fitted charge must also be stored. Magnetic-field mode rejects a
track when this state, covariance, or charge is missing.

The primary-vertex producer should store the diagonal errors corresponding to
each vertex entry:

```cpp
vertex.AddXErr(sigmaXcm);
vertex.AddYErr(sigmaYcm);
vertex.AddZErr(sigmaZcm);
```

and a rank with the same vector index. Full PV correlations are not present in
the current `CeeTpcVertex` schema; the adapter therefore builds a diagonal PV
covariance.

After changing persistent classes, rebuild the CEE library and ROOT dictionary
before compiling `CeeV0Analysis`. Old files can be read for zero-field studies,
but their new seed/error fields will normally be absent and use the documented
fallback. Magnetic-field validation requires newly produced files.
