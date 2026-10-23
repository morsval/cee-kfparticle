# CEE input adapter

`CeeV0LiteCeeAdapter.h` is the checked boundary between persisted CEE objects
and the detector-independent V0 finder.  All lengths are cm, momenta are
GeV/c, times are ns, and covariances follow the state order
`(x,y,z,px,py,pz)`.

## Track modes

- `kCeeZeroFieldTofMode`: the `GetP*` or vertex-seed momentum is used only as
  a direction.  The adapter sets `hasTrackingMomentum=false`, marks charge as
  unknown, and requires a valid `CeeTofMeasurement` to determine momentum.
- `kCeeMagneticFieldTrackingMode`: requires `HasVertexSeed()`, a finite seed
  position and momentum, a usable 6x6 seed covariance, and a measured fitted
  charge.  TOF is optional, so this mode also supports tracking-only input.

When a complete valid vertex seed exists, its position, momentum, and covariance
are validated and selected together. This prevents mixing a state at one
reference surface with a covariance from another. In magnetic tracking mode the
complete 6x6 seed covariance is passed to Finder. In zero-field TOF mode only the
3x3 position block and a scalar direction resolution compressed from the
momentum block are used; anisotropic transverse direction errors,
position--direction cross terms, and the old longitudinal momentum error are not
retained when TOF supplies a new momentum magnitude. Zero-field mode may fall
back as one unit to `GetX/Y/Z`, `GetPx/Py/Pz`, and configured resolutions;
magnetic-field mode never does so.

TOF matching must use the public `SelectTrackReferenceState` helper.  It makes
exactly the same seed/covariance validation and fallback decision later used by
`MakeTrackInput`; independently checking only `HasVertexSeed()` and seed
position/momentum can otherwise select a different reference surface.

```cpp
CeeV0Lite::CeeTrackReferenceState reference;
CeeV0Lite::CeeAdapterStatus referenceStatus;
if (!CeeV0Lite::SelectTrackReferenceState(
      track, CeeV0Lite::kCeeZeroFieldTofMode, reference,
      adapterConfig, &referenceStatus)) {
  // reject track
}
// Match hits to the line defined by reference.positionCm/reference.direction.
```

`reference.usedVertexSeed` records the decision directly.  Successful status
codes also distinguish `kCeeAdapterSuccessVertexSeed` from
`kCeeAdapterSuccessLegacyFallback`; `CeeAdapterStatusIsSuccess(status)` accepts
either one (as well as the generic vertex-conversion success status).

```cpp
CeeV0Lite::CeeTofMeasurement tof;
tof.available = true;
tof.positionCm = CeeV0Lite::Vec3(hitX, hitY, hitZ);
tof.timeNs = hitTimeNs; // absolute time; finder subtracts EventInput start time

CeeV0Lite::TrackInput input;
CeeV0Lite::CeeAdapterStatus status;
const bool ok = CeeV0Lite::MakeTrackInput(
  track, trackIndex, CeeV0Lite::kCeeZeroFieldTofMode, tof,
  CeeV0Lite::kPionHypothesis | CeeV0Lite::kProtonHypothesis,
  input, CeeV0Lite::CeeAdapterConfig(), &status);
```

## Primary vertex

`MakePrimaryVertex` checks the coordinate-vector lengths and finite values,
then selects the valid entry with the highest rank.  Its `XErr/YErr/ZErr` are
used when positive and finite.  Missing or invalid errors use configurable
fallback resolutions unless `allowVertexErrorFallback` is disabled.

```cpp
CeeV0Lite::PrimaryVertex primary;
std::size_t selected = 0;
const bool ok = CeeV0Lite::MakePrimaryVertex(
  ceeVertex, primary, &selected, config, &status);
```

Always check the returned `bool`.  `CeeAdapterStatusMessage(status)` provides a
short diagnostic suitable for a rejection counter or log.  The standalone
`test_cee_adapter` target uses mock CEE objects, so this boundary is compiled
and tested without ROOT/FairRoot.
