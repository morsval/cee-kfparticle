#ifndef CEE_PICO_DST_ADAPTER_EXAMPLE_H
#define CEE_PICO_DST_ADAPTER_EXAMPLE_H
#include "CeeKFParticleEvent.h"
#include <functional>
// A compileable adapter for an external reader. This does not claim a CEE
// picoDst schema already exists. The callback owns that format-specific mapping.
class PicoDstAdapterExample : public CeeV0Analysis::CeeKFParticleEventSource {
 public:
  typedef std::function<ReadStatus(CeeV0Analysis::AnalysisEvent&,std::string&)> ReadEvent;
  explicit PicoDstAdapterExample(ReadEvent read) : fRead(read) {}
  ReadStatus Next(CeeV0Analysis::AnalysisEvent& e,std::string& error) override {
    e=CeeV0Analysis::AnalysisEvent(); error.clear(); return fRead(e,error);
  }
 private:
  ReadEvent fRead;
};
// In the callback, fill eventId, input.primaryVertex (position + covariance),
// timing, and ONLY selected AnalysisTrack objects. Every magnetic TrackInput
// needs an ID, measured charge, physical momentum, full 6x6 Cartesian covariance
// and position at the SAME reference. Fill calibrated PID if available.
// Make no extra daughter-PV constraint. Maker/Finder applies its usual default
// track-PV DCA reference preparation after PID, once per selected track.
#endif
