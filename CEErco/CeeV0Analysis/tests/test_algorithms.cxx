#include "CeeV0AnalysisAlgorithms.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace
{

bool Close(double first, double second, double tolerance = 1.e-10)
{
  return std::fabs(first - second) < tolerance;
}

} // namespace

int main()
{
  using namespace CeeV0Analysis;

  ProjectionTrack sideTrack;
  sideTrack.id = 7;
  sideTrack.direction[0] = 1.;
  sideTrack.direction[1] = 0.5;
  sideTrack.direction[2] = -0.25;
  TofHitData sideHit;
  sideHit.id = 2;
  sideHit.detector = kItofDetector;
  sideHit.module = 3;
  sideHit.positionCm[0] = 10.;
  sideHit.positionCm[1] = 5.3;
  sideHit.positionCm[2] = -2.1;
  assert(Close(ProjectionDistanceCm(sideTrack, sideHit), 0.5));

  ProjectionTrack forwardTrack;
  forwardTrack.id = 8;
  forwardTrack.direction[0] = 0.1;
  forwardTrack.direction[1] = -0.2;
  forwardTrack.direction[2] = 1.;
  TofHitData forwardHit;
  forwardHit.id = 3;
  forwardHit.detector = kEtofDetector;
  forwardHit.positionCm[0] = 2.2;
  forwardHit.positionCm[1] = -4.1;
  forwardHit.positionCm[2] = 20.;
  assert(Close(ProjectionDistanceCm(forwardTrack, forwardHit),
               std::sqrt(0.05)));

  // Both tracks can use hit 0; the global assignment gives it to the exact
  // track and sends the second track to hit 1.
  ProjectionTrack secondTrack = forwardTrack;
  secondTrack.id = 9;
  secondTrack.positionCm[0] = 0.2;
  TofHitData exactHit = forwardHit;
  exactHit.positionCm[0] = 2.;
  exactHit.positionCm[1] = -4.;
  TofHitData alternateHit = exactHit;
  alternateHit.id = 4;
  alternateHit.positionCm[0] = 2.45;
  std::vector<ProjectionTrack> tracks;
  tracks.push_back(forwardTrack);
  tracks.push_back(secondTrack);
  std::vector<TofHitData> hits;
  hits.push_back(exactHit);
  hits.push_back(alternateHit);
  const std::vector<TofMatch> matches =
    AssignUniqueTofHits(tracks, hits, 1., 1.);
  assert(matches.size() == 2u);
  assert(matches[0].trackId == 8 && matches[0].hitId == 3);
  assert(matches[1].trackId == 9 && matches[1].hitId == 4);

  // Regression: a lowest-edge greedy matcher would choose A-X and leave B
  // unmatched. Maximum-cardinality min-cost matching must choose A-Y + B-X.
  ProjectionTrack trackA;
  trackA.id = 100;
  trackA.direction[2] = 1.;
  ProjectionTrack trackB = trackA;
  trackB.id = 101;
  trackB.positionCm[0] = -0.5;
  TofHitData hitX;
  hitX.id = 200;
  hitX.detector = kEtofDetector;
  hitX.positionCm[0] = 1.;
  hitX.positionCm[2] = 10.;
  TofHitData hitY = hitX;
  hitY.id = 201;
  hitY.positionCm[0] = 2.;
  std::vector<ProjectionTrack> counterTracks;
  counterTracks.push_back(trackA);
  counterTracks.push_back(trackB);
  std::vector<TofHitData> counterHits;
  counterHits.push_back(hitX);
  counterHits.push_back(hitY);
  const std::vector<TofMatch> counterMatches =
    AssignUniqueTofHits(counterTracks, counterHits, 2.1, 2.1);
  assert(counterMatches.size() == 2u);
  assert(counterMatches[0].trackId == 100 && counterMatches[0].hitId == 201);
  assert(counterMatches[1].trackId == 101 && counterMatches[1].hitId == 200);

  T0SelectionConfig t0;
  t0.source = kT0Hit;
  t0.aggregation = kMedianT0;
  t0.rawToNsScale = 1.e-3;
  t0.offsetNs = -0.2;
  std::vector<double> rawTimes;
  rawTimes.push_back(1000.);
  rawTimes.push_back(3000.);
  rawTimes.push_back(2000.);
  const T0SelectionResult selected = SelectEventT0(rawTimes, t0);
  assert(selected.valid);
  assert(Close(selected.selectedRawTime, 2000.));
  assert(Close(selected.timeNs, 1.8));

  t0.aggregation = kMeanT0;
  assert(Close(SelectEventT0(rawTimes, t0).selectedRawTime, 2000.));
  t0.aggregation = kEarliestT0;
  assert(Close(SelectEventT0(rawTimes, t0).selectedRawTime, 1000.));
  t0.aggregation = kMedianT0;
  rawTimes.push_back(4000.);
  assert(Close(SelectEventT0(rawTimes, t0).selectedRawTime, 2500.));

  t0.rawToNsScale = -1.;
  assert(!SelectEventT0(rawTimes, t0).valid);
  t0.rawToNsScale = 1.e-3;
  const std::vector<double> emptyTimes;
  assert(!SelectEventT0(emptyTimes, t0).valid);

  t0.source = kFixedT0;
  t0.fixedTimeNs = 4.5;
  const T0SelectionResult fixed = SelectEventT0(rawTimes, t0);
  assert(fixed.valid && Close(fixed.timeNs, 4.5));

  std::cout << "CeeV0Analysis algorithm tests passed" << std::endl;
  return 0;
}
