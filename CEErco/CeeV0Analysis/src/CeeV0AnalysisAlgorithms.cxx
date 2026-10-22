#include "CeeV0AnalysisAlgorithms.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>

namespace CeeV0Analysis
{

namespace
{

const double kInvalidDistanceCm = 1.e30;

bool IsFinite3(const double value[3])
{
  return std::isfinite(value[0]) && std::isfinite(value[1]) &&
         std::isfinite(value[2]);
}

struct MatchEdge
{
  MatchEdge(int destination, int reverseIndex, int capacityValue,
            double costValue, int trackValue, int hitValue,
            double distanceValue)
    : to(destination), reverse(reverseIndex), capacity(capacityValue),
      cost(costValue), trackIndex(trackValue), hitIndex(hitValue),
      distanceCm(distanceValue) {}

  int to;
  int reverse;
  int capacity;
  double cost;
  int trackIndex;
  int hitIndex;
  double distanceCm;
};

void AddFlowEdge(std::vector<std::vector<MatchEdge> >& graph,
                 int from, int to, double cost,
                 int trackIndex = -1, int hitIndex = -1,
                 double distanceCm = 0.)
{
  const int forwardReverse = static_cast<int>(graph[to].size());
  const int backwardReverse = static_cast<int>(graph[from].size());
  graph[from].push_back(MatchEdge(to, forwardReverse, 1, cost,
                                  trackIndex, hitIndex, distanceCm));
  graph[to].push_back(MatchEdge(from, backwardReverse, 0, -cost,
                                -1, -1, distanceCm));
}

struct TrackOrder
{
  bool operator()(const TofMatch& first, const TofMatch& second) const
  {
    return first.trackIndex < second.trackIndex;
  }
};

} // namespace

ProjectionTrack::ProjectionTrack() : id(-1)
{
  for (int i = 0; i < 3; ++i) {
    positionCm[i] = 0.;
    direction[i] = 0.;
  }
}

TofHitData::TofHitData()
  : id(-1), detector(kNoTofDetector), module(-1), rawTime(0.)
{
  for (int i = 0; i < 3; ++i) positionCm[i] = 0.;
}

TofMatch::TofMatch()
  : trackIndex(0), hitIndex(0), trackId(-1), hitId(-1),
    detector(kNoTofDetector), distanceCm(kInvalidDistanceCm)
{
}

double ProjectionDistanceCm(const ProjectionTrack& track,
                            const TofHitData& hit)
{
  if (!IsFinite3(track.positionCm) || !IsFinite3(track.direction) ||
      !IsFinite3(hit.positionCm)) return kInvalidDistanceCm;

  const double x0 = track.positionCm[0];
  const double y0 = track.positionCm[1];
  const double z0 = track.positionCm[2];
  const double dx = track.direction[0];
  const double dy = track.direction[1];
  const double dz = track.direction[2];

  if (hit.detector == kItofDetector &&
      ((hit.module >= 0 && hit.module <= 8) ||
       (hit.module >= 15 && hit.module <= 24))) {
    if (std::fabs(dx) < 1.e-12) return kInvalidDistanceCm;
    const double path = (hit.positionCm[0] - x0) / dx;
    const double residualY = hit.positionCm[1] - (y0 + path * dy);
    const double residualZ = hit.positionCm[2] - (z0 + path * dz);
    return std::sqrt(residualY * residualY + residualZ * residualZ);
  }

  if (hit.detector == kEtofDetector ||
      (hit.detector == kItofDetector && hit.module >= 9 &&
       hit.module <= 14)) {
    if (std::fabs(dz) < 1.e-12) return kInvalidDistanceCm;
    const double path = (hit.positionCm[2] - z0) / dz;
    const double residualX = hit.positionCm[0] - (x0 + path * dx);
    const double residualY = hit.positionCm[1] - (y0 + path * dy);
    return std::sqrt(residualX * residualX + residualY * residualY);
  }

  // Unknown iTOF module: retain the old AnaMaker xy-line fallback.
  if (hit.detector == kItofDetector) {
    const double norm2 = dx * dx + dy * dy;
    if (norm2 < 1.e-24) return kInvalidDistanceCm;
    const double rx = x0 - hit.positionCm[0];
    const double ry = y0 - hit.positionCm[1];
    return std::fabs(rx * dy - ry * dx) / std::sqrt(norm2);
  }

  return kInvalidDistanceCm;
}

std::vector<TofMatch> AssignUniqueTofHits(
    const std::vector<ProjectionTrack>& tracks,
    const std::vector<TofHitData>& hits,
    double maximumItofDistanceCm,
    double maximumEtofDistanceCm)
{
  const int source = 0;
  const int firstTrackNode = 1;
  const int firstHitNode = firstTrackNode + static_cast<int>(tracks.size());
  const int sink = firstHitNode + static_cast<int>(hits.size());
  const int numberOfNodes = sink + 1;
  std::vector<std::vector<MatchEdge> > graph(numberOfNodes);
  for (std::size_t trackIndex = 0; trackIndex < tracks.size(); ++trackIndex)
    AddFlowEdge(graph, source, firstTrackNode + static_cast<int>(trackIndex), 0.);
  for (std::size_t hitIndex = 0; hitIndex < hits.size(); ++hitIndex)
    AddFlowEdge(graph, firstHitNode + static_cast<int>(hitIndex), sink, 0.);

  for (std::size_t trackIndex = 0; trackIndex < tracks.size(); ++trackIndex) {
    for (std::size_t hitIndex = 0; hitIndex < hits.size(); ++hitIndex) {
      const double maximumDistance = hits[hitIndex].detector == kItofDetector
                                   ? maximumItofDistanceCm
                                   : maximumEtofDistanceCm;
      if (!(maximumDistance >= 0.) || !std::isfinite(maximumDistance))
        continue;
      const double distance = ProjectionDistanceCm(tracks[trackIndex],
                                                   hits[hitIndex]);
      if (!std::isfinite(distance) || distance > maximumDistance) continue;
      // Detector-normalized cost. A tiny stable term resolves exact ties
      // without changing the physical ordering at meaningful precision.
      const double normalized = maximumDistance > 0.
                              ? distance / maximumDistance : 0.;
      const double tie = 1.e-12 *
        (static_cast<double>(trackIndex) * (hits.size() + 1u) + hitIndex);
      AddFlowEdge(graph,
                  firstTrackNode + static_cast<int>(trackIndex),
                  firstHitNode + static_cast<int>(hitIndex),
                  normalized + tie, static_cast<int>(trackIndex),
                  static_cast<int>(hitIndex), distance);
    }
  }

  // Successive shortest augmenting paths on the unit-capacity residual graph.
  // Running until no path remains gives maximum cardinality; shortest paths,
  // including negative reverse arcs, give minimum cost at that cardinality.
  const double infinity = std::numeric_limits<double>::infinity();
  while (true) {
    std::vector<double> distance(numberOfNodes, infinity);
    std::vector<int> previousNode(numberOfNodes, -1);
    std::vector<int> previousEdge(numberOfNodes, -1);
    std::vector<bool> inQueue(numberOfNodes, false);
    std::deque<int> queue;
    distance[source] = 0.;
    queue.push_back(source);
    inQueue[source] = true;
    while (!queue.empty()) {
      const int node = queue.front();
      queue.pop_front();
      inQueue[node] = false;
      for (std::size_t edgeIndex = 0; edgeIndex < graph[node].size(); ++edgeIndex) {
        const MatchEdge& edge = graph[node][edgeIndex];
        if (edge.capacity <= 0) continue;
        const double candidate = distance[node] + edge.cost;
        if (candidate + 1.e-15 >= distance[edge.to]) continue;
        distance[edge.to] = candidate;
        previousNode[edge.to] = node;
        previousEdge[edge.to] = static_cast<int>(edgeIndex);
        if (!inQueue[edge.to]) {
          queue.push_back(edge.to);
          inQueue[edge.to] = true;
        }
      }
    }
    if (previousNode[sink] < 0) break;
    for (int node = sink; node != source; node = previousNode[node]) {
      MatchEdge& edge = graph[previousNode[node]][previousEdge[node]];
      --edge.capacity;
      ++graph[node][edge.reverse].capacity;
    }
  }

  std::vector<TofMatch> selected;
  for (std::size_t trackIndex = 0; trackIndex < tracks.size(); ++trackIndex) {
    const int node = firstTrackNode + static_cast<int>(trackIndex);
    for (std::size_t edgeIndex = 0; edgeIndex < graph[node].size(); ++edgeIndex) {
      const MatchEdge& edge = graph[node][edgeIndex];
      if (edge.trackIndex < 0 || edge.hitIndex < 0 || edge.capacity != 0) continue;
      const std::size_t hitIndex = static_cast<std::size_t>(edge.hitIndex);
      TofMatch match;
      match.trackIndex = trackIndex;
      match.hitIndex = hitIndex;
      match.trackId = tracks[trackIndex].id;
      match.hitId = hits[hitIndex].id;
      match.detector = hits[hitIndex].detector;
      match.distanceCm = edge.distanceCm;
      selected.push_back(match);
      break;
    }
  }
  std::sort(selected.begin(), selected.end(), TrackOrder());
  return selected;
}

T0SelectionConfig::T0SelectionConfig()
  : source(kFixedT0), aggregation(kMedianT0), fixedTimeNs(0.),
    rawToNsScale(1.), offsetNs(0.), sigmaNs(0.10)
{
}

T0SelectionResult::T0SelectionResult()
  : valid(false), numberOfHits(0),
    selectedRawTime(std::numeric_limits<double>::quiet_NaN()),
    timeNs(std::numeric_limits<double>::quiet_NaN()), sigmaNs(0.)
{
}

T0SelectionResult SelectEventT0(const std::vector<double>& rawHitTimes,
                                const T0SelectionConfig& config)
{
  T0SelectionResult result;
  result.numberOfHits = static_cast<int>(rawHitTimes.size());
  result.sigmaNs = config.sigmaNs;
  if (!(config.sigmaNs >= 0.) || !std::isfinite(config.sigmaNs)) return result;

  if (config.source == kFixedT0) {
    if (!std::isfinite(config.fixedTimeNs)) return result;
    result.valid = true;
    result.timeNs = config.fixedTimeNs;
    return result;
  }

  std::vector<double> finiteTimes;
  for (std::size_t i = 0; i < rawHitTimes.size(); ++i)
    if (std::isfinite(rawHitTimes[i])) finiteTimes.push_back(rawHitTimes[i]);
  if (finiteTimes.empty() || !(config.rawToNsScale > 0.) ||
      !std::isfinite(config.rawToNsScale) ||
      !std::isfinite(config.offsetNs)) return result;

  double selected = finiteTimes.front();
  if (config.aggregation == kEarliestT0) {
    selected = *std::min_element(finiteTimes.begin(), finiteTimes.end());
  }
  else if (config.aggregation == kMeanT0) {
    selected = 0.;
    for (std::size_t i = 0; i < finiteTimes.size(); ++i)
      selected += finiteTimes[i];
    selected /= static_cast<double>(finiteTimes.size());
  }
  else if (config.aggregation == kMedianT0) {
    std::sort(finiteTimes.begin(), finiteTimes.end());
    const std::size_t middle = finiteTimes.size() / 2;
    selected = finiteTimes[middle];
    if ((finiteTimes.size() % 2u) == 0u)
      selected = 0.5 * (finiteTimes[middle - 1] + finiteTimes[middle]);
  }

  result.selectedRawTime = selected;
  result.timeNs = selected * config.rawToNsScale + config.offsetNs;
  result.valid = std::isfinite(result.timeNs);
  return result;
}

const char* TofDetectorName(TofDetector detector)
{
  if (detector == kItofDetector) return "iTOF";
  if (detector == kEtofDetector) return "eTOF";
  return "none";
}

} // namespace CeeV0Analysis
