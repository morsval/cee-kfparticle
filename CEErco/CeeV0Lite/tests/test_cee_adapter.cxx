#include "CeeV0LiteCeeAdapter.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace
{

struct MockVector3
{
  MockVector3() : x(0.), y(0.), z(0.) {}
  MockVector3(double inX, double inY, double inZ) : x(inX), y(inY), z(inZ) {}
  double X() const { return x; }
  double Y() const { return y; }
  double Z() const { return z; }
  double x;
  double y;
  double z;
};

class MockMatrix6
{
 public:
  MockMatrix6()
  {
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j) value[i][j] = 0.;
  }
  int GetNrows() const { return 6; }
  double& operator()(int row, int column) { return value[row][column]; }
  double operator()(int row, int column) const { return value[row][column]; }
  double value[6][6];
};

struct MockTrack
{
  MockTrack()
    : hasSeed(false), converged(true), x(1.), y(2.), z(3.),
      px(0.), py(1.), pz(0.), charge(-1.), chi2(12.), ndf(6.)
  {
    for (int i = 0; i < 6; ++i) seedCovariance(i, i) = 0.01 * (i + 1);
  }

  bool HasVertexSeed() const { return hasSeed; }
  const MockVector3& GetVertexSeedPos() const { return seedPosition; }
  const MockVector3& GetVertexSeedMom() const { return seedMomentum; }
  void GetVertexSeedCov(MockMatrix6& output) const { output = seedCovariance; }
  bool GetIsFitConverged() const { return converged; }
  double GetX() const { return x; }
  double GetY() const { return y; }
  double GetZ() const { return z; }
  double GetPx() const { return px; }
  double GetPy() const { return py; }
  double GetPz() const { return pz; }
  double GetFitCharge() const { return charge; }
  double GetFitChi2() const { return chi2; }
  double GetFitNdf() const { return ndf; }

  bool hasSeed;
  bool converged;
  MockVector3 seedPosition;
  MockVector3 seedMomentum;
  MockMatrix6 seedCovariance;
  double x;
  double y;
  double z;
  double px;
  double py;
  double pz;
  double charge;
  double chi2;
  double ndf;
};

struct MockVertex
{
  const std::vector<double>& GetX() const { return x; }
  const std::vector<double>& GetY() const { return y; }
  const std::vector<double>& GetZ() const { return z; }
  const std::vector<double>& GetXErr() const { return xError; }
  const std::vector<double>& GetYErr() const { return yError; }
  const std::vector<double>& GetZErr() const { return zError; }
  const std::vector<int>& GetRank() const { return rank; }
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> xError;
  std::vector<double> yError;
  std::vector<double> zError;
  std::vector<int> rank;
};

bool Near(double a, double b, double tolerance = 1.e-12)
{
  return std::fabs(a - b) < tolerance;
}

void TestPrimaryVertexSelection()
{
  MockVertex source;
  source.x.push_back(1.);
  source.x.push_back(4.);
  source.x.push_back(std::numeric_limits<double>::quiet_NaN());
  source.y.push_back(2.);
  source.y.push_back(5.);
  source.y.push_back(8.);
  source.z.push_back(3.);
  source.z.push_back(6.);
  source.z.push_back(9.);
  source.rank.push_back(2);
  source.rank.push_back(8);
  source.rank.push_back(100);
  source.xError.push_back(0.01);
  source.xError.push_back(0.02);
  source.yError.push_back(0.03);
  source.yError.push_back(0.04);
  source.zError.push_back(0.05); // index 1 deliberately falls back

  CeeV0Lite::PrimaryVertex output;
  CeeV0Lite::CeeAdapterStatus status = CeeV0Lite::kCeeAdapterNoValidVertex;
  std::size_t index = 99u;
  assert(CeeV0Lite::MakePrimaryVertex(source, output, &index,
                                      CeeV0Lite::CeeAdapterConfig(), &status));
  assert(status == CeeV0Lite::kCeeAdapterSuccess);
  assert(index == 1u);
  assert(Near(output.positionCm[0], 4.));
  assert(Near(output.positionCm[1], 5.));
  assert(Near(output.positionCm[2], 6.));
  assert(Near(output.covarianceCm2.m[0][0], 0.02 * 0.02));
  assert(Near(output.covarianceCm2.m[1][1], 0.04 * 0.04));
  assert(Near(output.covarianceCm2.m[2][2], 0.15 * 0.15));

  MockVertex empty;
  assert(!CeeV0Lite::MakePrimaryVertex(empty, output, 0,
                                       CeeV0Lite::CeeAdapterConfig(),
                                       &status));
  assert(status == CeeV0Lite::kCeeAdapterNoValidVertex);
}

void TestZeroFieldTrackInput()
{
  MockTrack source;
  source.hasSeed = true;
  source.seedPosition = MockVector3(10., 20., 30.);
  source.seedMomentum = MockVector3(3., 4., 0.);

  CeeV0Lite::CeeTofMeasurement tof;
  tof.available = true;
  tof.positionCm = CeeV0Lite::Vec3(100., 200., 300.);
  tof.timeNs = 12.;
  CeeV0Lite::TrackInput output;
  CeeV0Lite::CeeAdapterStatus status;

  CeeV0Lite::CeeTrackReferenceState reference;
  assert(CeeV0Lite::SelectTrackReferenceState(
    source, CeeV0Lite::kCeeZeroFieldTofMode, reference,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(status == CeeV0Lite::kCeeAdapterSuccessVertexSeed);
  assert(CeeV0Lite::CeeAdapterStatusIsSuccess(status));
  assert(reference.usedVertexSeed);
  assert(reference.hasTrackingCovariance);
  assert(Near(reference.positionCm[0], 10.));
  assert(Near(reference.direction[1], 4.));

  assert(CeeV0Lite::MakeTrackInput(
    source, 17, CeeV0Lite::kCeeZeroFieldTofMode, tof,
    CeeV0Lite::kProtonHypothesis, output,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(status == CeeV0Lite::kCeeAdapterSuccessVertexSeed);
  assert(output.id == 17);
  assert(Near(output.positionCm[0], 10.));
  assert(Near(output.direction[0], 3.));
  assert(!output.hasTrackingMomentum);
  assert(!output.hasTrackingCovariance);
  assert(output.hasTof);
  assert(!output.chargeKnown && output.charge == 0);
  assert(Near(output.qualityChi2, 12.));
  assert(output.qualityNdf == 6);

  // No vertex seed is allowed only in zero-field mode, and the complete
  // legacy state plus configured covariance is used (no reference mixing).
  source.hasSeed = false;
  assert(CeeV0Lite::SelectTrackReferenceState(
    source, CeeV0Lite::kCeeZeroFieldTofMode, reference,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(status == CeeV0Lite::kCeeAdapterSuccessLegacyFallback);
  assert(CeeV0Lite::CeeAdapterStatusIsSuccess(status));
  assert(!reference.usedVertexSeed);
  assert(!reference.hasTrackingCovariance);
  assert(Near(reference.positionCm[0], source.x));
  assert(Near(reference.direction[1], source.py));

  assert(CeeV0Lite::MakeTrackInput(
    source, 18, CeeV0Lite::kCeeZeroFieldTofMode, tof,
    CeeV0Lite::kPionHypothesis, output,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(Near(output.positionCm[0], source.x));
  assert(Near(output.direction[1], source.py));
  assert(status == CeeV0Lite::kCeeAdapterSuccessLegacyFallback);

  tof.available = false;
  assert(!CeeV0Lite::MakeTrackInput(
    source, 18, CeeV0Lite::kCeeZeroFieldTofMode, tof,
    CeeV0Lite::kPionHypothesis, output,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(status == CeeV0Lite::kCeeAdapterMissingTof);
}

void TestMagneticFieldTrackInput()
{
  MockTrack source;
  CeeV0Lite::CeeTofMeasurement noTof;
  CeeV0Lite::TrackInput output;
  CeeV0Lite::CeeAdapterStatus status;

  assert(!CeeV0Lite::MakeTrackInput(
    source, 1, CeeV0Lite::kCeeMagneticFieldTrackingMode, noTof,
    CeeV0Lite::kPionHypothesis, output,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(status == CeeV0Lite::kCeeAdapterMissingVertexSeed);

  source.hasSeed = true;
  source.seedPosition = MockVector3(0.1, 0.2, 0.3);
  source.seedMomentum = MockVector3(0.4, 0.5, 0.6);
  assert(CeeV0Lite::MakeTrackInput(
    source, 2, CeeV0Lite::kCeeMagneticFieldTrackingMode, noTof,
    CeeV0Lite::kPionHypothesis, output,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(output.hasTrackingMomentum);
  assert(output.hasTrackingCovariance);
  assert(!output.hasTof);
  assert(Near(output.trackingMomentumGeV[2], 0.6));
  assert(Near(output.trackingCovariance.m[5][5], 0.06));

  source.seedCovariance(5, 5) = -1.;
  assert(!CeeV0Lite::MakeTrackInput(
    source, 3, CeeV0Lite::kCeeMagneticFieldTrackingMode, noTof,
    CeeV0Lite::kPionHypothesis, output,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(status == CeeV0Lite::kCeeAdapterInvalidTrackingCovariance);

  // The zero-field helper and MakeTrackInput must make the same fallback
  // decision when a seed exists but its covariance is unusable.  This is the
  // state TOF matching must consume.
  CeeV0Lite::CeeTrackReferenceState reference;
  assert(CeeV0Lite::SelectTrackReferenceState(
    source, CeeV0Lite::kCeeZeroFieldTofMode, reference,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(!reference.usedVertexSeed);
  assert(status == CeeV0Lite::kCeeAdapterSuccessLegacyFallback);

  CeeV0Lite::CeeTofMeasurement tof;
  tof.available = true;
  tof.positionCm = CeeV0Lite::Vec3(100., 0., 0.);
  tof.timeNs = 10.;
  assert(CeeV0Lite::MakeTrackInput(
    source, 4, CeeV0Lite::kCeeZeroFieldTofMode, tof,
    CeeV0Lite::kPionHypothesis, output,
    CeeV0Lite::CeeAdapterConfig(), &status));
  assert(status == CeeV0Lite::kCeeAdapterSuccessLegacyFallback);
  assert(Near(output.positionCm[0], reference.positionCm[0]));
  assert(Near(output.positionCm[1], reference.positionCm[1]));
  assert(Near(output.positionCm[2], reference.positionCm[2]));
  assert(Near(output.direction[0], reference.direction[0]));
  assert(Near(output.direction[1], reference.direction[1]));
  assert(Near(output.direction[2], reference.direction[2]));
}

} // namespace

int main()
{
  TestPrimaryVertexSelection();
  TestZeroFieldTrackInput();
  TestMagneticFieldTrackInput();
  std::cout << "CEE adapter tests passed\n";
  return 0;
}
