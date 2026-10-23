#include "CeeV0LiteVertexFitter.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <iostream>

namespace
{
bool Close(double a, double b, double tolerance)
{
  return std::fabs(a - b) <= tolerance;
}

double Uniform01()
{
  static std::uint64_t state = 0x1234abcd9876ULL;
  state = state * 6364136223846793005ULL + 1442695040888963407ULL;
  return static_cast<double>((state >> 11) & ((1ULL << 53) - 1)) /
         static_cast<double>(1ULL << 53);
}

double StandardNormal()
{
  static bool hasSpare = false;
  static double spare = 0.;
  if (hasSpare) {
    hasSpare = false;
    return spare;
  }
  const double u1 = std::max(1.e-15, Uniform01());
  const double u2 = Uniform01();
  const double radius = std::sqrt(-2. * std::log(u1));
  const double angle = 6.2831853071795864769 * u2;
  spare = radius * std::sin(angle);
  hasSpare = true;
  return radius * std::cos(angle);
}

void SmearState(CeeV0Lite::TrackState& state)
{
  const double sigmaPosition[3] = {0.03, 0.04, 0.05};
  const double sigmaMomentum[3] = {0.004, 0.005, 0.006};
  const double correlation = 0.25;
  const double independentScale = std::sqrt(1. - correlation * correlation);
  for (int i = 0; i < 3; ++i) {
    const double common = StandardNormal();
    const double independent = StandardNormal();
    state.positionCm[i] += sigmaPosition[i] * common;
    state.momentumGeV[i] += sigmaMomentum[i] *
      (correlation * common + independentScale * independent);
  }
}

CeeV0Lite::TrackState MakeState(const CeeV0Lite::Vec3& vertex,
                                const CeeV0Lite::Vec3& directionInput,
                                double distanceCm,
                                double momentumGeV,
                                double massGeV,
                                int charge)
{
  CeeV0Lite::Vec3 direction;
  assert(CeeV0Lite::Normalize(directionInput, direction));
  CeeV0Lite::TrackState state;
  state.positionCm = vertex + distanceCm * direction;
  state.momentumGeV = momentumGeV * direction;
  state.massGeV = massGeV;
  state.charge = charge;

  const double sigmaPosition[3] = {0.03, 0.04, 0.05};
  const double sigmaMomentum[3] = {0.004, 0.005, 0.006};
  for (int i = 0; i < 3; ++i) {
    state.covariance.m[i][i] = sigmaPosition[i] * sigmaPosition[i];
    state.covariance.m[3 + i][3 + i] =
      sigmaMomentum[i] * sigmaMomentum[i];
    // Non-zero x-p correlations make the common-SV update feed back into the
    // fitted daughter momenta, exercising the correlation that the old
    // independent sigma(M) formula could not retain.
    const double cross = 0.25 * sigmaPosition[i] * sigmaMomentum[i];
    state.covariance.m[i][3 + i] = cross;
    state.covariance.m[3 + i][i] = cross;
  }
  return state;
}

void AssertMotherEqual(const CeeV0Lite::MotherState& first,
                       const CeeV0Lite::MotherState& second,
                       double covarianceTolerance = 1.e-10)
{
  assert(first.covarianceValid && second.covarianceValid);
  assert(CeeV0Lite::Norm(first.positionCm - second.positionCm) < 1.e-10);
  assert(CeeV0Lite::Norm(first.momentumGeV - second.momentumGeV) < 1.e-10);
  assert(Close(first.energyGeV, second.energyGeV, 1.e-11));
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      assert(Close(first.covariance.m[i][j],
                   second.covariance.m[i][j], covarianceTolerance));
}

}

int main()
{
  const double protonMass = 0.9382720813;
  const double pionMass = 0.13957039;
  const CeeV0Lite::Vec3 vertex(0.4, -0.2, 4.0);
  const CeeV0Lite::TrackState proton = MakeState(
    vertex, CeeV0Lite::Vec3(0.3, 0.1, 1.), 20., 1.0,
    protonMass, +1);
  CeeV0Lite::TrackState pion = MakeState(
    vertex, CeeV0Lite::Vec3(-0.5, 0.2, 1.), 15., 0.35,
    pionMass, -1);
  // Deliberately make the measured lines slightly skew; an exactly
  // intersecting toy would hide sequential/order and DCA-residual bugs.
  pion.positionCm[1] += 0.012;

  CeeV0Lite::VertexFitter fitter;
  CeeV0Lite::VertexFitResult fit;
  assert(fitter.Fit(proton, pion, 0., fit));
  assert(fitter.RecalculateMotherCovariance(proton, pion, 0., fit.mother));
  assert(fit.mother.covarianceValid);
  assert(fit.mother.covarianceStatus == CeeV0Lite::kMotherCovarianceValid);
  assert(fit.mother.jacobianFitsAttempted == 24);
  assert(fit.mother.jacobianFitsSucceeded == 24);
  assert(fit.mother.jacobianFailureCount == 0);
  assert(fit.mother.jacobianPlusFailureMask == 0u);
  assert(fit.mother.jacobianMinusFailureMask == 0u);
  assert(fit.mother.jacobianFailedColumn == -1);
  assert(fit.mother.jacobianFailedSide == 0);
  assert(CeeV0Lite::IsFinite(fit.mother.covariance));
  for (int i = 0; i < 7; ++i) {
    assert(fit.mother.covariance.m[i][i] >= -1.e-14);
    for (int j = 0; j < 7; ++j)
      assert(Close(fit.mother.covariance.m[i][j],
                   fit.mother.covariance.m[j][i], 1.e-14));
  }

  double mass = 0.;
  double sigmaMass = 0.;
  assert(CeeV0Lite::MotherMassAndError(fit.mother, mass, sigmaMass));
  assert(mass > 0. && sigmaMass > 0.);
  const double expectedMass = std::sqrt(
    fit.mother.energyGeV * fit.mother.energyGeV -
    CeeV0Lite::Norm2(fit.mother.momentumGeV));
  assert(Close(mass, expectedMass, 1.e-14));

  // A nominal vertex can be valid while one central-difference refit is not.
  // With axis-aligned unit momenta and h_p=1 GeV, the -h perturbation of px1
  // (column 3) and py2 (column 10) makes that daughter's momentum exactly
  // zero.  The implementation must try all 24 sides, record both failures and
  // invalidate the final covariance without erasing the nominal mother state.
  const CeeV0Lite::TrackState axisFirst = MakeState(
    vertex, CeeV0Lite::Vec3(1., 0., 0.), 10., 1., protonMass, +1);
  const CeeV0Lite::TrackState axisSecond = MakeState(
    vertex, CeeV0Lite::Vec3(0., 1., 0.), 10., 1., pionMass, -1);
  CeeV0Lite::VertexFitterConfig failureConfig;
  failureConfig.motherJacobianMinimumMomentumStepGeV = 1.;
  CeeV0Lite::VertexFitter failureFitter(failureConfig);
  CeeV0Lite::VertexFitResult failureFit;
  assert(failureFitter.Fit(axisFirst, axisSecond, 0., failureFit));
  const CeeV0Lite::Vec3 nominalFailurePosition =
    failureFit.mother.positionCm;
  const CeeV0Lite::Vec3 nominalFailureMomentum =
    failureFit.mother.momentumGeV;
  const double nominalFailureEnergy = failureFit.mother.energyGeV;
  assert(!failureFitter.RecalculateMotherCovariance(
    axisFirst, axisSecond, 0., failureFit.mother));
  assert(failureFit.mother.covarianceStatus ==
         CeeV0Lite::kMotherCovarianceJacobianPerturbationFailed);
  assert(!failureFit.mother.covarianceValid);
  assert(failureFit.mother.jacobianFitsAttempted == 24);
  assert(failureFit.mother.jacobianFitsSucceeded == 22);
  assert(failureFit.mother.jacobianFailureCount == 2);
  assert((failureFit.mother.jacobianMinusFailureMask & (1u << 3)) != 0u);
  assert((failureFit.mother.jacobianMinusFailureMask & (1u << 10)) != 0u);
  assert(failureFit.mother.jacobianPlusFailureMask == 0u);
  assert(failureFit.mother.jacobianFailedColumn == 3);
  assert(failureFit.mother.jacobianFailedSide == -1);
  assert(failureFit.mother.jacobianFirstFailureFitStatus ==
         static_cast<int>(CeeV0Lite::kVertexFitNumericalFailure));
  assert(CeeV0Lite::Norm(failureFit.mother.positionCm -
                         nominalFailurePosition) == 0.);
  assert(CeeV0Lite::Norm(failureFit.mother.momentumGeV -
                         nominalFailureMomentum) == 0.);
  assert(failureFit.mother.energyGeV == nominalFailureEnergy);
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      assert(std::isnan(failureFit.mother.covariance.m[i][j]));

  // Swapping complete daughter states must not change the mother.  This is a
  // useful regression test for sequential-update order bias.
  CeeV0Lite::VertexFitResult swapped;
  assert(fitter.Fit(pion, proton, 0., swapped));
  assert(fitter.RecalculateMotherCovariance(pion, proton, 0.,
                                             swapped.mother));
  AssertMotherEqual(fit.mother, swapped.mother);

  // Moving a zero-field reference point along the same geometrical line must
  // only change the fitted path parameter, not the mother estimator or its
  // outer-Jacobian covariance.
  CeeV0Lite::TrackState slidProton;
  CeeV0Lite::TransportConfig transportConfig;
  assert(CeeV0Lite::PropagateState(proton, 7., 0., transportConfig,
                                   slidProton));
  CeeV0Lite::VertexFitResult slid;
  assert(fitter.Fit(slidProton, pion, 0., slid));
  assert(fitter.RecalculateMotherCovariance(slidProton, pion, 0.,
                                             slid.mother));
  AssertMotherEqual(fit.mother, slid.mother, 5.e-7);

  // The symmetric linearization must also remain daughter-order invariant in
  // the finite-field helix path.
  const CeeV0Lite::TrackState protonAtVertex = MakeState(
    vertex, CeeV0Lite::Vec3(0.3, 0.1, 1.), 0., 1.0,
    protonMass, +1);
  const CeeV0Lite::TrackState pionAtVertex = MakeState(
    vertex, CeeV0Lite::Vec3(-0.5, 0.2, 1.), 0., 0.35,
    pionMass, -1);
  CeeV0Lite::TrackState fieldProton;
  CeeV0Lite::TrackState fieldPion;
  assert(CeeV0Lite::PropagateState(protonAtVertex, 20., 0.5,
                                   transportConfig, fieldProton));
  assert(CeeV0Lite::PropagateState(pionAtVertex, 15., 0.5,
                                   transportConfig, fieldPion));
  fieldPion.positionCm[1] += 0.005;
  CeeV0Lite::VertexFitResult fieldFit;
  CeeV0Lite::VertexFitResult fieldSwapped;
  assert(fitter.Fit(fieldProton, fieldPion, 0.5, fieldFit));
  assert(fitter.RecalculateMotherCovariance(fieldProton, fieldPion, 0.5,
                                             fieldFit.mother));
  assert(fitter.Fit(fieldPion, fieldProton, 0.5, fieldSwapped));
  assert(fitter.RecalculateMotherCovariance(fieldPion, fieldProton, 0.5,
                                             fieldSwapped.mother));
  AssertMotherEqual(fieldFit.mother, fieldSwapped.mother, 1.e-8);

  // The same reference-surface invariance must hold on a helix, not only on a
  // straight zero-field line.
  CeeV0Lite::TrackState slidFieldProton;
  assert(CeeV0Lite::PropagateState(fieldProton, 7., 0.5, transportConfig,
                                   slidFieldProton));
  CeeV0Lite::VertexFitResult slidFieldFit;
  assert(fitter.Fit(slidFieldProton, fieldPion, 0.5, slidFieldFit));
  assert(fitter.RecalculateMotherCovariance(slidFieldProton, fieldPion, 0.5,
                                             slidFieldFit.mother));
  AssertMotherEqual(fieldFit.mother, slidFieldFit.mother, 5.e-8);

  // A common scalar nuisance adds a rank-one matrix exactly once.
  CeeV0Lite::MotherState withCommon = fit.mother;
  const double derivative[7] = {0.01, -0.02, 0.03,
                                0.2, -0.1, 0.3, 0.05};
  const double variance = 0.04;
  assert(CeeV0Lite::AddCommonNuisanceCovariance(
    withCommon, derivative, variance));
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      assert(Close(withCommon.covariance.m[i][j] -
                   fit.mother.covariance.m[i][j],
                   derivative[i] * variance * derivative[j], 1.e-14));

  // Significantly negative propagated variance is an error, not something to
  // hide with an unconditional max(0, variance).
  CeeV0Lite::MotherState invalid = fit.mother;
  invalid.covariance = CeeV0Lite::Mat7();
  invalid.covariance.m[3][3] = -1.;
  invalid.covarianceValid = true;
  assert(!CeeV0Lite::MotherMassAndError(invalid, mass, sigmaMass));

  CeeV0Lite::MotherState indefinite = fit.mother;
  indefinite.covariance = CeeV0Lite::Mat7();
  indefinite.covariance.m[3][3] = 1.;
  indefinite.covariance.m[6][6] = 1.;
  indefinite.covariance.m[3][6] = 2.;
  indefinite.covariance.m[6][3] = 2.;
  indefinite.covarianceValid = true;
  assert(!CeeV0Lite::IsPositiveSemidefinite(indefinite.covariance));
  assert(!CeeV0Lite::MotherMassAndError(indefinite, mass, sigmaMass));

  // Small-error toy-MC closure of the complete estimator.  This test would
  // fail for the former fixed-path/independent-daughter sigma(M), which
  // noticeably underestimated the sampled width.
  CeeV0Lite::VertexFitterConfig closureConfig;
  closureConfig.measurementVarianceFloorCm2 = 0.;
  CeeV0Lite::VertexFitter closureFitter(closureConfig, transportConfig);
  CeeV0Lite::VertexFitResult closureNominal;
  assert(closureFitter.Fit(proton, pion, 0., closureNominal));
  assert(closureFitter.RecalculateMotherCovariance(
    proton, pion, 0., closureNominal.mother));
  double predictedMass = 0.;
  double predictedSigma = 0.;
  assert(CeeV0Lite::MotherMassAndError(closureNominal.mother,
                                       predictedMass, predictedSigma));
  const int samples = 500;
  double sumMass = 0.;
  double sumMass2 = 0.;
  for (int sample = 0; sample < samples; ++sample) {
    CeeV0Lite::TrackState smearedProton = proton;
    CeeV0Lite::TrackState smearedPion = pion;
    SmearState(smearedProton);
    SmearState(smearedPion);
    CeeV0Lite::VertexFitResult sampledFit;
    assert(closureFitter.Fit(smearedProton, smearedPion, 0., sampledFit));
    const double sampledMass2 = sampledFit.mother.energyGeV *
      sampledFit.mother.energyGeV -
      CeeV0Lite::Norm2(sampledFit.mother.momentumGeV);
    assert(sampledMass2 > 0.);
    const double sampledMass = std::sqrt(sampledMass2);
    sumMass += sampledMass;
    sumMass2 += sampledMass * sampledMass;
  }
  const double sampleMean = sumMass / samples;
  const double sampleVariance = (sumMass2 - samples * sampleMean * sampleMean) /
                                (samples - 1);
  assert(sampleVariance > 0.);
  const double sampledSigma = std::sqrt(sampleVariance);
  assert(std::fabs(sampledSigma / predictedSigma - 1.) < 0.20);

  // Independent process-noise modes are not part of the outer 12D input yet;
  // the API must reject that configuration explicitly rather than understate
  // the covariance.
  CeeV0Lite::TransportConfig noisyTransport = transportConfig;
  noisyTransport.positionProcessNoiseCm2PerCm = 1.e-6;
  CeeV0Lite::VertexFitter noisyFitter(closureConfig, noisyTransport);
  CeeV0Lite::VertexFitResult noisyFit;
  assert(noisyFitter.Fit(proton, pion, 0., noisyFit));
  assert(!noisyFitter.RecalculateMotherCovariance(
    proton, pion, 0., noisyFit.mother));
  assert(noisyFit.mother.covarianceStatus ==
         CeeV0Lite::kMotherCovarianceProcessNoiseUnsupported);

  std::cout << "STAR-like mother covariance and shared-nuisance tests passed"
            << std::endl;
  return 0;
}
