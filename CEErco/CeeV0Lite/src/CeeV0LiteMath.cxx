#include "CeeV0LiteMath.h"

#include <algorithm>
#include <cmath>

namespace CeeV0Lite
{

Vec3::Vec3()
  : v{0., 0., 0.}
{
}

Vec3::Vec3(double x, double y, double z)
  : v{x, y, z}
{
}

Mat3::Mat3()
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      m[i][j] = 0.;
}

Mat3 Mat3::Identity(double scale)
{
  Mat3 out;
  for (int i = 0; i < 3; ++i) out.m[i][i] = scale;
  return out;
}

Mat6::Mat6()
{
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      m[i][j] = 0.;
}

Mat6 Mat6::Identity(double scale)
{
  Mat6 out;
  for (int i = 0; i < 6; ++i) out.m[i][i] = scale;
  return out;
}

Mat7::Mat7()
{
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      m[i][j] = 0.;
}

Mat7 Mat7::Identity(double scale)
{
  Mat7 out;
  for (int i = 0; i < 7; ++i) out.m[i][i] = scale;
  return out;
}

Vec3 operator+(const Vec3& a, const Vec3& b)
{
  return Vec3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vec3 operator-(const Vec3& a, const Vec3& b)
{
  return Vec3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vec3 operator*(double scale, const Vec3& value)
{
  return Vec3(scale * value[0], scale * value[1], scale * value[2]);
}

Vec3 operator*(const Vec3& value, double scale)
{
  return scale * value;
}

Vec3 operator/(const Vec3& value, double scale)
{
  return Vec3(value[0] / scale, value[1] / scale, value[2] / scale);
}

double Dot(const Vec3& a, const Vec3& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 Cross(const Vec3& a, const Vec3& b)
{
  return Vec3(a[1] * b[2] - a[2] * b[1],
              a[2] * b[0] - a[0] * b[2],
              a[0] * b[1] - a[1] * b[0]);
}

double Norm2(const Vec3& value)
{
  return Dot(value, value);
}

double Norm(const Vec3& value)
{
  return std::sqrt(Norm2(value));
}

bool Normalize(const Vec3& input, Vec3& output)
{
  const double norm = Norm(input);
  if (!(norm > 1.e-14) || !std::isfinite(norm)) return false;
  output = input / norm;
  return true;
}

bool IsFinite(const Vec3& value)
{
  return std::isfinite(value[0]) && std::isfinite(value[1]) &&
         std::isfinite(value[2]);
}

Vec3 Multiply(const Mat3& matrix, const Vec3& vector)
{
  Vec3 out;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out[i] += matrix.m[i][j] * vector[j];
  return out;
}

Mat3 Multiply(const Mat3& a, const Mat3& b)
{
  Mat3 out;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        out.m[i][j] += a.m[i][k] * b.m[k][j];
  return out;
}

Mat3 Transpose(const Mat3& matrix)
{
  Mat3 out;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out.m[i][j] = matrix.m[j][i];
  return out;
}

Mat3 operator+(const Mat3& a, const Mat3& b)
{
  Mat3 out;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out.m[i][j] = a.m[i][j] + b.m[i][j];
  return out;
}

Mat3 operator-(const Mat3& a, const Mat3& b)
{
  Mat3 out;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out.m[i][j] = a.m[i][j] - b.m[i][j];
  return out;
}

Mat3 operator*(double scale, const Mat3& matrix)
{
  Mat3 out;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out.m[i][j] = scale * matrix.m[i][j];
  return out;
}

double QuadraticForm(const Vec3& vector, const Mat3& matrix)
{
  return Dot(vector, Multiply(matrix, vector));
}

bool InvertSymmetric3(const Mat3& a, Mat3& inverse, double minimumDeterminant)
{
  const double c00 = a.m[1][1] * a.m[2][2] - a.m[1][2] * a.m[2][1];
  const double c01 = a.m[0][2] * a.m[2][1] - a.m[0][1] * a.m[2][2];
  const double c02 = a.m[0][1] * a.m[1][2] - a.m[0][2] * a.m[1][1];
  const double determinant = a.m[0][0] * c00 + a.m[1][0] * c01 + a.m[2][0] * c02;
  if (!std::isfinite(determinant) || std::fabs(determinant) < minimumDeterminant)
    return false;

  const double invDet = 1. / determinant;
  inverse.m[0][0] = c00 * invDet;
  inverse.m[0][1] = c01 * invDet;
  inverse.m[0][2] = c02 * invDet;
  inverse.m[1][0] = (a.m[1][2] * a.m[2][0] - a.m[1][0] * a.m[2][2]) * invDet;
  inverse.m[1][1] = (a.m[0][0] * a.m[2][2] - a.m[0][2] * a.m[2][0]) * invDet;
  inverse.m[1][2] = (a.m[0][2] * a.m[1][0] - a.m[0][0] * a.m[1][2]) * invDet;
  inverse.m[2][0] = (a.m[1][0] * a.m[2][1] - a.m[1][1] * a.m[2][0]) * invDet;
  inverse.m[2][1] = (a.m[0][1] * a.m[2][0] - a.m[0][0] * a.m[2][1]) * invDet;
  inverse.m[2][2] = (a.m[0][0] * a.m[1][1] - a.m[0][1] * a.m[1][0]) * invDet;
  Symmetrize(inverse);
  return true;
}

bool InvertPositiveDefinite3(const Mat3& input, Mat3& inverse,
                             double relativePivotTolerance)
{
  if (!(relativePivotTolerance >= 0.) ||
      !std::isfinite(relativePivotTolerance)) return false;
  double scale = 0.;
  for (int i = 0; i < 3; ++i) {
    if (!std::isfinite(input.m[i][i])) return false;
    scale = std::max(scale, std::fabs(input.m[i][i]));
  }
  if (!(scale > 0.)) return false;

  double lower[3][3] = {{0.}};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j <= i; ++j) {
      double value = input.m[i][j];
      for (int k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
      if (i == j) {
        if (!(value > relativePivotTolerance * scale) ||
            !std::isfinite(value)) return false;
        lower[i][j] = std::sqrt(value);
      }
      else {
        lower[i][j] = value / lower[j][j];
      }
    }
  }

  // Solve L L^T X = I, one column at a time.
  inverse = Mat3();
  for (int column = 0; column < 3; ++column) {
    double forward[3] = {0.};
    for (int i = 0; i < 3; ++i) {
      double value = (i == column) ? 1. : 0.;
      for (int k = 0; k < i; ++k) value -= lower[i][k] * forward[k];
      forward[i] = value / lower[i][i];
    }
    double solution[3] = {0.};
    for (int i = 2; i >= 0; --i) {
      double value = forward[i];
      for (int k = i + 1; k < 3; ++k)
        value -= lower[k][i] * solution[k];
      solution[i] = value / lower[i][i];
    }
    for (int i = 0; i < 3; ++i) inverse.m[i][column] = solution[i];
  }
  Symmetrize(inverse);
  return true;
}

bool InvertSymmetric2(const double input[2][2], double inverse[2][2],
                      double minimumDeterminant)
{
  const double determinant = input[0][0] * input[1][1] -
                             input[0][1] * input[1][0];
  if (!std::isfinite(determinant) || std::fabs(determinant) < minimumDeterminant)
    return false;
  const double invDet = 1. / determinant;
  inverse[0][0] = input[1][1] * invDet;
  inverse[0][1] = -input[0][1] * invDet;
  inverse[1][0] = -input[1][0] * invDet;
  inverse[1][1] = input[0][0] * invDet;
  return true;
}

Mat3 PositionBlock(const Mat6& covariance)
{
  Mat3 out;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out.m[i][j] = covariance.m[i][j];
  return out;
}

Mat3 MomentumBlock(const Mat6& covariance)
{
  Mat3 out;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out.m[i][j] = covariance.m[i + 3][j + 3];
  return out;
}

void SetPositionBlock(Mat6& covariance, const Mat3& block)
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      covariance.m[i][j] = block.m[i][j];
}

void SetMomentumBlock(Mat6& covariance, const Mat3& block)
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      covariance.m[i + 3][j + 3] = block.m[i][j];
}

void Symmetrize(Mat3& matrix)
{
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < i; ++j) {
      const double value = 0.5 * (matrix.m[i][j] + matrix.m[j][i]);
      matrix.m[i][j] = value;
      matrix.m[j][i] = value;
    }
  }
}

void Symmetrize(Mat6& matrix)
{
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < i; ++j) {
      const double value = 0.5 * (matrix.m[i][j] + matrix.m[j][i]);
      matrix.m[i][j] = value;
      matrix.m[j][i] = value;
    }
  }
}

void Symmetrize(Mat7& matrix)
{
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < i; ++j) {
      const double value = 0.5 * (matrix.m[i][j] + matrix.m[j][i]);
      matrix.m[i][j] = value;
      matrix.m[j][i] = value;
    }
  }
}

bool IsPositiveSemidefinite(const Mat6& matrix, double relativeTolerance)
{
  if (!(relativeTolerance >= 0.) || !std::isfinite(relativeTolerance))
    return false;
  double scale[6] = {0.};
  for (int i = 0; i < 6; ++i) {
    if (!std::isfinite(matrix.m[i][i]) || matrix.m[i][i] < 0.) return false;
    scale[i] = std::sqrt(matrix.m[i][i]);
  }
  double normalized[6][6] = {{0.}};
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      if (!std::isfinite(matrix.m[i][j])) return false;
      if (scale[i] == 0. || scale[j] == 0.) {
        if (std::fabs(matrix.m[i][j]) > relativeTolerance) return false;
      }
      else {
        normalized[i][j] = matrix.m[i][j] / (scale[i] * scale[j]);
      }
    }
  }
  double lower[6][6] = {{0.}};
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j <= i; ++j) {
      double value = normalized[i][j];
      for (int k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
      if (i == j) {
        if (value < -relativeTolerance) return false;
        lower[i][j] = value > relativeTolerance ? std::sqrt(value) : 0.;
      }
      else if (lower[j][j] > 0.) {
        lower[i][j] = value / lower[j][j];
      }
      else if (std::fabs(value) > relativeTolerance) {
        return false;
      }
    }
  }
  return true;
}

bool IsFinite(const Mat7& matrix)
{
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      if (!std::isfinite(matrix.m[i][j])) return false;
  return true;
}

bool IsPositiveSemidefinite(const Mat7& matrix, double relativeTolerance)
{
  if (!IsFinite(matrix) || !(relativeTolerance >= 0.) ||
      !std::isfinite(relativeTolerance)) return false;

  // Scale to a correlation-like matrix before checking, because Mat7 mixes
  // cm and GeV.  A zero-variance coordinate must also have a zero row/column.
  double scale[7] = {0.};
  for (int i = 0; i < 7; ++i) {
    if (matrix.m[i][i] < 0.) return false;
    scale[i] = std::sqrt(matrix.m[i][i]);
  }
  double normalized[7][7] = {{0.}};
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      if (scale[i] == 0. || scale[j] == 0.) {
        if (std::fabs(matrix.m[i][j]) > relativeTolerance) return false;
        normalized[i][j] = 0.;
      }
      else {
        normalized[i][j] = matrix.m[i][j] / (scale[i] * scale[j]);
      }
    }
  }

  // Semidefinite Cholesky with a small negative-pivot tolerance.  A zero
  // pivot may only have numerically zero remaining couplings.
  double lower[7][7] = {{0.}};
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j <= i; ++j) {
      double value = normalized[i][j];
      for (int k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
      if (i == j) {
        if (value < -relativeTolerance) return false;
        lower[i][j] = value > relativeTolerance ? std::sqrt(value) : 0.;
      }
      else if (lower[j][j] > 0.) {
        lower[i][j] = value / lower[j][j];
      }
      else if (std::fabs(value) > relativeTolerance) {
        return false;
      }
    }
  }
  return true;
}

double QuadraticForm(const double vector[7], const Mat7& matrix)
{
  double value = 0.;
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      value += vector[i] * matrix.m[i][j] * vector[j];
  return value;
}

Mat8::Mat8()
{
  for (int i=0;i<8;++i)
    for (int j=0;j<8;++j) m[i][j]=0.;
}

void Symmetrize(Mat8& matrix)
{
  for (int i=0;i<8;++i)
    for (int j=0;j<i;++j)
      matrix.m[i][j]=matrix.m[j][i]=0.5*(matrix.m[i][j]+matrix.m[j][i]);
}

bool IsFinite(const Mat8& matrix)
{
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j)
      if (!std::isfinite(matrix.m[i][j])) return false;
  return true;
}

bool IsPositiveSemidefinite(const Mat8& matrix, double relativeTolerance)
{
  if (!IsFinite(matrix) || !(relativeTolerance >= 0.) ||
      !std::isfinite(relativeTolerance)) return false;

  // Scale to a correlation-like matrix before checking, because Mat8 mixes
  // cm and GeV.  A zero-variance coordinate must also have a zero row/column.
  double scale[8] = {0.};
  for (int i = 0; i < 8; ++i) {
    if (matrix.m[i][i] < 0.) return false;
    scale[i] = std::sqrt(matrix.m[i][i]);
  }
  double normalized[8][8] = {{0.}};
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      if (scale[i] == 0. || scale[j] == 0.) {
        if (std::fabs(matrix.m[i][j]) > relativeTolerance) return false;
        normalized[i][j] = 0.;
      }
      else {
        normalized[i][j] = matrix.m[i][j] / (scale[i] * scale[j]);
      }
    }
  }

  // Semidefinite Cholesky with a small negative-pivot tolerance.  A zero
  // pivot may only have numerically zero remaining couplings.
  double lower[8][8] = {{0.}};
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j <= i; ++j) {
      double value = normalized[i][j];
      for (int k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
      if (i == j) {
        if (value < -relativeTolerance) return false;
        lower[i][j] = value > relativeTolerance ? std::sqrt(value) : 0.;
      }
      else if (lower[j][j] > 0.) {
        lower[i][j] = value / lower[j][j];
      }
      else if (std::fabs(value) > relativeTolerance) {
        return false;
      }
    }
  }
  return true;
}

double QuadraticForm(const double vector[8], const Mat8& matrix)
{
  double value = 0.;
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j)
      value += vector[i] * matrix.m[i][j] * vector[j];
  return value;
}

void MakeTransverseBasis(const Vec3& inputDirection, Vec3& first, Vec3& second)
{
  Vec3 direction;
  if (!Normalize(inputDirection, direction)) {
    first = Vec3(1., 0., 0.);
    second = Vec3(0., 1., 0.);
    return;
  }

  const Vec3 reference = (std::fabs(direction[2]) < 0.8)
    ? Vec3(0., 0., 1.) : Vec3(1., 0., 0.);
  Normalize(Cross(direction, reference), first);
  Normalize(Cross(direction, first), second);
}

} // namespace CeeV0Lite
