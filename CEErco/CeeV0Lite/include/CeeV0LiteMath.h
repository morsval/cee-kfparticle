#ifndef CEE_V0_LITE_MATH_H
#define CEE_V0_LITE_MATH_H

#include <cstddef>

namespace CeeV0Lite
{

struct Vec3
{
  Vec3();
  Vec3(double x, double y, double z);

  double& operator[](std::size_t i) { return v[i]; }
  const double& operator[](std::size_t i) const { return v[i]; }

  double v[3];
};

struct Mat3
{
  Mat3();
  static Mat3 Identity(double scale = 1.);
  double m[3][3];
};

struct Mat6
{
  Mat6();
  static Mat6 Identity(double scale = 1.);
  double m[6][6];
};

// Covariance of the STAR-like mother state
//   (x, y, z, px, py, pz, E).
// The mixed units are therefore cm^2, cm*GeV and GeV^2.
struct Mat7
{
  Mat7();
  static Mat7 Identity(double scale = 1.);
  double m[7][7];
};

// Neutral production/decay state (x,y,z,px,py,pz,E,S), S in cm/(GeV/c).
struct Mat8
{
  Mat8();
  double m[8][8];
};

Vec3 operator+(const Vec3& a, const Vec3& b);
Vec3 operator-(const Vec3& a, const Vec3& b);
Vec3 operator*(double scale, const Vec3& value);
Vec3 operator*(const Vec3& value, double scale);
Vec3 operator/(const Vec3& value, double scale);

double Dot(const Vec3& a, const Vec3& b);
Vec3 Cross(const Vec3& a, const Vec3& b);
double Norm2(const Vec3& value);
double Norm(const Vec3& value);
bool Normalize(const Vec3& input, Vec3& output);
bool IsFinite(const Vec3& value);

Vec3 Multiply(const Mat3& matrix, const Vec3& vector);
Mat3 Multiply(const Mat3& a, const Mat3& b);
Mat3 Transpose(const Mat3& matrix);
Mat3 operator+(const Mat3& a, const Mat3& b);
Mat3 operator-(const Mat3& a, const Mat3& b);
Mat3 operator*(double scale, const Mat3& matrix);

double QuadraticForm(const Vec3& vector, const Mat3& matrix);
bool InvertSymmetric3(const Mat3& input, Mat3& inverse,
                      double minimumDeterminant = 1.e-30);
bool InvertPositiveDefinite3(const Mat3& input, Mat3& inverse,
                             double relativePivotTolerance = 1.e-14);
bool InvertSymmetric2(const double input[2][2], double inverse[2][2],
                      double minimumDeterminant = 1.e-30);

Mat3 PositionBlock(const Mat6& covariance);
Mat3 MomentumBlock(const Mat6& covariance);
void SetPositionBlock(Mat6& covariance, const Mat3& block);
void SetMomentumBlock(Mat6& covariance, const Mat3& block);
void Symmetrize(Mat3& matrix);
void Symmetrize(Mat6& matrix);
void Symmetrize(Mat7& matrix);
void Symmetrize(Mat8& matrix);
bool IsFinite(const Mat8& matrix);
bool IsPositiveSemidefinite(const Mat8& matrix,
                            double relativeTolerance = 1.e-10);
double QuadraticForm(const double vector[8], const Mat8& matrix);
bool IsPositiveSemidefinite(const Mat6& matrix,
                            double relativeTolerance = 1.e-10);
bool IsFinite(const Mat7& matrix);
bool IsPositiveSemidefinite(const Mat7& matrix,
                            double relativeTolerance = 1.e-10);
double QuadraticForm(const double vector[7], const Mat7& matrix);

void MakeTransverseBasis(const Vec3& direction, Vec3& first, Vec3& second);

} // namespace CeeV0Lite

#endif
