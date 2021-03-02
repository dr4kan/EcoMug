////////////////////////////////////////////////////////////////////////////
// EcoMug: Efficient COsmic MUon Generator                                //
// Copyright (C) 2021 Davide Pagano <davide.pagano@unibs.it>              //
// EcoMug is based on the following work:                                 //
//                 arXiv:XXXX.XXXX                                        //
//                                                                        //
// This program is free software: you can redistribute it and/or modify   //
// it under the terms of the GNU General Public License as published by   //
// the Free Software Foundation, either version 3 of the License, or      //
// (at your option) any later version.                                    //
//                                                                        //
// This program is distributed in the hope that it will be useful,        //
// but WITHOUT ANY WARRANTY; without even the implied warranty of         //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
// GNU General Public License for more details.                           //
//                                                                        //
// You should have received a copy of the GNU General Public License      //
// along with this program.  If not, see <https://www.gnu.org/licenses/>. //
////////////////////////////////////////////////////////////////////////////

#ifndef EcoMug_H
#define EcoMug_H

#include <math.h>
#include <array>
#include <random>

///////////////////////////////////////////////////////////////
// Fast generation of random numbers, based on xoroshiro128+.
///////////////////////////////////////////////////////////////
class EMRandom {
public:
  EMRandom() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(0, std::numeric_limits<uint64_t>::max());
    s[0] = dis(gen);
    s[1] = dis(gen);
  };

  double GenerateRandomDouble() {
    uint64_t x = next();
    return to_double(x);
  };

private:
  int64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
  };

  uint64_t next() {
    const uint64_t s0 = s[0];
    uint64_t s1 = s[1];
    const uint64_t result = s0 + s1;
    s1 ^= s0;
    s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14);
    s[1] = rotl(s1, 36);
    return result;
  };

  double to_double(uint64_t x) {
    const union { uint64_t i; double d; } u = { .i = UINT64_C(0x3FF) << 52 | x >> 12 };
    return u.d - 1.0;
  };

  uint64_t s[2];
};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////
// Generation of cosmic muons position, direction and momentum
// from a plane, cylinder and half-sphere
///////////////////////////////////////////////////////////////
class EcoMug {
public:
  enum EMGeometry {Sky, Cylinder, HSphere};

private:
  EMGeometry mGenMethod;
  std::array<double, 3> mGenerationPosition;
  double mGenerationTheta;
  double mGenerationPhi;
  double mGenerationMomentum;
  double mMinimumMomentum;
  double mMaximumMomentum;
  double mMinimumTheta;
  double mMaximumTheta;
  double mMinimumPhi;
  double mMaximumPhi;
  std::array<double, 2> mSkySize;
  std::array<double, 3> mSkyCenterPosition;
  double mCylinderHeight;
  double mCylinderRadius;
  std::array<double, 3> mCylinderCenterPosition;
  double mHSphereRadius;
  double mMaxFuncSkyCylinder;
  std::array<double, 3> mHSphereCenterPosition;
  EMRandom mRandom;

public:
  // Default constructor
  EcoMug() : mGenMethod(Sky),
  mGenerationPosition({{0., 0., 0.}}), mGenerationTheta(0.), mGenerationPhi(0.),
  mGenerationMomentum(0.), mMinimumMomentum(0.01), mMaximumMomentum(1000.),
  mMinimumTheta(0.), mMaximumTheta(M_PI/2.), mMinimumPhi(0.), mMaximumPhi(2.*M_PI),
  mSkySize({{0., 0.}}), mSkyCenterPosition({{0., 0., 0.}}), mCylinderHeight(0.),
  mCylinderRadius(0.), mCylinderCenterPosition({{0., 0., 0.}}), mHSphereCenterPosition({{0., 0., 0.}}),
  mMaxFuncSkyCylinder(5.3176)
  {};

  ///////////////////////////////////////////////////////////////
  // Methods to access the parameters of the generated muon
  ///////////////////////////////////////////////////////////////
  // Get the generation position
  const std::array<double, 3>& GetGenerationPosition() const {
    return mGenerationPosition;
  }
  // Get the generation momentum
  double GetGenerationMomentum() const {
    return mGenerationMomentum;
  }
  // Get the generation theta
  double GetGenerationTheta() const {
    return mGenerationTheta;
  }
  // Get the generation phi
  double GetGenerationPhi() const {
    return mGenerationPhi;
  }
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Methods for the geometry of the generation
  ///////////////////////////////////////////////////////////////
  // Set generation from sky
  void SetUseSky() {
    mGenMethod = Sky;
  };
  // Set cylindrical generation
  void SetUseCylinder() {
    mGenMethod = Cylinder;
  };
  // Set half-sphere generation
  void SetUseHSphere() {
    mGenMethod = HSphere;
  };
  // Set the generation method (Sky, Cylinder or HSphere)
  void SetGenerationMethod(EMGeometry genM) {
    mGenMethod = genM;
  };

  // Get the generation method (Sky, Cylinder or HSphere)
  EMGeometry GetGenerationMethod() const {
    return mGenMethod;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Common methods to all geometries
  ///////////////////////////////////////////////////////////////
  // Set minimum generation Momentum
  void SetMinimumMomentum(double momentum) {
    mMinimumMomentum = momentum;
    ComputeMaxSkyCylinder();
  };
  // Set minimum generation Momentum
  void SetMaximumMomentum(double momentum) {
    mMaximumMomentum = momentum;
  };
  // Set minimum generation Theta
  void SetMinimumTheta(double theta) {
    mMinimumTheta = theta;
  };
  // Set minimum generation Theta
  void SetMaximumTheta(double theta) {
    mMaximumTheta = theta;
  };
  // Set minimum generation Phi
  void SetMinimumPhi(double phi) {
    mMinimumPhi = phi;
  };
  // Set minimum generation Theta
  void SetMaximumPhi(double phi) {
    mMaximumPhi = phi;
  };

  // Get minimum generation Momentum
  double GetMinimumMomentum() const {
    return mMinimumMomentum;
  };
  // Get minimum generation Momentum
  double GetMaximumMomentum() const {
    return mMaximumMomentum;
  };
  // Set minimum generation Theta
  double GetMinimumTheta() const {
    return mMinimumTheta;
  };
  // Set minimum generation Theta
  double GetMaximumTheta() const {
    return mMaximumTheta;
  };
  // Set minimum generation Phi
  double GetMinimumPhi() const {
    return mMinimumPhi;
  };
  // Set minimum generation Theta
  double GetMaximumPhi() const {
    return mMaximumPhi;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Methods for the plane-based generation
  ///////////////////////////////////////////////////////////////
  // Set sky size
  void SetSkySize(const std::array<double, 2>& size) {
    mSkySize = size;
  };

  // Set sky center position
  void SetSkyCenterPosition(const std::array<double, 3>& position) {
    mSkyCenterPosition = position;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Methods for the cylinder-based generation
  ///////////////////////////////////////////////////////////////
  // Set cylinder radius
  void SetCylinderRadius(double radius) {
    mCylinderRadius = radius;
  };
  // Set cylinder height
  void SetCylinderHeight(double height) {
    mCylinderHeight = height;
  };
  // Set cylinder center position
  void SetCylinderCenterPosition(const std::array<double, 3>& position) {
    mCylinderCenterPosition = position;
  };

  // Get cylinder radius
  double GetCylinderRadius() const {
    return mCylinderRadius;
  };
  // Get cylinder height
  double GetCylinderHeight() const {
    return mCylinderHeight;
  };
  // Get cylinder center position
  const std::array<double, 3>& GetCylinderCenterPosition() const {
    return mCylinderCenterPosition;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Methods for the half sphere-based generation
  ///////////////////////////////////////////////////////////////
  // Set sphere radius
  void SetHSphereRadius(double radius) {
    mHSphereRadius = radius;
  };

  // Set half-sphere center position
  void SetHSphereCenterPosition(const std::array<double, 3>& position) {
    mHSphereCenterPosition = position;
  };

  // Get the sphere radius
  double GetHSphereRadius() const {
    return mHSphereRadius;
  };

  // Get half-sphere center position
  const std::array<double, 3>& GetHSphereCenterPosition() const {
    return mHSphereCenterPosition;
  };
  ///////////////////////////////////////////////////////////////


private:
  void ComputeMaxSkyCylinder() {
    if (mMinimumMomentum < 0.448) return;
    double n = 2.856-0.655*log(mMinimumMomentum);
    mMaxFuncSkyCylinder = 1600*Pow(mMinimumMomentum+2.68, -3.175)*Pow(mMinimumMomentum, 0.279)*Pow(Cos(0.6968), n+1)*Sin(0.6968);
  };

  double Sin(double x) const {
    double x2 = x * x;
    return (((((-2.05342856289746600727e-08*x2 + 2.70405218307799040084e-06)*x2
    - 1.98125763417806681909e-04)*x2 + 8.33255814755188010464e-03)*x2
    - 1.66665772196961623983e-01)*x2 + 9.99999707044156546685e-01)*x;
  };

  double Cos(double x) const {
    double x2 = x * x;
    return ((((-2.21941782786353727022e-07*x2 + 2.42532401381033027481e-05)*x2
    - 1.38627507062573673756e-03)*x2 + 4.16610337354021107429e-02)*x2
    - 4.99995582499065048420e-01)*x2 + 9.99999443739537210853e-01;
  };

  double Pow(double a, double b) const {
    union {
      double n;
      int c[2];
    } unum = { a };
    unum.c[1] = (int)(b * (unum.c[1] - 1072632447) + 1072632447);
    unum.c[0] = 0;
    return unum.n;
  };

  void GeneratePositionSky() {
    mGenerationPosition[0] = mSkySize[0]*mRandom.GenerateRandomDouble() + mSkyCenterPosition[0]-mSkySize[0]/2.;
    mGenerationPosition[1] = mSkySize[1]*mRandom.GenerateRandomDouble() + mSkyCenterPosition[1]-mSkySize[1]/2.;
    mGenerationPosition[2] = mSkyCenterPosition[2];
  };

  double GeneratePositionCylinder() {
    double phi0            = 2.*M_PI*mRandom.GenerateRandomDouble();
    mGenerationPosition[0] = mCylinderCenterPosition[0] + mCylinderRadius*Cos(phi0);
    mGenerationPosition[1] = mCylinderCenterPosition[1] + mCylinderRadius*Sin(phi0);
    mGenerationPosition[2] = mCylinderCenterPosition[2] + mCylinderHeight*mRandom.GenerateRandomDouble() - mCylinderHeight/2.;
    return phi0;
  };

public:
  ///////////////////////////////////////////////////////////////
  // Generate a cosmic muon
  ///////////////////////////////////////////////////////////////
  void Generate() {
    bool accepted = false;
    double r1t    = 0.;
    double r1m    = 0.;
    double r1p    = 0.;
    double r2     = 0.;
    double n      = 0.;
    double ftheta = 0.;
    double fphi   = 0.;

    // Sky or cylinder generation
    if (mGenMethod == Sky || mGenMethod == Cylinder) {
      // Generation of the momentum and theta angle
      while (!accepted) {
        r1t = mRandom.GenerateRandomDouble();
        r1m = mRandom.GenerateRandomDouble();
        r2  = mRandom.GenerateRandomDouble();
        mGenerationTheta = (mMaximumTheta - mMinimumTheta)*r1t + mMinimumTheta;
        mGenerationMomentum = (mMaximumMomentum - mMinimumMomentum)*r1m + mMinimumMomentum;
        n = 2.856-0.655*log(mGenerationMomentum);
        if (n < 0.1) n = 0.1;

        if (mGenMethod == Sky) {
          ftheta = 1600*Pow(mGenerationMomentum+2.68, -3.175)*Pow(mGenerationMomentum, 0.279)*Pow(Cos(mGenerationTheta), n+1)*Sin(mGenerationTheta);
          if (mMaxFuncSkyCylinder*r2 < ftheta) accepted = true;
        }

        if(mGenMethod == Cylinder)  {
          ftheta = 1600*Pow(mGenerationMomentum+2.68, -3.175)*Pow(mGenerationMomentum, 0.279)*Pow(Cos(mGenerationTheta), n)*Pow(Sin(mGenerationTheta), 2);
          if (mMaxFuncSkyCylinder*r2 < ftheta) accepted = true;
        }
      }
      mGenerationTheta = M_PI - mGenerationTheta;

      // Generation of the position and phi angle
      if (mGenMethod == Sky) {
        GeneratePositionSky();
        r1p = mRandom.GenerateRandomDouble();
        mGenerationPhi = (mMaximumPhi - mMinimumPhi)*r1p + mMinimumPhi;
      }
      if (mGenMethod == Cylinder) {
        double phi0 = GeneratePositionCylinder();
        accepted = false;
        while (!accepted) {
          r1p = mRandom.GenerateRandomDouble();
          r2  = mRandom.GenerateRandomDouble();
          mGenerationPhi = (mMaximumPhi - mMinimumPhi)*r1p + mMinimumPhi;
          fphi = fabs(Cos(mGenerationPhi));
          if (r2 < fphi) {
            accepted = true;
          }
        }
        mGenerationPhi = mGenerationPhi + phi0;
        if (mGenerationPhi >= 2.*M_PI) mGenerationPhi -= 2.*M_PI;

        // Check if the muon is inward
        if (Sin(mGenerationTheta)*Cos(mGenerationPhi)*mGenerationPosition[0] + Sin(mGenerationTheta)*Sin(mGenerationPhi)*mGenerationPosition[1] > 0) Generate();
      }
    }

    // Half-sphere generation
    if (mGenMethod == HSphere) {
      // Generation point on the half-sphere
      double phi0   = 2.*M_PI*mRandom.GenerateRandomDouble();
      double theta0 = 0.;
      double r0t    = 0.;
      while (!accepted) {
        r0t = mRandom.GenerateRandomDouble();
        theta0 = acos(r0t);
        r1t = mRandom.GenerateRandomDouble();
        mGenerationTheta = (mMaximumTheta - mMinimumTheta)*r1t + mMinimumTheta;
        r1m = mRandom.GenerateRandomDouble();
        mGenerationMomentum = (mMaximumMomentum - mMinimumMomentum)*r1m + mMinimumMomentum;
        r1p = mRandom.GenerateRandomDouble();
        mGenerationPhi = (mMaximumPhi - mMinimumPhi)*r1p + mMinimumPhi;
        r2  = mRandom.GenerateRandomDouble();
        n = 2.856-0.655*log(mGenerationMomentum);
        if (n < 0.1) n = 0.1;
        ftheta = 1600*Pow(mGenerationMomentum+2.68, -3.175)*Pow(mGenerationMomentum, 3.175-2.896)*Pow(Cos(mGenerationTheta), n)*fabs(Sin(mGenerationTheta)*Sin(theta0)*Cos(mGenerationPhi)+Cos(mGenerationTheta)*Cos(theta0))*Sin(mGenerationTheta);
        if (11*r2 < ftheta) accepted = true;
      }

      mGenerationPosition[0] = mHSphereCenterPosition[0] + mHSphereRadius*Sin(theta0)*Sin(phi0);
      mGenerationPosition[1] = mHSphereCenterPosition[1] + mHSphereRadius*Cos(theta0)-mHSphereRadius/2.;
      mGenerationPosition[2] = mHSphereCenterPosition[2] + mHSphereRadius*Sin(theta0)*Cos(phi0);
      mGenerationTheta = M_PI - mGenerationTheta;
      mGenerationPhi = mGenerationPhi + phi0;
      if (mGenerationPhi >= 2.*M_PI) mGenerationPhi -= 2.*M_PI;
      mGenerationPhi += M_PI;
      if (mGenerationPhi >= 2.*M_PI) mGenerationPhi -= 2.*M_PI;
    }
  };
  ///////////////////////////////////////////////////////////////

};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


#endif
