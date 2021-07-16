/////////////////////////////////////////////////////////////////////////////////////
// EcoMug: Efficient COsmic MUon Generator                                         //
// Copyright (C) 2021 Davide Pagano <davide.pagano@unibs.it>                       //
// EcoMug is based on the following work:                                          //
// authors,                                                                        //
// "EcoMug: an Efficient COsmic MUon Generator for cosmic-ray muons applications", //
// arXiv:XXXX.XXXX                                                                 //
//                                                                                 //
// This program is free software: you can redistribute it and/or modify            //
// it under the terms of the GNU General Public License as published by            //
// the Free Software Foundation, either version 3 of the License, or               //
// (at your option) any later version.                                             //
//                                                                                 //
// This program is distributed in the hope that it will be useful,                 //
// but WITHOUT ANY WARRANTY; without even the implied warranty of                  //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   //
// GNU General Public License for more details.                                    //
//                                                                                 //
// You should have received a copy of the GNU General Public License               //
// along with this program.  If not, see <https://www.gnu.org/licenses/>.          //
/////////////////////////////////////////////////////////////////////////////////////

#ifndef EcoMug_H
#define EcoMug_H

#include <math.h>
#include <array>
#include <random>

//! Fast generation of random numbers
//! This class is based on the xoroshiro128+ generator.
//! https://prng.di.unimi.it/
class EMRandom {
public:
  EMRandom() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(0, std::numeric_limits<uint64_t>::max());
    s[0] = dis(gen);
    s[1] = dis(gen);
  };

  void SetSeed(uint64_t seed) {
    s[0] = seed;
    s[1] = seed;
  };

  double GenerateRandomDouble() {
    uint64_t x = next();
    return to_double(x);
  };

  double GenerateRandomDouble(double x1, double x2) {
    return (x2-x1)*GenerateRandomDouble()+x1;
  };

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
    union U {
      uint64_t i;
      double d;
    };
    U u = { UINT64_C(0x3FF) << 52 | x >> 12 };
    return u.d - 1.0;
  };

  uint64_t s[2];
};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


//! Class for the generation of cosmic muons
class EcoMug {
  friend class EMRandom;

public:
  /// Possible generation methods
  enum EMGeometry {
    Sky,      ///< generation from a plane (flat sky)
    Cylinder, ///< generation from a cylinder
    HSphere   ///< generation from a half-sphere
  };

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
  int    mCharge;
  std::array<double, 2> mSkySize;
  std::array<double, 3> mSkyCenterPosition;
  double mCylinderHeight;
  double mCylinderRadius;
  std::array<double, 3> mCylinderCenterPosition;
  double mHSphereRadius;
  double mMaxFuncSkyCylinder;
  std::array<double, 3> mHSphereCenterPosition;
  EMRandom mRandom;
  std::default_random_engine mEngineC;
  std::discrete_distribution<int> mDiscDistC;

public:
  // Default constructor
  EcoMug() : mGenMethod(Sky),
  mGenerationPosition({{0., 0., 0.}}), mGenerationTheta(0.), mGenerationPhi(0.),
  mGenerationMomentum(0.), mMinimumMomentum(0.01), mMaximumMomentum(1000.),
  mMinimumTheta(0.), mMaximumTheta(M_PI/2.), mMinimumPhi(0.), mMaximumPhi(2.*M_PI),
  mCharge(1), mSkySize({{0., 0.}}), mSkyCenterPosition({{0., 0., 0.}}), mCylinderHeight(0.),
  mCylinderRadius(0.), mCylinderCenterPosition({{0., 0., 0.}}), mHSphereRadius(0.),
  mMaxFuncSkyCylinder(5.3176), mHSphereCenterPosition({{0., 0., 0.}}),
  mEngineC(std::random_device{}()) {
    mDiscDistC = std::discrete_distribution<int>({128, 100});
  };

  ///////////////////////////////////////////////////////////////
  // Methods to access the parameters of the generated muon
  ///////////////////////////////////////////////////////////////
  /// Get the generation position
  const std::array<double, 3>& GetGenerationPosition() const {
    return mGenerationPosition;
  };
  /// Get the generation momentum
  double GetGenerationMomentum() const {
    return mGenerationMomentum;
  };
  /// Get the generation theta
  double GetGenerationTheta() const {
    return mGenerationTheta;
  };
  /// Get the generation phi
  double GetGenerationPhi() const {
    return mGenerationPhi;
  };
  /// Get charge
  int GetCharge() const {
    return mCharge;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Methods for the geometry of the generation
  ///////////////////////////////////////////////////////////////
  /// Set generation from sky
  void SetUseSky() {
    mGenMethod = Sky;
  };
  /// Set cylindrical generation
  void SetUseCylinder() {
    mGenMethod = Cylinder;
  };
  /// Set half-sphere generation
  void SetUseHSphere() {
    mGenMethod = HSphere;
  };
  /// Set the generation method (Sky, Cylinder or HSphere)
  void SetGenerationMethod(EMGeometry genM) {
    mGenMethod = genM;
  };

  /// Get the generation method (Sky, Cylinder or HSphere)
  EMGeometry GetGenerationMethod() const {
    return mGenMethod;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Common methods to all geometries
  ///////////////////////////////////////////////////////////////
  /// Set the seed for the internal PRNG (if 0 a random seed is used)
  void SetSeed(uint64_t seed) {
    if (seed > 0) mRandom.SetSeed(seed);
  };

  /// Set minimum generation Momentum
  void SetMinimumMomentum(double momentum) {
    mMinimumMomentum = momentum;
  };
  /// Set maximum generation Momentum
  void SetMaximumMomentum(double momentum) {
    mMaximumMomentum = momentum;
  };
  /// Set minimum generation Theta
  void SetMinimumTheta(double theta) {
    mMinimumTheta = theta;
  };
  /// Set maximum generation Theta
  void SetMaximumTheta(double theta) {
    mMaximumTheta = theta;
  };
  /// Set minimum generation Phi
  void SetMinimumPhi(double phi) {
    mMinimumPhi = phi;
  };
  /// Set maximum generation Phi
  void SetMaximumPhi(double phi) {
    mMaximumPhi = phi;
  };

  /// Get minimum generation Momentum
  double GetMinimumMomentum() const {
    return mMinimumMomentum;
  };
  /// Get maximum generation Momentum
  double GetMaximumMomentum() const {
    return mMaximumMomentum;
  };
  /// Get minimum generation Theta
  double GetMinimumTheta() const {
    return mMinimumTheta;
  };
  /// Get maximum generation Theta
  double GetMaximumTheta() const {
    return mMaximumTheta;
  };
  /// Get minimum generation Phi
  double GetMinimumPhi() const {
    return mMinimumPhi;
  };
  /// Get maximum generation Phi
  double GetMaximumPhi() const {
    return mMaximumPhi;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Methods for the plane-based generation
  ///////////////////////////////////////////////////////////////
  /// Set sky size
  void SetSkySize(const std::array<double, 2>& size) {
    mSkySize = size;
  };

  /// Set sky center position
  void SetSkyCenterPosition(const std::array<double, 3>& position) {
    mSkyCenterPosition = position;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Methods for the cylinder-based generation
  ///////////////////////////////////////////////////////////////
  /// Set cylinder radius
  void SetCylinderRadius(double radius) {
    mCylinderRadius = radius;
  };
  /// Set cylinder height
  void SetCylinderHeight(double height) {
    mCylinderHeight = height;
  };
  /// Set cylinder center position
  void SetCylinderCenterPosition(const std::array<double, 3>& position) {
    mCylinderCenterPosition = position;
  };

  /// Get cylinder radius
  double GetCylinderRadius() const {
    return mCylinderRadius;
  };
  /// Get cylinder height
  double GetCylinderHeight() const {
    return mCylinderHeight;
  };
  /// Get cylinder center position
  const std::array<double, 3>& GetCylinderCenterPosition() const {
    return mCylinderCenterPosition;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  // Methods for the half sphere-based generation
  ///////////////////////////////////////////////////////////////
  /// Set half-sphere radius
  void SetHSphereRadius(double radius) {
    mHSphereRadius = radius;
  };
  /// Set half-sphere center position
  void SetHSphereCenterPosition(const std::array<double, 3>& position) {
    mHSphereCenterPosition = position;
  };
  /// Get half-sphere radius
  double GetHSphereRadius() const {
    return mHSphereRadius;
  };
  /// Get half-sphere center position
  const std::array<double, 3>& GetHSphereCenterPosition() const {
    return mHSphereCenterPosition;
  };
  ///////////////////////////////////////////////////////////////


private:
  double F1Cumulative(double x) {
    return 1. - 8.534790171171021/pow(x + 2.68, 87./40.);
  };
  double F1Inverse(double x) {
    return (2.68 - 2.68*pow(1. - x, 40./87.))/pow(1. - x, 40./87.);
  };

  double GenerateMomentumF1() {
    double z = mRandom.GenerateRandomDouble(F1Cumulative(mMinimumMomentum), F1Cumulative(mMaximumMomentum));
    return F1Inverse(z);
  };

  void GeneratePositionSky() {
    mGenerationPosition[0] = mRandom.GenerateRandomDouble(mSkyCenterPosition[0]-mSkySize[0]/2., mSkyCenterPosition[0]+mSkySize[0]/2.);
    mGenerationPosition[1] = mRandom.GenerateRandomDouble(mSkyCenterPosition[1]-mSkySize[1]/2., mSkyCenterPosition[1]+mSkySize[1]/2.);
    mGenerationPosition[2] = mSkyCenterPosition[2];
  };

  double GeneratePositionCylinder() {
    double phi0            = 2.*M_PI*mRandom.GenerateRandomDouble();
    mGenerationPosition[0] = mCylinderCenterPosition[0] + mCylinderRadius*cos(phi0);
    mGenerationPosition[1] = mCylinderCenterPosition[1] + mCylinderRadius*sin(phi0);
    mGenerationPosition[2] = mRandom.GenerateRandomDouble(mCylinderCenterPosition[2]-mCylinderHeight/2., mCylinderCenterPosition[2]+mCylinderHeight/2.);
    return phi0;
  };

public:
  ///////////////////////////////////////////////////////////////
  /// Generate a cosmic muon
  ///////////////////////////////////////////////////////////////
  void Generate() {
    bool accepted = false;
    double r2     = 0.;
    double n      = 0.;
    double ftheta = 0.;
    double fphi   = 0.;

    // Sky or cylinder generation
    if (mGenMethod == Sky || mGenMethod == Cylinder) {
      // Generation of the momentum and theta angle
      while (!accepted) {
        r2  = mRandom.GenerateRandomDouble();
        mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationMomentum = GenerateMomentumF1();
        n = 2.856-0.655*log(mGenerationMomentum);
        if (n < 0.1) n = 0.1;

        if (mGenMethod == Sky) {
          ftheta = 1600*pow(mGenerationMomentum, 0.279)*pow(cos(mGenerationTheta), n+1)*sin(mGenerationTheta);
          if (2797*r2 < ftheta) accepted = true;
        }

        if(mGenMethod == Cylinder)  {
          ftheta = 1600*pow(mGenerationMomentum, 0.279)*pow(cos(mGenerationTheta), n)*pow(sin(mGenerationTheta), 2);
          if (4730*r2 < ftheta) accepted = true;
        }
      }
      mGenerationTheta = M_PI - mGenerationTheta;

      // Generation of the position and phi angle
      if (mGenMethod == Sky) {
        GeneratePositionSky();
        mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      }
      if (mGenMethod == Cylinder) {
        double phi0 = GeneratePositionCylinder();
        accepted = false;
        while (!accepted) {
          r2  = mRandom.GenerateRandomDouble();
          mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
          fphi = fabs(cos(mGenerationPhi));
          if (r2 < fphi) accepted = true;
        }
        mGenerationPhi = mGenerationPhi + phi0;
        if (mGenerationPhi >= 2.*M_PI) mGenerationPhi -= 2.*M_PI;

        // Check if the muon is inward
        if (sin(mGenerationTheta)*cos(mGenerationPhi)*mGenerationPosition[0] + sin(mGenerationTheta)*sin(mGenerationPhi)*mGenerationPosition[1] > 0) Generate();
      }
    }

    // Half-sphere generation
    if (mGenMethod == HSphere) {
      // Generation point on the half-sphere
      double phi0   = 2.*M_PI*mRandom.GenerateRandomDouble();
      double theta0 = 0.;
      double r0t    = 0.;
      while (!accepted) {
        r0t                 = mRandom.GenerateRandomDouble();
        r2                  = mRandom.GenerateRandomDouble();
        theta0              = acos(r0t);
        mGenerationTheta    = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationPhi      = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
        mGenerationMomentum = GenerateMomentumF1();
        n                   = 2.856-0.655*log(mGenerationMomentum);
        if (n < 0.1) n = 0.1;

        ftheta = 1600*pow(mGenerationMomentum, 0.279)*pow(cos(mGenerationTheta), n)*(sin(mGenerationTheta)*sin(theta0)*cos(mGenerationPhi)+cos(mGenerationTheta)*cos(theta0))*sin(mGenerationTheta);
        if (4891*r2 < ftheta) accepted = true;
      }

      mGenerationPosition[0] = mHSphereRadius*sin(theta0)*cos(phi0) + mHSphereCenterPosition[0];
      mGenerationPosition[1] = mHSphereRadius*sin(theta0)*sin(phi0) + mHSphereCenterPosition[1];
      mGenerationPosition[2] = mHSphereRadius*cos(theta0) + mHSphereCenterPosition[2];

      mGenerationTheta = M_PI - mGenerationTheta;
      mGenerationPhi = mGenerationPhi + phi0;
      if (mGenerationPhi >= 2*M_PI) mGenerationPhi -= 2*M_PI;

      mGenerationPhi += M_PI;
      if (mGenerationPhi >= 2*M_PI) mGenerationPhi -= 2*M_PI;
    }

    // Generate the charge
    mCharge = (mDiscDistC(mEngineC) == 0) ? 1 : -1;
  };
  ///////////////////////////////////////////////////////////////

};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#endif
