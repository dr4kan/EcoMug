/////////////////////////////////////////////////////////////////////////////////////
// EcoMug: Efficient COsmic MUon Generator                                         //
// Copyright (C) 2023 Davide Pagano <davide.pagano@unibs.it>                       //
//                                                                                 //
// EcoMug is based on the following work:                                          //
// D. Pagano, G. Bonomi, A. Donzella, A. Zenoni, G. Zumerle, N. Zurlo,             //
// "EcoMug: an Efficient COsmic MUon Generator for cosmic-ray muons applications", //
// doi:10.1016/j.nima.2021.165732                                                  //
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

#include <cmath>
#include <array>
#include <random>
#include <functional>
#include <iostream>
#include <initializer_list>
#include <sstream>  

#define ECOMUG_VERSION "2.1"

#ifndef M_PI
# define M_PI_NOT_DEFINED
# define M_PI 3.14159265358979323846
#endif

namespace EMUnits {
  // Default units:
  // meter              (m)
  // second             (s)
  // Giga electron Volt (GeV)
  // radian             (rad)

  // Lengths and areas
  static const double m     = 1.; 
  static const double cm    = 1.e-2*m;
  static const double mm    = 1.e-3*m;
  static const double km    = 1000.*m;                 
  static const double mm2   = mm*mm;
  static const double cm2   = cm*cm;
  static const double m2    = m*m;
  static const double km2   = km*km;    

  // Angles 
  static const double rad   = 1.;                  
  static const double mrad  = 1.e-3*rad;
  static const double deg   = (M_PI/180.0)*rad;

  // Time 
  static const double s     = 1.;
  static const double ms    = 1.e-3*s;
  static const double us    = 1.e-6*s;
  static const double ns    = 1.e-9*s;
  static const double min   = 60.*s;
  static const double hour  = 60.*min;
  static const double day   = 24.*hour;
  static const double hertz = 1./s;

  // Energy/momentum
  static const double GeV = 1.;
  static const double MeV = 1.e-3*GeV;
  static const double keV = 1.e-3*MeV;
  static const double TeV = 1.e+6*MeV;
  static const double  eV = 1.e-6*MeV;
};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


//! Class for the logging system
class EMLog {
public:
  enum TLogLevel {ERROR, WARNING, INFO, DEBUG};

  enum TMsgType {EcoMug, EMRandom, EMMultiGen, EMMaximization, UNKNOWN};

  EMLog() {};
  virtual ~EMLog() {
    os << std::endl;
    fprintf(stderr, "%s", os.str().c_str());
    fflush(stderr);
  };

  std::ostringstream& Get(TLogLevel level = INFO);
  void Get(TLogLevel level, const std::string& msg, TMsgType type = EcoMug) {
    std::cout << "[EcoMug v" << ECOMUG_VERSION << "] [";
    std::cout << ToString(level);
    if (type != UNKNOWN) std::cout << " in " << ToString(type);
    std::cout << "]: ";
    std::cout << std::string(level > DEBUG ? level - DEBUG : 0, '\t');
    std::cout << msg << "\n";
  };

public:
  static std::string ToString(TLogLevel level) {
    static const char* const buffer[] = {"ERROR", "WARNING", "INFO", "DEBUG"};
    return buffer[level];
  };

  static std::string ToString(TMsgType type) {
    static const char* const buffer[] = {"EcoMug", "EMRandom", "EMMultiGen", "EMMaximization", "UNKNOWN"};
    return buffer[type];
  };

  static TLogLevel FromString(const std::string& level) {
    if (level == "DEBUG")
    return DEBUG;
    if (level == "INFO")
    return INFO;
    if (level == "WARNING")
    return WARNING;
    if (level == "ERROR")
    return ERROR;
    EMLog().Get(WARNING) << "Unknown logging level '" << level << "'. Using INFO level as default.";
    return INFO;
  };

  static EMLog::TLogLevel ReportingLevel;
private:
  EMLog(const EMLog&);
  EMLog& operator =(const EMLog&);
  std::ostringstream os;
};

#define EMLogger(level, msg, type) \
if (level > EMLog::ReportingLevel) ; \
else EMLog().Get(level, msg, type)

//EMLog::TLogLevel EMLog::ReportingLevel = WARNING;
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


//! Fast generation of random numbers
//! This class is based on the xoroshiro128+ generator.
//! https://prng.di.unimi.it/
class EMRandom {
public:
  EMRandom() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::uint64_t> dis(0, std::numeric_limits<std::uint64_t>::max());
    s[0] = dis(gen);
    s[1] = dis(gen);
  };

  void SetSeed(std::uint64_t seed) {
    s[0] = seed;
    s[1] = seed;
  };

  double GenerateRandomDouble() {
    std::uint64_t x = next();
    return to_double(x);
  };

  double GenerateRandomDouble(double x1, double x2) {
    return (x2-x1)*GenerateRandomDouble()+x1;
  };

  int64_t rotl(const std::uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
  };

  std::uint64_t next() {
    const std::uint64_t s0 = s[0];
    std::uint64_t s1 = s[1];
    const std::uint64_t result = s0 + s1;
    s1 ^= s0;
    s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14);
    s[1] = rotl(s1, 36);
    return result;
  };

  double to_double(std::uint64_t x) {
    union U {
      std::uint64_t i;
      double d;
    };
    U u = { UINT64_C(0x3FF) << 52 | x >> 12 };
    return u.d - 1.0;
  };

  std::uint64_t s[2];
};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


//! Class for maximization based on "Whale Optimization Algorithm"
//! doi:10.1016/j.advengsoft.2016.01.008
class EMMaximization {
private:
  EMRandom mRandom;
  std::size_t mPopSize;
  std::size_t mNIter;
  int    mGenMethod; // 0 = sky, 1 = cylinder, 2 = hspere
  double m_a;
  double m_a2;
  std::vector<std::vector<double> > mRanges;
  std::vector<std::vector<double> > mPopulation;
  double mBestCost;
  std::vector<double> mBestSolution;
  std::function<double(double, double)> mFunc;

public:
  EMMaximization(const EMRandom& random, int genMethod) : mRandom(random), mPopSize(200),
  mNIter(500), mGenMethod(genMethod), m_a(0.), m_a2(0.), mBestCost(-1.) {
    mFunc = &DefaultJ;
  };
  ///////////////////////////////////////////////////////////////

  static double DefaultJ(double p, double theta) {
    double n = std::max(0.1, 2.856-0.655*log(p));
    return 1600*pow(p, 0.279)*pow(cos(theta), n);
  }
  ///////////////////////////////////////////////////////////////

  void SetParameters(double minP, double maxP, double minTheta, double maxTheta) {
    mRanges.push_back({minP, maxP});
    mRanges.push_back({minTheta, maxTheta});
  }
  ///////////////////////////////////////////////////////////////

  void SetParameters(double minP, double maxP, double minTheta, double maxTheta, double minPhi, double maxPhi) {
    mRanges.push_back({minP, maxP});
    mRanges.push_back({minTheta, maxTheta});
    mRanges.push_back({minPhi, maxPhi});
    mRanges.push_back({0, M_PI/2.});
  }
  ///////////////////////////////////////////////////////////////

  void SetFunction(std::function<double(double, double)> func) {
    mFunc = func;
  }
  ///////////////////////////////////////////////////////////////

  double SkyFunc(double p, double theta) {
    return mFunc(p, theta)*cos(theta)*sin(theta);
  }
  ///////////////////////////////////////////////////////////////

  double CylFunc(double p, double theta) {
    return mFunc(p, theta)*pow(sin(theta), 2);
  }
  ///////////////////////////////////////////////////////////////

  double HSFunc(double p, double theta, double phi, double theta0) {
    return mFunc(p, theta)*(sin(theta0)*sin(theta)*cos(phi) + cos(theta0)*cos(theta))*sin(theta);
  }
  ///////////////////////////////////////////////////////////////

  double Evaluate(std::vector<double> &v) {
    if (mGenMethod == 0) {
      return SkyFunc(v[0], v[1]);
    } else if (mGenMethod == 1) {
      return CylFunc(v[0], v[1]);
    } else {
      return HSFunc(v[0], v[1], v[2], v[3]);
    }
    return -1;
  }
  ///////////////////////////////////////////////////////////////

  void Evaluate() {
    double value;
    for (std::size_t i = 0; i < mPopSize; ++i) {
      value = Evaluate(mPopulation[i]);
      if (value > mBestCost) {
        mBestCost = value;
        mBestSolution = mPopulation[i];
      }
    }
  }
  ///////////////////////////////////////////////////////////////

  void Init() {
    std::size_t dim = mRanges.size();
    mPopulation.resize(mPopSize);
    for (std::size_t i = 0; i < mPopSize; ++i) {
      mPopulation[i].resize(dim);
      for (std::size_t j = 0; j < dim; ++j) {
        mPopulation[i][j] = mRandom.GenerateRandomDouble(mRanges[j][0], mRanges[j][1]);
      }
    }
  }
  ///////////////////////////////////////////////////////////////

  void UpdateParameters(std::size_t t) {
    m_a  = 2. - t*(2./mNIter);
    m_a2 = -1. + t*((-1.)/mNIter);
  }
  ///////////////////////////////////////////////////////////////

  void Move() {
    double r1, r2, A, C, b, l, rw, p, D_tmp, D_best, distance;
    std::vector<double> tmp;
    for (std::size_t i = 0; i < mPopulation.size(); ++i) {
      r1 = mRandom.GenerateRandomDouble();
      r2 = mRandom.GenerateRandomDouble();
      A  = 2*m_a*r1-m_a;
      C  = 2*r2;
      b  = 1.;
      l  = (m_a2-1)*mRandom.GenerateRandomDouble()+1;
      p  = mRandom.GenerateRandomDouble();

      for (std::size_t j = 0; j < mPopulation[0].size(); ++j) {
        if (p < 0.5) {
          if (fabs(A) >= 1) {
            rw = floor(mRandom.GenerateRandomDouble()*mPopulation.size());
            tmp = mPopulation[rw];
            D_tmp = fabs(C*tmp[j] - mPopulation[i][j]);
            mPopulation[i][j] = tmp[j] - A*D_tmp;
          } else {
            D_best = fabs(C*mBestSolution[j] - mPopulation[i][j]);
            mPopulation[i][j] = mBestSolution[j]-A*D_best;
          }
        } else {
          distance = fabs(mBestSolution[j] - mPopulation[i][j]);
          mPopulation[i][j] = distance*exp(b*l)*cos(l*2*M_PI) + mBestSolution[j];
        }
        if (mPopulation[i][j] < mRanges[j][0]) mPopulation[i][j] = mRanges[j][0];
        if (mPopulation[i][j] > mRanges[j][1]) mPopulation[i][j] = mRanges[j][1];
      }
    }
  }
  ///////////////////////////////////////////////////////////////

  double Maximize() {
    Init();
    Evaluate();
    for (std::size_t iter = 1; iter < mNIter; ++iter) {
      UpdateParameters(iter);
      Move();
      Evaluate();
    }
    return mBestCost;
  }
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
  double mHorizontalRate;
  double mCylinderMinPositionPhi;
  double mCylinderMaxPositionPhi;
  double mHSphereMinPositionPhi;
  double mHSphereMaxPositionPhi;
  double mHSphereMinPositionTheta;
  double mHSphereMaxPositionTheta;
  double mHSphereCosMinPositionTheta;
  double mHSphereCosMaxPositionTheta;
  double mJPrime;
  double mN;
  double mRandAccRej;
  double mPhi0;
  double mTheta0;
  bool   mAccepted;
  std::array<double, 2> mSkySize;
  std::array<double, 3> mSkyCenterPosition;
  double mCylinderHeight;
  double mCylinderRadius;
  std::array<double, 3> mCylinderCenterPosition;
  double mHSphereRadius;
  double mMaxFuncSkyCylinder;
  std::array<double, 3> mHSphereCenterPosition;
  bool mCustomJ;
  EMRandom mRandom;
  std::default_random_engine mEngineC;
  std::discrete_distribution<int> mDiscDistC;
  std::array<double, 3> mMaxJ;
  std::array<double, 3> mMaxCustomJ;
  std::function<double(double, double)> mJ;

public:
  // Default constructor
  EcoMug() : mGenMethod(Sky),
  mGenerationPosition({{0., 0., 0.}}), mGenerationTheta(0.), mGenerationPhi(0.),
  mGenerationMomentum(0.), mMinimumMomentum(0.01), mMaximumMomentum(1000.),
  mMinimumTheta(0.), mMaximumTheta(M_PI/2.), mMinimumPhi(0.), mMaximumPhi(2.*M_PI),
  mCharge(1), mHorizontalRate(129*EMUnits::hertz/EMUnits::m2), mCylinderMinPositionPhi(0.), mCylinderMaxPositionPhi(2.*M_PI),
  mHSphereMinPositionPhi(0.), mHSphereMaxPositionPhi(2.*M_PI), mHSphereMinPositionTheta(0.),
  mHSphereMaxPositionTheta(M_PI/2.), mHSphereCosMinPositionTheta(1.), mHSphereCosMaxPositionTheta(0.),
  mJPrime(0.), mN(0.), mRandAccRej(0.), mPhi0(0.), mTheta0(0.), mAccepted(false),
  mSkySize({{0., 0.}}), mSkyCenterPosition({{0., 0., 0.}}), mCylinderHeight(0.),
  mCylinderRadius(0.), mCylinderCenterPosition({{0., 0., 0.}}), mHSphereRadius(0.),
  mMaxFuncSkyCylinder(5.3176), mHSphereCenterPosition({{0., 0., 0.}}), mCustomJ(false),
  mEngineC(std::random_device{}()) {
    mDiscDistC = std::discrete_distribution<int>({128, 100});
    mMaxJ = {-1., -1., -1.};
    mMaxCustomJ = {-1., -1., -1.};
  };

    // Copy constructor
    EcoMug(const EcoMug& t) {
    mGenMethod = t.mGenMethod;
    mGenerationPosition = t.mGenerationPosition;
    mGenerationTheta = t.mGenerationTheta;
    mGenerationPhi = t.mGenerationPhi;
    mGenerationMomentum = t.mGenerationMomentum;
    mMinimumMomentum = t.mMinimumMomentum;
    mMaximumMomentum = t.mMaximumMomentum;
    mMinimumTheta = t.mMinimumTheta;
    mMaximumTheta = t.mMaximumTheta;
    mMinimumPhi = t.mMinimumPhi;
    mMaximumPhi = t.mMaximumPhi;
    mCharge = t.mCharge;
    mHorizontalRate = t.mHorizontalRate;
    mCylinderMinPositionPhi = t.mCylinderMinPositionPhi;
    mCylinderMaxPositionPhi = t.mCylinderMaxPositionPhi;
    mHSphereMinPositionPhi = t.mHSphereMinPositionPhi;
    mHSphereMaxPositionPhi = t.mHSphereMaxPositionPhi;
    mHSphereMinPositionTheta = t.mHSphereMinPositionTheta;
    mHSphereMaxPositionTheta = t.mHSphereMaxPositionTheta;
    mHSphereCosMinPositionTheta = t.mHSphereCosMinPositionTheta;
    mHSphereCosMaxPositionTheta = t.mHSphereCosMaxPositionTheta;
    mJPrime = t.mJPrime;
    mN = t.mN;
    mRandAccRej = t.mRandAccRej;
    mPhi0 = t.mPhi0;
    mTheta0 = t.mTheta0;
    mAccepted = t.mAccepted;
    mSkySize = t.mSkySize;
    mSkyCenterPosition = t.mSkyCenterPosition;
    mCylinderHeight = t.mCylinderHeight;
    mCylinderRadius = t.mCylinderRadius;
    mCylinderCenterPosition = t.mCylinderCenterPosition;
    mHSphereRadius = t.mHSphereRadius;
    mMaxFuncSkyCylinder = t.mMaxFuncSkyCylinder;
    mHSphereCenterPosition = t.mHSphereCenterPosition;
    mCustomJ = t.mCustomJ;
    mRandom = t.mRandom;
    mEngineC = t.mEngineC;
    mDiscDistC = t.mDiscDistC;
    mMaxJ = t.mMaxJ;
    mMaxCustomJ = t.mMaxCustomJ;
    mJ = t.mJ;
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
  /// Get the generation momentum
  void GetGenerationMomentum(std::array<double, 3>& momentum) const {
    momentum  = {
      mGenerationMomentum*sin(mGenerationTheta)*cos(mGenerationPhi),
      mGenerationMomentum*sin(mGenerationTheta)*sin(mGenerationPhi),
      mGenerationMomentum*cos(mGenerationTheta)
    };
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
  /// Set the differential flux J. Accepted functions are like
  /// double J(double momentum, double theta)
  /// momentum has to be in GeV/c and theta in radians
  void SetDifferentialFlux(std::function<double(double, double)> J) {
    mJ = J;
    mCustomJ = true;
  };
  /// Set the seed for the internal PRNG (if 0 a random seed is used)
  void SetSeed(std::uint64_t seed) {
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
  /// Set the rate of cosmic ray muons per square unit are through
  /// a horizontal surface. Default value is 129 Hz/m^2.
  void SetHorizontalRate(double rate) {
    mHorizontalRate = rate;
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
  /// Get horizontal rate of cosmic muons
  double GetHorizontalRate() const {
    return mHorizontalRate;
  };
  /// Get the generation surface area
  double GetGenSurfaceArea() const {
    double area = 0.;
    if (mGenMethod == Sky) {
      area = mSkySize[0]*mSkySize[1];
    } else if (mGenMethod == Cylinder) {
      // A = \Delta\phi r h
      area = (mCylinderMaxPositionPhi-mCylinderMinPositionPhi)*mCylinderRadius*mCylinderHeight;
    } else {
      // A = \Delta\phi r^2\left(\cos\theta_{min} - \cos\theta_{max}\right)
      area = (mHSphereMaxPositionPhi-mHSphereMinPositionPhi)*pow(mHSphereRadius, 2)*(cos(mHSphereMinPositionTheta) - cos(mHSphereMaxPositionTheta));
    }
    return area;
  };
  /// Get the average rate of cosmic muons as well the error on it
  /// The optional parameter npoints defines the number
  /// of points to be used in the MC integration of the J'
  void GetAverageGenRateAndError(double &rate, double &error, int npoints = 1e7) {
    // 129.0827 is the integral of the J'prime for the sky in
    // the full range of theta, phi and up to 3 TeV in energy.
    // For custom J it is the user who should account for the correction.
    double k = mHorizontalRate/129.0827; 
    if (mGenMethod == Sky) {
      if (mCustomJ) {
        MCJprimeCustomSkyIntegration(rate, error, npoints);
        return;
      } else MCJprimeSkyIntegration(rate, error, npoints);
    } else if (mGenMethod == Cylinder) {
      if (mCustomJ) {
        MCJprimeCustomCylinderIntegration(rate, error, npoints);
        return;
      } else MCJprimeCylinderIntegration(rate, error, npoints);
    } else {
      if (mCustomJ) {
        MCJprimeCustomHSphereIntegration(rate, error, npoints);
        return;
      } else MCJprimeHSphereIntegration(rate, error, npoints);
    }
    rate *= k;
    error *= k;
  };
  /// Get the average rate of cosmic muons 
  /// The optional parameter npoints defines the number
  /// of points to be used in the MC integration of the J'
  double GetAverageGenRate(int npoints = 1e7) {
    double rate, error;
    GetAverageGenRateAndError(rate, error, npoints);
    return rate;
  };
  /// Get the estimated corresponding to the provided statistics
  double GetEstimatedTime(int nmuons) {
    if (mCustomJ) return 0.;
    return (nmuons/(GetGenSurfaceArea()/EMUnits::m2))/(GetAverageGenRate()/EMUnits::hertz*EMUnits::m2);
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
  double GetSkySize(unsigned int index) {
    return mSkySize[index];
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
  void SetCylinderMinPositionPhi(double phi) {
    mCylinderMinPositionPhi = phi;
  };
  void SetCylinderMaxPositionPhi(double phi) {
    mCylinderMaxPositionPhi = phi;
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
  void SetHSphereMinPositionPhi(double phi) {
    mHSphereMinPositionPhi = phi;
  };
  void SetHSphereMaxPositionPhi(double phi) {
    mHSphereMaxPositionPhi = phi;
  };
  void SetHSphereMinPositionTheta(double theta) {
    mHSphereMinPositionTheta = theta;
    mHSphereCosMinPositionTheta = cos(mHSphereMinPositionTheta);
  };
  void SetHSphereMaxPositionTheta(double theta) {
    mHSphereMaxPositionTheta = theta;
    mHSphereCosMaxPositionTheta = cos(mHSphereMaxPositionTheta);
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

  double maxSkyJFunc() {
    return 1600*pow(mMaximumMomentum, 0.279)*pow(cos(0.76158), 1.1)*sin(0.76158);
  };

  double maxCylJFunc() {
    return 1600*pow(mMaximumMomentum, 0.279)*pow(cos(1.35081), 0.1)*pow(sin(1.35081), 2);
  };

  double maxHSJFunc() {
    return 1600*pow(mMaximumMomentum, 0.279)*pow(cos(1.26452), 0.1)*(sin(1.26452)*sin(1.26452)+cos(1.26452)*cos(1.26452))*sin(1.26452);
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

  void GeneratePositionCylinder() {
    mPhi0                  = mRandom.GenerateRandomDouble(mCylinderMinPositionPhi, mCylinderMaxPositionPhi);
    mGenerationPosition[0] = mCylinderCenterPosition[0] + mCylinderRadius*cos(mPhi0);
    mGenerationPosition[1] = mCylinderCenterPosition[1] + mCylinderRadius*sin(mPhi0);
    mGenerationPosition[2] = mRandom.GenerateRandomDouble(mCylinderCenterPosition[2]-mCylinderHeight/2., mCylinderCenterPosition[2]+mCylinderHeight/2.);
  };

  void ComputeMaximumCustomJ() {
    EMMaximization maximizer(mRandom, mGenMethod);
    maximizer.SetFunction(mJ);
    if (mGenMethod == 0 || mGenMethod == 1) {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta);
    } else {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta, mMinimumPhi, mMaximumPhi);
    }
    mMaxCustomJ[mGenMethod] = maximizer.Maximize();
  };

  void ComputeMaximum() {
    EMMaximization maximizer(mRandom, mGenMethod);
    if (mGenMethod == 0 || mGenMethod == 1) {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta);
    } else {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta, mMinimumPhi, mMaximumPhi);
    }
    mMaxJ[mGenMethod] = maximizer.Maximize();
  };

  void MCJprimeCustomSkyIntegration(double &rate, double &error, int npoints) { 
    double I = 0., I2 = 0., value = 0.;
    for (auto i = 0; i < npoints; ++i) {
      mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      value = mJ(mGenerationMomentum, mGenerationTheta)*cos(mGenerationTheta)*sin(mGenerationTheta);
      I  += value;
      I2 += pow(value, 2);
    }
    double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    double expected = I/npoints;
    double expectedSquare = I2/npoints;
    rate = V*I/npoints;
    error = V*pow((expectedSquare-pow(expected,2))/(npoints-1), 0.5);
  };

  void MCJprimeCustomCylinderIntegration(double &rate, double &error, int npoints) { 
    double I = 0., I2 = 0., value = 0.;
    for (auto i = 0; i < npoints; ++i) {
      mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      value = mJ(mGenerationMomentum, mGenerationTheta)*pow(sin(mGenerationTheta), 2)*cos(mGenerationPhi);
      if (value < 0) value = 0;
      I  += value;
      I2 += pow(value, 2);
    }
    double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    double expected = I/npoints;
    double expectedSquare = I2/npoints;
    rate = V*I/npoints;
    error = V*pow((expectedSquare-pow(expected,2))/(npoints-1), 0.5);
  };

    void MCJprimeCustomHSphereIntegration(double &rate, double &error, int npoints) { 
    double I = 0., I2 = 0., value = 0.;
    for (auto i = 0; i < npoints; ++i) {
      mTheta0 = mRandom.GenerateRandomDouble(mHSphereMinPositionTheta, mHSphereMaxPositionTheta);
      mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      value = mJ(mGenerationMomentum, mGenerationTheta)*(sin(mTheta0)*sin(mGenerationTheta)*sin(mGenerationTheta)*cos(mGenerationPhi)+cos(mTheta0)*cos(mGenerationTheta)*sin(mGenerationTheta))*sin(mTheta0);
      if (value < 0) value = 0;
      I  += value;
      I2 += pow(value, 2);
    }
    double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi)*(mHSphereMaxPositionTheta-mHSphereMinPositionTheta);
    double expected = I/npoints;
    double expectedSquare = I2/npoints;
    rate = V*I/npoints;
    error = V*pow((expectedSquare-pow(expected,2))/(npoints-1), 0.5);
  };

  void MCJprimeSkyIntegration(double &rate, double &error, int npoints) { 
    double I = 0., I2 = 0., value = 0.;
    for (auto i = 0; i < npoints; ++i) {
      mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      mN = 2.856-0.655*log(mGenerationMomentum);
      if (mN < 0.1) mN = 0.1;
      value = 1600*pow(mGenerationMomentum+2.68, -3.175)*pow(mGenerationMomentum, 0.279)*pow(cos(mGenerationTheta), mN+1)*sin(mGenerationTheta);
      I  += value;
      I2 += pow(value, 2);
    }
    double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    double expected = I/npoints;
    double expectedSquare = I2/npoints;
    rate = V*I/npoints;
    error = V*pow((expectedSquare-pow(expected,2))/(npoints-1), 0.5);
  };

  void MCJprimeCylinderIntegration(double &rate, double &error, int npoints) { 
    double I = 0., I2 = 0., value = 0.;
    for (auto i = 0; i < npoints; ++i) {
      mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      mN = 2.856-0.655*log(mGenerationMomentum);
      if (mN < 0.1) mN = 0.1;
      value = 1600*pow(mGenerationMomentum+2.68, -3.175)*pow(mGenerationMomentum, 0.279)*pow(cos(mGenerationTheta), mN)*pow(sin(mGenerationTheta), 2)*cos(mGenerationPhi);
      if (value < 0) value = 0;
      I  += value;
      I2 += pow(value, 2);
    }
    double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    double expected = I/npoints;
    double expectedSquare = I2/npoints;
    rate = V*I/npoints;
    error = V*pow((expectedSquare-pow(expected,2))/(npoints-1), 0.5);
  };

  void MCJprimeHSphereIntegration(double &rate, double &error, int npoints) { 
    double I = 0., I2 = 0., value = 0.;
    for (auto i = 0; i < npoints; ++i) {
      mTheta0 = mRandom.GenerateRandomDouble(mHSphereMinPositionTheta, mHSphereMaxPositionTheta);
      mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      mN = 2.856-0.655*log(mGenerationMomentum);
      if (mN < 0.1) mN = 0.1;
      value = 1600*pow(mGenerationMomentum+2.68, -3.175)*pow(mGenerationMomentum, 0.279)*pow(cos(mGenerationTheta), mN)*(sin(mTheta0)*sin(mGenerationTheta)*sin(mGenerationTheta)*cos(mGenerationPhi)+cos(mTheta0)*cos(mGenerationTheta)*sin(mGenerationTheta))*sin(mTheta0);
      if (value < 0) value = 0;
      I  += value;
      I2 += pow(value, 2);
    }
    double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi)*(mHSphereMaxPositionTheta-mHSphereMinPositionTheta);
    double expected = I/npoints;
    double expectedSquare = I2/npoints;
    rate = V*I/npoints;
    error = V*pow((expectedSquare-pow(expected,2))/(npoints-1), 0.5);
  };

public:
  ///////////////////////////////////////////////////////////////
  /// Generate a cosmic muon from the pre-defined J
  ///////////////////////////////////////////////////////////////
  void Generate() {
    mAccepted = false;

    if (mMaxJ[mGenMethod] < 0) ComputeMaximum();

    // Sky or cylinder generation
    if (mGenMethod == Sky || mGenMethod == Cylinder) {
      // Generation of the momentum and theta angle
      while (!mAccepted) {
        mRandAccRej  = mRandom.GenerateRandomDouble();
        mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationMomentum = GenerateMomentumF1();
        mN = 2.856-0.655*log(mGenerationMomentum);
        if (mN < 0.1) mN = 0.1;

        if (mGenMethod == Sky) {
          mJPrime = 1600*pow(mGenerationMomentum, 0.279)*pow(cos(mGenerationTheta), mN+1)*sin(mGenerationTheta);
          if (mMaxJ[mGenMethod]*mRandAccRej < mJPrime) mAccepted = true;
        }

        if(mGenMethod == Cylinder)  {
          mJPrime = 1600*pow(mGenerationMomentum, 0.279)*pow(cos(mGenerationTheta), mN)*pow(sin(mGenerationTheta), 2);
          if (mMaxJ[mGenMethod]*mRandAccRej < mJPrime) mAccepted = true;
        }
      }
      mGenerationTheta = M_PI - mGenerationTheta;

      // Generation of the position and phi angle
      if (mGenMethod == Sky) {
        GeneratePositionSky();
        mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      }
      if (mGenMethod == Cylinder) {
        mAccepted = false;
        GeneratePositionCylinder();
        while (!mAccepted) {
          mRandAccRej  = mRandom.GenerateRandomDouble();
          mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
          if (mRandAccRej < fabs(cos(mGenerationPhi))) mAccepted = true;
        }
        mGenerationPhi = mGenerationPhi + mPhi0;
        if (mGenerationPhi >= 2.*M_PI) mGenerationPhi -= 2.*M_PI;

        // Check if the muon is inward
        if (sin(mGenerationTheta)*cos(mGenerationPhi)*mGenerationPosition[0] + sin(mGenerationTheta)*sin(mGenerationPhi)*mGenerationPosition[1] > 0) Generate();
      }
    }

    // Half-sphere generation
    if (mGenMethod == HSphere) {
      // Generation point on the half-sphere
      mPhi0      = mRandom.GenerateRandomDouble(mHSphereMinPositionPhi, mHSphereMaxPositionPhi);
      while (!mAccepted) {
        mRandAccRej         = mRandom.GenerateRandomDouble();
        mTheta0             = acos(mRandom.GenerateRandomDouble(mHSphereCosMaxPositionTheta, mHSphereCosMinPositionTheta));
        mGenerationTheta    = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationPhi      = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
        mGenerationMomentum = GenerateMomentumF1();
        mN                  = 2.856-0.655*log(mGenerationMomentum);
        if (mN < 0.1) mN = 0.1;

        mJPrime = 1600*pow(mGenerationMomentum, 0.279)*pow(cos(mGenerationTheta), mN)*(sin(mGenerationTheta)*sin(mTheta0)*cos(mGenerationPhi)+cos(mGenerationTheta)*cos(mTheta0))*sin(mGenerationTheta);
        if (mJPrime > 0 && mMaxJ[mGenMethod]*mRandAccRej < mJPrime) mAccepted = true;
      }

      mGenerationPosition[0] = mHSphereRadius*sin(mTheta0)*cos(mPhi0) + mHSphereCenterPosition[0];
      mGenerationPosition[1] = mHSphereRadius*sin(mTheta0)*sin(mPhi0) + mHSphereCenterPosition[1];
      mGenerationPosition[2] = mHSphereRadius*cos(mTheta0) + mHSphereCenterPosition[2];

      mGenerationTheta = M_PI - mGenerationTheta;
      mGenerationPhi = mGenerationPhi + mPhi0;
      if (mGenerationPhi >= 2*M_PI) mGenerationPhi -= 2*M_PI;

      mGenerationPhi += M_PI;
      if (mGenerationPhi >= 2*M_PI) mGenerationPhi -= 2*M_PI;
    }

    // Generate the charge
    mCharge = (mDiscDistC(mEngineC) == 0) ? 1 : -1;
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  /// Generate a cosmic muon for the user-defined J
  ///////////////////////////////////////////////////////////////
  void GenerateFromCustomJ() {
    mAccepted = false;

    std::cout << mMaxCustomJ[mGenMethod] << std::endl;

    if (mMaxCustomJ[mGenMethod] < 0) ComputeMaximumCustomJ();

    // Sky or cylinder generation
    if (mGenMethod == Sky || mGenMethod == Cylinder) {
      // Generation of the momentum and theta angle
      while (!mAccepted) {
        mRandAccRej  = mRandom.GenerateRandomDouble();
        mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);

        if (mGenMethod == Sky) {
          mJPrime = mJ(mGenerationMomentum, mGenerationTheta)*cos(mGenerationTheta)*sin(mGenerationTheta);
          if (mMaxCustomJ[mGenMethod]*mRandAccRej < mJPrime) mAccepted = true;
        }

        if(mGenMethod == Cylinder)  {
          mJPrime = mJ(mGenerationMomentum, mGenerationTheta)*pow(sin(mGenerationTheta), 2)*cos(mGenerationPhi);
          if (mMaxCustomJ[mGenMethod]*mRandAccRej < mJPrime) mAccepted = true;
        }
      }
      mGenerationTheta = M_PI - mGenerationTheta;

      // Generation of the position and phi angle
      if (mGenMethod == Sky) {
        GeneratePositionSky();
        mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      }
      if (mGenMethod == Cylinder) {
        mAccepted = false;
        GeneratePositionCylinder();
        while (!mAccepted) {
          mRandAccRej  = mRandom.GenerateRandomDouble();
          mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
          if (mRandAccRej < fabs(cos(mGenerationPhi))) mAccepted = true;
        }
        mGenerationPhi = mGenerationPhi + mPhi0;
        if (mGenerationPhi >= 2.*M_PI) mGenerationPhi -= 2.*M_PI;

        // Check if the muon is inward
        if (sin(mGenerationTheta)*cos(mGenerationPhi)*mGenerationPosition[0] + sin(mGenerationTheta)*sin(mGenerationPhi)*mGenerationPosition[1] > 0) Generate();
      }
    }

    // Half-sphere generation
    if (mGenMethod == HSphere) {
      // Generation point on the half-sphere
      mPhi0                 = mRandom.GenerateRandomDouble(mHSphereMinPositionPhi, mHSphereMaxPositionPhi);
      while (!mAccepted) {
        mRandAccRej         = mRandom.GenerateRandomDouble();
        mTheta0             = acos(mRandom.GenerateRandomDouble(mHSphereCosMaxPositionTheta, mHSphereCosMinPositionTheta));
        mGenerationTheta    = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationPhi      = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
        mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);

        mJPrime = mJ(mGenerationMomentum, mGenerationTheta)*(sin(mTheta0)*sin(mGenerationTheta)*cos(mGenerationPhi) + cos(mTheta0)*cos(mGenerationTheta))*sin(mGenerationTheta);
        if (mMaxCustomJ[mGenMethod]*mRandAccRej < mJPrime) mAccepted = true;
      }

      mGenerationPosition[0] = mHSphereRadius*sin(mTheta0)*cos(mPhi0) + mHSphereCenterPosition[0];
      mGenerationPosition[1] = mHSphereRadius*sin(mTheta0)*sin(mPhi0) + mHSphereCenterPosition[1];
      mGenerationPosition[2] = mHSphereRadius*cos(mTheta0) + mHSphereCenterPosition[2];

      mGenerationTheta = M_PI - mGenerationTheta;
      mGenerationPhi = mGenerationPhi + mPhi0;
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


//! Class to handle the generation from multiple distributions,
//! for example to take into account backgound sources.
class EMMultiGen {
public:
EMMultiGen(const EcoMug& signal, const std::vector<EcoMug>& backgrounds) : 
  mIndex(-1), mSigInstance(signal), mBckInstances{backgrounds}, mLimits(backgrounds.size()+2), mWeights(backgrounds.size()+1, 1.), 
  mPID(backgrounds.size()+1), mRd(std::random_device{}()) {
  for (std::size_t i = 0; i < mLimits.size(); ++i) mLimits[i] = i;
  mDd = std::piecewise_constant_distribution<>(mLimits.begin(), mLimits.end(), mWeights.begin());
};

/// Set the weights for all EcoMug background instance. The number of elements
/// must be equal to the background instances defined
void SetBckWeights(const std::vector<double>& weights) {
  if (mBckInstances.size() != weights.size()) {
    EMLogger(EMLog::ERROR, "Expected " + std::to_string(mBckInstances.size()) + " weights, but " + std::to_string(weights.size()) + " were provided. Setting them to 1.", EMLog::EMMultiGen);
    std::fill(mWeights.begin(), mWeights.end(), 1);
  } else {
    for (std::size_t i = 0; i < weights.size(); ++i) mWeights[i+1] = weights[i];
  }
  mDd = std::piecewise_constant_distribution<>(mLimits.begin(), mLimits.end(), mWeights.begin());
};

/// Set the PID for all background instances. 
void SetBckPID(const std::vector<int>& values) {
  if (mBckInstances.size() != values.size()) {
    EMLogger(EMLog::ERROR, "Expected " + std::to_string(mBckInstances.size()) + " PID, but " + std::to_string(values.size()) + " were provided. Setting them to 0.", EMLog::EMMultiGen);
    std::fill(mPID.begin(), mPID.end(), 0);
    return;
  }
  for (std::size_t i = 0; i < values.size(); ++i) mPID[i+1] = values[i];
};

/// Get the generation position
const std::array<double, 3>& GetGenerationPosition() const {
  return mBckInstances[mIndex].GetGenerationPosition();
};

/// Get the generation momentum
double GetGenerationMomentum() const {
  return mBckInstances[mIndex].GetGenerationMomentum();
};

/// Get the generation momentum
void GetGenerationMomentum(std::array<double, 3>& momentum) const {
  mBckInstances[mIndex].GetGenerationMomentum(momentum);
};

/// Get the generation theta
double GetGenerationTheta() const {
  return mBckInstances[mIndex].GetGenerationTheta();
};

/// Get the generation phi
double GetGenerationPhi() const {
  return mBckInstances[mIndex].GetGenerationPhi();
};

/// Get PID
int GetPID() const {
  // muon case
  if (mPID[mIndex] == 0) {
    if (mSigInstance.GetCharge() < 0) return 13;
    else return -13;
  }
  return mPID[mIndex];
};

void Generate() {
  mIndex = (int) mDd(mRd);
  if (mIndex == 0) mSigInstance.Generate(); 
  else mBckInstances[mIndex-1].Generate();
};

private:
int mIndex;
EcoMug mSigInstance;
std::vector<EcoMug> mBckInstances;
std::vector<double> mLimits;
std::vector<double> mWeights;
std::vector<int> mPID;
std::default_random_engine mRd;
std::piecewise_constant_distribution<> mDd;
};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#ifdef ECOMUG_VERSION 
#undef ECOMUG_VERSION
#endif

#ifdef M_PI_NOT_DEFINED
#undef M_PI
#undef M_PI_NOT_DEFINED
#endif

#endif
