/////////////////////////////////////////////////////////////////////////////////////
// EcoMug: Efficient COsmic MUon Generator                                         //
// Copyright (C) 2022 Davide Pagano <davide.pagano@unibs.it>                       //
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
#include <cstdint>
#include <cstring>
#include <random>
#include <functional>
#include <iostream>
#include <initializer_list>
#include <sstream>  

#define ECOMUG_VERSION "3.0"

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

  inline static TLogLevel ReportingLevel = WARNING;
private:
  EMLog(const EMLog&);
  EMLog& operator =(const EMLog&);
  std::ostringstream os;
};

#define EMLogger(level, msg, type) \
if (level > EMLog::ReportingLevel) ; \
else EMLog().Get(level, msg, type)

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

  //! Scramble a seed value with the SplitMix64 generator, advancing state.
  //! The authors of xoroshiro128+ recommend initialising its state this way:
  //! a state with few bits set (as obtained by copying a small seed into both
  //! words) needs many iterations before the output looks random.
  static std::uint64_t SplitMix64(std::uint64_t& state) {
    std::uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30))*UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27))*UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
  };

  void SetSeed(std::uint64_t seed) {
    std::uint64_t state = seed;
    s[0] = SplitMix64(state);
    s[1] = SplitMix64(state);
  };

  double GenerateRandomDouble() {
    std::uint64_t x = next();
    return to_double(x);
  };

  double GenerateRandomDouble(double x1, double x2) {
    return (x2-x1)*GenerateRandomDouble()+x1;
  };

  std::uint64_t rotl(const std::uint64_t x, int k) {
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

  double to_double(std::uint64_t x) const {
    // Take the top 52 bits: the low bits of xoroshiro128+ are the weakest.
    // std::memcpy rather than a union, because type-punning through a union is
    // undefined behaviour in C++ (it is only legal in C). Both compile to a
    // single register move.
    const std::uint64_t bits = UINT64_C(0x3FF) << 52 | x >> 12;
    double d;
    std::memcpy(&d, &bits, sizeof d);
    return d - 1.0;
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
  int    mGenMethod; // mirrors EcoMug::EMGeometry: 0 = sky, 1 = cylinder,
                     // 2 = half-sphere, 3 = target sphere
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
    return mFunc(p, theta)*sin(theta)*sin(theta);
  }
  ///////////////////////////////////////////////////////////////

  double HSFunc(double p, double theta, double phi, double theta0) {
    return mFunc(p, theta)*(sin(theta0)*sin(theta)*cos(phi) + cos(theta0)*cos(theta))*sin(theta);
  }
  ///////////////////////////////////////////////////////////////

  /// Target-sphere mode. The generation disc is perpendicular to the muon, so
  /// there is no projection cosine here: only the solid-angle Jacobian.
  double TargetFunc(double p, double theta) {
    return mFunc(p, theta)*sin(theta);
  }
  ///////////////////////////////////////////////////////////////

  double Evaluate(std::vector<double> &v) {
    if (mGenMethod == 0) {
      return SkyFunc(v[0], v[1]);
    } else if (mGenMethod == 1) {
      return CylFunc(v[0], v[1]);
    } else if (mGenMethod == 3) {
      return TargetFunc(v[0], v[1]);
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

  // Parameters of the F1 momentum distribution used as the built-in momentum
  // spectrum (and, when it pays off, as the proposal for a user-supplied flux):
  //   CDF(p) = 1 - kF1Norm*(p + kF1Shift)^-kF1Index
  static constexpr double kF1Norm  = 8.534790171171021;
  static constexpr double kF1Shift = 2.68;
  static constexpr double kF1Index = 87./40.;

public:
  /// Possible generation methods
  enum EMGeometry {
    Sky,          ///< generation from a plane (flat sky)
    Cylinder,     ///< generation from a cylinder
    HSphere,      ///< generation from a half-sphere
    TargetSphere  ///< generation aimed at a sphere enclosing the detector
  };
  /// Number of generation geometries, i.e. the size of the per-geometry caches
  static const int kNGeometries = 4;

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
  double mTargetSphereRadius;
  std::array<double, 3> mTargetSphereCenterPosition;
  std::array<double, kNGeometries> mMaxJ;
  std::array<double, kNGeometries> mMaxCustomJ;
  std::array<bool, kNGeometries> mCustomJUseF1;
  double mF1Min;
  double mF1Max;
  std::function<double(double, double)> mJ;

  /// F1Cumulative at the momentum limits. Cached because GenerateMomentumF1 needs
  /// both on every single call while they only change when the limits change.
  void UpdateF1Range() {
    mF1Min = F1Cumulative(mMinimumMomentum);
    mF1Max = F1Cumulative(mMaximumMomentum);
  };

  /// Drop the cached accept-reject envelopes. Must be called by every setter that
  /// changes the domain the envelope was computed over, otherwise a stale (too
  /// small) envelope silently clips the generated distribution.
  void InvalidateMaxima() {
    mMaxJ.fill(-1.);
    mMaxCustomJ.fill(-1.);
  };

  /// Generate the muon charge, with a mu+/mu- ratio of 128/100.
  /// It is drawn from mRandom so that it follows the seed set by SetSeed and
  /// is reproducible across standard library implementations.
  void GenerateCharge() {
    constexpr double posFraction = 128./(128.+100.);
    mCharge = (mRandom.GenerateRandomDouble() < posFraction) ? 1 : -1;
  };

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
  mTargetSphereRadius(0.), mTargetSphereCenterPosition({{0., 0., 0.}}),
  mF1Min(0.), mF1Max(0.) {
    InvalidateMaxima();
    mCustomJUseF1.fill(false);
    UpdateF1Range();
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
    mTargetSphereRadius = t.mTargetSphereRadius;
    mTargetSphereCenterPosition = t.mTargetSphereCenterPosition;
    mCustomJ = t.mCustomJ;
    mRandom = t.mRandom;
    mMaxJ = t.mMaxJ;
    mMaxCustomJ = t.mMaxCustomJ;
    mCustomJUseF1 = t.mCustomJUseF1;
    mF1Min = t.mF1Min;
    mF1Max = t.mF1Max;
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
  /// True if a user-supplied differential flux was set with SetDifferentialFlux,
  /// i.e. if this instance should be driven with GenerateFromCustomJ()
  bool UsesCustomFlux() const {
    return mCustomJ;
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
  /// Aim the generation at a sphere enclosing the detector.
  ///
  /// This replaces the generation surface: instead of covering a plane, a
  /// cylinder or a dome and throwing away the muons that miss the apparatus,
  /// every muon is generated on the sphere of radius SetTargetSphereRadius()
  /// around SetTargetSphereCenterPosition(), already pointing through it.
  /// Typical detector setups keep fewer than one generated muon in a thousand,
  /// so this is usually the fastest way to feed a transport code.
  ///
  /// The construction is exact, not an approximation: for a given direction the
  /// muons crossing a sphere of radius R are exactly those crossing the disc of
  /// radius R perpendicular to that direction, so the muon is started uniformly
  /// on that disc and traced back onto the sphere. Because the disc is normal to
  /// the muon there is no projection cosine, and the generation area is pi*R^2
  /// independently of the direction. Rates and GetEstimatedTime() are normalised
  /// accordingly and stay comparable with the other geometries.
  void SetUseTargetSphere() {
    mGenMethod = TargetSphere;
  };
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
    InvalidateMaxima();
  };
  /// Set the seed for the internal PRNG (if 0 a random seed is used)
  void SetSeed(std::uint64_t seed) {
    if (seed > 0) mRandom.SetSeed(seed);
  };
  /// Set minimum generation Momentum
  void SetMinimumMomentum(double momentum) {
    mMinimumMomentum = momentum;
    UpdateF1Range();
    InvalidateMaxima();
  };
  /// Set maximum generation Momentum
  void SetMaximumMomentum(double momentum) {
    mMaximumMomentum = momentum;
    UpdateF1Range();
    InvalidateMaxima();
  };
  /// Set minimum generation Theta
  void SetMinimumTheta(double theta) {
    mMinimumTheta = theta;
    InvalidateMaxima();
  };
  /// Set maximum generation Theta
  void SetMaximumTheta(double theta) {
    mMaximumTheta = theta;
    InvalidateMaxima();
  };
  /// Set minimum generation Phi
  void SetMinimumPhi(double phi) {
    mMinimumPhi = phi;
    InvalidateMaxima();
  };
  /// Set maximum generation Phi
  void SetMaximumPhi(double phi) {
    mMaximumPhi = phi;
    InvalidateMaxima();
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
    if (mGenMethod == TargetSphere) {
      // the effective area is the disc the muons are generated on, the same for
      // every direction because the disc is perpendicular to the muon
      area = M_PI*mTargetSphereRadius*mTargetSphereRadius;
    } else if (mGenMethod == Sky) {
      area = mSkySize[0]*mSkySize[1];
    } else if (mGenMethod == Cylinder) {
      // A = \Delta\phi r h
      area = (mCylinderMaxPositionPhi-mCylinderMinPositionPhi)*mCylinderRadius*mCylinderHeight;
    } else {
      // A = \Delta\phi r^2\left(\cos\theta_{min} - \cos\theta_{max}\right)
      area = (mHSphereMaxPositionPhi-mHSphereMinPositionPhi)*mHSphereRadius*mHSphereRadius*(cos(mHSphereMinPositionTheta) - cos(mHSphereMaxPositionTheta));
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
    if (mGenMethod == TargetSphere) {
      if (mCustomJ) {
        MCJprimeCustomTargetIntegration(rate, error, npoints);
        return;
      } else MCJprimeTargetIntegration(rate, error, npoints);
      rate *= k;
      error *= k;
      return;
    }
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
    if (mCustomJ) {
      // A user-supplied flux carries its own normalisation, so EcoMug cannot turn
      // a number of muons into a live time. Say so rather than returning a silent
      // zero that the caller is likely to divide by.
      EMLogger(EMLog::WARNING, "GetEstimatedTime() returns 0 for a user-supplied differential flux: the normalisation of that flux is not known to EcoMug, so the live time has to be computed by the user.", EMLog::EcoMug);
      return 0.;
    }
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
  // Methods for the target-sphere generation
  ///////////////////////////////////////////////////////////////
  /// Set the radius of the sphere enclosing the detector
  void SetTargetSphereRadius(double radius) {
    if (radius <= 0.) {
      EMLogger(EMLog::ERROR, "The target sphere radius must be positive.", EMLog::EcoMug);
      return;
    }
    mTargetSphereRadius = radius;
  };
  /// Set the centre of the sphere enclosing the detector
  void SetTargetSphereCenterPosition(const std::array<double, 3>& position) {
    mTargetSphereCenterPosition = position;
  };
  /// Get the target sphere radius
  double GetTargetSphereRadius() const {
    return mTargetSphereRadius;
  };
  /// Get the target sphere centre
  const std::array<double, 3>& GetTargetSphereCenterPosition() const {
    return mTargetSphereCenterPosition;
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
  double F1Cumulative(double x) const {
    return 1. - kF1Norm/pow(x + kF1Shift, kF1Index);
  };

  double F1Inverse(double x) const {
    const double w = pow(1. - x, 1./kF1Index);
    return (kF1Shift - kF1Shift*w)/w;
  };

  /// Probability density of the momentum proposal produced by GenerateMomentumF1,
  /// normalised over [mMinimumMomentum, mMaximumMomentum]
  double F1ProposalDensity(double x) const {
    return kF1Index*kF1Norm/(pow(x + kF1Shift, kF1Index + 1.)*(mF1Max - mF1Min));
  };

  double maxSkyJFunc() {
    return 1600*pow(mMaximumMomentum, 0.279)*pow(cos(0.76158), 1.1)*sin(0.76158);
  };

  double maxCylJFunc() {
    return 1600*pow(mMaximumMomentum, 0.279)*pow(cos(1.35081), 0.1)*sin(1.35081)*sin(1.35081);
  };

  double maxHSJFunc() {
    return 1600*pow(mMaximumMomentum, 0.279)*pow(cos(1.26452), 0.1)*(sin(1.26452)*sin(1.26452)+cos(1.26452)*cos(1.26452))*sin(1.26452);
  };

  double GenerateMomentumF1() {
    return F1Inverse(mRandom.GenerateRandomDouble(mF1Min, mF1Max));
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

  /// Place the muon on the target sphere: uniformly over the disc of radius R
  /// perpendicular to its direction, then traced back onto the sphere surface so
  /// that it starts just outside the detector.
  ///
  /// mGenerationTheta and mGenerationPhi must already hold the final (downward)
  /// direction when this is called.
  void GeneratePositionTargetSphere() {
    const double sinTheta = sin(mGenerationTheta), cosTheta = cos(mGenerationTheta);
    const double sinPhi   = sin(mGenerationPhi),   cosPhi   = cos(mGenerationPhi);
    // direction of flight
    const double dx = sinTheta*cosPhi, dy = sinTheta*sinPhi, dz = cosTheta;
    // orthonormal basis of the plane perpendicular to it (spherical unit vectors)
    const double e1x = cosTheta*cosPhi, e1y = cosTheta*sinPhi, e1z = -sinTheta;
    const double e2x = -sinPhi,         e2y = cosPhi,          e2z = 0.;
    // uniform point in the disc of radius R
    const double r     = mTargetSphereRadius*sqrt(mRandom.GenerateRandomDouble());
    const double alpha = mRandom.GenerateRandomDouble(0., 2.*M_PI);
    const double u = r*cos(alpha), v = r*sin(alpha);
    // step back onto the sphere along -d
    const double back = sqrt(std::max(0., mTargetSphereRadius*mTargetSphereRadius - r*r));
    mGenerationPosition[0] = mTargetSphereCenterPosition[0] + u*e1x + v*e2x - back*dx;
    mGenerationPosition[1] = mTargetSphereCenterPosition[1] + u*e1y + v*e2y - back*dy;
    mGenerationPosition[2] = mTargetSphereCenterPosition[2] + u*e1z + v*e2z - back*dz;
  };

  /// True if the generated muon points away from the cylinder axis, in which
  /// case the event has to be resampled.
  bool IsOutwardCylinder() const {
    const double sinTheta = sin(mGenerationTheta);
    return sinTheta*cos(mGenerationPhi)*mGenerationPosition[0]
         + sinTheta*sin(mGenerationPhi)*mGenerationPosition[1] > 0;
  };

  void SetMaximizerRanges(EMMaximization& maximizer) const {
    if (mGenMethod == Sky || mGenMethod == Cylinder || mGenMethod == TargetSphere) {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta);
    } else {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta, mMinimumPhi, mMaximumPhi);
    }
  };

  /// Pick the momentum proposal for the user-supplied flux and compute the
  /// corresponding accept-reject envelope.
  ///
  /// Rejection sampling accepts with probability I/M, where M is the supremum of
  /// target/proposal-density over the sampled domain, so the cheapest proposal is
  /// the one with the smallest M. Two candidates are compared:
  ///   - uniform in [pmin, pmax]: what EcoMug always used. For a realistic, steeply
  ///     falling cosmic-ray flux this is very inefficient (of order 1000 trials per
  ///     accepted muon for p in [0.5, 200] GeV/c), because the target spans several
  ///     decades over the momentum range while the proposal is flat.
  ///   - the built-in F1 spectrum, i.e. the same inverse-CDF sampling Generate()
  ///     uses. Any flux with a roughly power-law momentum dependence is close to
  ///     F1, so the ratio is nearly flat and the acceptance is high.
  /// Both suprema come from the same maximiser, once per geometry, so a flux that
  /// does not look like F1 at all simply keeps the uniform proposal.
  void ComputeMaximumCustomJ() {
    EMMaximization uniformMaximizer(mRandom, mGenMethod);
    uniformMaximizer.SetFunction(mJ);
    SetMaximizerRanges(uniformMaximizer);
    const double supUniform = uniformMaximizer.Maximize();

    EMMaximization f1Maximizer(mRandom, mGenMethod);
    f1Maximizer.SetFunction([this](double p, double theta) {
      return mJ(p, theta)/F1ProposalDensity(p);
    });
    SetMaximizerRanges(f1Maximizer);
    const double supF1 = f1Maximizer.Maximize();

    // Express the uniform envelope per unit proposal density too, so that the two
    // are comparable: the uniform density is 1/(pmax - pmin).
    const double envUniform = supUniform*(mMaximumMomentum - mMinimumMomentum);

    mCustomJUseF1[mGenMethod] = (supF1 > 0. && supF1 < envUniform);
    mMaxCustomJ[mGenMethod]   = mCustomJUseF1[mGenMethod] ? supF1 : supUniform;
  };

  void ComputeMaximum() {
    EMMaximization maximizer(mRandom, mGenMethod);
    SetMaximizerRanges(maximizer);
    mMaxJ[mGenMethod] = maximizer.Maximize();
  };

  /// Draw the momentum for the custom-J accept-reject, and return the value the
  /// envelope has to be weighted with (1 for the uniform proposal).
  double GenerateMomentumCustomJ() {
    if (mCustomJUseF1[mGenMethod]) {
      mGenerationMomentum = GenerateMomentumF1();
      return F1ProposalDensity(mGenerationMomentum);
    }
    mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
    return 1.;
  };

  /// Accumulate a Monte Carlo estimate of the rate integral.
  ///
  /// Every integration below draws from a COPY of mRandom and keeps its samples in
  /// local variables. Previously they consumed the generator's own stream and wrote
  /// into mGenerationTheta / mGenerationMomentum / mTheta0 / ..., so a single
  /// GetAverageGenRate() call both overwrote the last generated muon and shifted
  /// every subsequent muon of a seeded run.
  struct MCAccumulator {
    double I = 0., I2 = 0.;
    void Add(double value) {
      if (value < 0.) value = 0.;
      I += value;
      I2 += value*value;
    };
    void Finish(double V, int npoints, double &rate, double &error) const {
      rate = V*I/npoints;
      if (npoints < 2) { error = 0.; return; }
      const double expected = I/npoints;
      const double variance = I2/npoints - expected*expected;
      error = V*sqrt(variance > 0. ? variance/(npoints-1) : 0.);
    };
  };

  /// The built-in flux is J(p,theta) = 1600*(p+2.68)^-3.175*p^0.279*cos^n(theta).
  /// Two different momentum factors are needed depending on how p was drawn, and
  /// mixing them up double-counts the (p+2.68)^-3.175 term:
  ///
  ///   BuiltinFluxComplete  -- the whole thing. Use where p is sampled UNIFORMLY,
  ///                           i.e. in the rate integrals.
  ///   BuiltinFluxResidual  -- only 1600*p^0.279. Use where p was drawn from F1,
  ///                           i.e. in the accept-reject of Generate(), because
  ///                           the F1 density already supplies (p+2.68)^-3.175.
  ///
  /// Both share the single logarithm with the zenith exponent n(p).
  void BuiltinFluxComplete(double p, double &jP, double &n) const {
    const double logP = log(p);
    jP = 1600.*exp(-3.175*log(p + kF1Shift) + 0.279*logP);
    n  = 2.856 - 0.655*logP;
    if (n < 0.1) n = 0.1;
  };

  void BuiltinFluxResidual(double p, double &jP, double &n) const {
    const double logP = log(p);
    jP = 1600.*exp(0.279*logP);
    n  = 2.856 - 0.655*logP;
    if (n < 0.1) n = 0.1;
  };

  void MCJprimeCustomSkyIntegration(double &rate, double &error, int npoints) {
    EMRandom rng = mRandom;
    MCAccumulator acc;
    for (int i = 0; i < npoints; ++i) {
      const double theta = rng.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      const double p     = rng.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      acc.Add(mJ(p, theta)*cos(theta)*sin(theta));
    }
    const double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    acc.Finish(V, npoints, rate, error);
  };

  void MCJprimeCustomCylinderIntegration(double &rate, double &error, int npoints) {
    EMRandom rng = mRandom;
    MCAccumulator acc;
    for (int i = 0; i < npoints; ++i) {
      const double theta = rng.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      const double phi   = rng.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      const double p     = rng.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      const double sinTheta = sin(theta);
      acc.Add(mJ(p, theta)*sinTheta*sinTheta*cos(phi));
    }
    const double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    acc.Finish(V, npoints, rate, error);
  };

  /// Rate per unit area of the half-sphere.
  ///
  /// The local rate at a point of the dome is int dp dOmega J(p,theta)*cos(psi),
  /// with cos(psi) = sin(theta0)sin(theta)cos(phi) + cos(theta0)cos(theta) the
  /// projection on the inward normal and phi the azimuth RELATIVE to that point,
  /// exactly as in Generate(). The rate per unit area is the average of that over
  /// the dome, i.e. over cos(theta0) uniform, so no sin(theta0) Jacobian and no
  /// position-angle factors belong in V.
  void MCJprimeCustomHSphereIntegration(double &rate, double &error, int npoints) {
    EMRandom rng = mRandom;
    MCAccumulator acc;
    for (int i = 0; i < npoints; ++i) {
      const double theta0 = acos(rng.GenerateRandomDouble(mHSphereCosMaxPositionTheta, mHSphereCosMinPositionTheta));
      const double theta  = rng.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      const double phi    = rng.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      const double p      = rng.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      const double sinTheta = sin(theta);
      acc.Add(mJ(p, theta)*(sin(theta0)*sinTheta*cos(phi) + cos(theta0)*cos(theta))*sinTheta);
    }
    const double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    acc.Finish(V, npoints, rate, error);
  };

  /// Rate per unit area of the generation disc, for the target-sphere mode.
  /// The disc is perpendicular to the muon, so the integrand carries only the
  /// solid-angle Jacobian and no projection cosine.
  void MCJprimeCustomTargetIntegration(double &rate, double &error, int npoints) {
    EMRandom rng = mRandom;
    MCAccumulator acc;
    for (int i = 0; i < npoints; ++i) {
      const double theta = rng.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      const double p     = rng.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      acc.Add(mJ(p, theta)*sin(theta));
    }
    const double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    acc.Finish(V, npoints, rate, error);
  };

  /// See MCJprimeCustomTargetIntegration.
  void MCJprimeTargetIntegration(double &rate, double &error, int npoints) {
    EMRandom rng = mRandom;
    MCAccumulator acc;
    double jP, n;
    for (int i = 0; i < npoints; ++i) {
      const double theta = rng.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      const double p     = rng.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      BuiltinFluxComplete(p, jP, n);
      acc.Add(jP*pow(cos(theta), n)*sin(theta));
    }
    const double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    acc.Finish(V, npoints, rate, error);
  };

  void MCJprimeSkyIntegration(double &rate, double &error, int npoints) {
    EMRandom rng = mRandom;
    MCAccumulator acc;
    double jP, n;
    for (int i = 0; i < npoints; ++i) {
      const double theta = rng.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      const double p     = rng.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      BuiltinFluxComplete(p, jP, n);
      acc.Add(jP*pow(cos(theta), n)*cos(theta)*sin(theta));
    }
    const double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    acc.Finish(V, npoints, rate, error);
  };

  void MCJprimeCylinderIntegration(double &rate, double &error, int npoints) {
    EMRandom rng = mRandom;
    MCAccumulator acc;
    double jP, n;
    for (int i = 0; i < npoints; ++i) {
      const double theta = rng.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      const double phi   = rng.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      const double p     = rng.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      BuiltinFluxComplete(p, jP, n);
      const double sinTheta = sin(theta);
      acc.Add(jP*pow(cos(theta), n)*sinTheta*sinTheta*cos(phi));
    }
    const double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    acc.Finish(V, npoints, rate, error);
  };

  /// See MCJprimeCustomHSphereIntegration for the geometry.
  void MCJprimeHSphereIntegration(double &rate, double &error, int npoints) {
    EMRandom rng = mRandom;
    MCAccumulator acc;
    double jP, n;
    for (int i = 0; i < npoints; ++i) {
      const double theta0 = acos(rng.GenerateRandomDouble(mHSphereCosMaxPositionTheta, mHSphereCosMinPositionTheta));
      const double theta  = rng.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
      const double phi    = rng.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      const double p      = rng.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);
      BuiltinFluxComplete(p, jP, n);
      const double sinTheta = sin(theta);
      acc.Add(jP*pow(cos(theta), n)*(sin(theta0)*sinTheta*cos(phi) + cos(theta0)*cos(theta))*sinTheta);
    }
    const double V = (mMaximumMomentum-mMinimumMomentum)*(mMaximumTheta-mMinimumTheta)*(mMaximumPhi-mMinimumPhi);
    acc.Finish(V, npoints, rate, error);
  };

public:
  ///////////////////////////////////////////////////////////////
  /// Generate a cosmic muon from the pre-defined J
  ///////////////////////////////////////////////////////////////
  void Generate() {
    if (mMaxJ[mGenMethod] < 0) ComputeMaximum();

    // Target-sphere generation: only the direction is accept-rejected, the
    // position then follows from it
    if (mGenMethod == TargetSphere) {
      mAccepted = false;
      double jP, n;
      while (!mAccepted) {
        mRandAccRej         = mRandom.GenerateRandomDouble();
        mGenerationTheta    = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationMomentum = GenerateMomentumF1();
        BuiltinFluxResidual(mGenerationMomentum, jP, n);
        mJPrime = jP*pow(cos(mGenerationTheta), n)*sin(mGenerationTheta);
        if (mMaxJ[mGenMethod]*mRandAccRej < mJPrime) mAccepted = true;
      }
      mGenerationTheta = M_PI - mGenerationTheta;
      mGenerationPhi   = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      GeneratePositionTargetSphere();
      GenerateCharge();
      return;
    }

    // Sky or cylinder generation
    if (mGenMethod == Sky || mGenMethod == Cylinder) {
      // For the cylinder a muon that ends up pointing outwards is rejected and the
      // whole event is resampled. That retry used to be a recursive call to
      // Generate(), which drew one extra charge per nesting level and could nest
      // arbitrarily deep; it is a loop now.
      for (;;) {
        // Generation of the momentum and theta angle
        mAccepted = false;
        while (!mAccepted) {
          mRandAccRej  = mRandom.GenerateRandomDouble();
          mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
          mGenerationMomentum = GenerateMomentumF1();
          const double logP = log(mGenerationMomentum);
          mN = 2.856-0.655*logP;
          if (mN < 0.1) mN = 0.1;
          // 1600*pow(p, 0.279) == 1600*exp(0.279*log(p)), reusing the log above
          const double jP = 1600*exp(0.279*logP);

          if (mGenMethod == Sky) {
            mJPrime = jP*pow(cos(mGenerationTheta), mN+1)*sin(mGenerationTheta);
            if (mMaxJ[mGenMethod]*mRandAccRej < mJPrime) mAccepted = true;
          }

          if(mGenMethod == Cylinder)  {
            const double sinTheta = sin(mGenerationTheta);
            mJPrime = jP*pow(cos(mGenerationTheta), mN)*sinTheta*sinTheta;
            if (mMaxJ[mGenMethod]*mRandAccRej < mJPrime) mAccepted = true;
          }
        }
        mGenerationTheta = M_PI - mGenerationTheta;

        // Generation of the position and phi angle
        if (mGenMethod == Sky) {
          GeneratePositionSky();
          mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
          break;
        }

        mAccepted = false;
        GeneratePositionCylinder();
        while (!mAccepted) {
          mRandAccRej  = mRandom.GenerateRandomDouble();
          mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
          if (mRandAccRej < fabs(cos(mGenerationPhi))) mAccepted = true;
        }
        mGenerationPhi = mGenerationPhi + mPhi0;
        if (mGenerationPhi >= 2.*M_PI) mGenerationPhi -= 2.*M_PI;

        // Keep the muon only if it points into the cylinder
        if (!IsOutwardCylinder()) break;
      }
    }

    // Half-sphere generation
    if (mGenMethod == HSphere) {
      mAccepted = false;
      // Generation point on the half-sphere
      mPhi0      = mRandom.GenerateRandomDouble(mHSphereMinPositionPhi, mHSphereMaxPositionPhi);
      while (!mAccepted) {
        mRandAccRej         = mRandom.GenerateRandomDouble();
        mTheta0             = acos(mRandom.GenerateRandomDouble(mHSphereCosMaxPositionTheta, mHSphereCosMinPositionTheta));
        mGenerationTheta    = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationPhi      = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
        mGenerationMomentum = GenerateMomentumF1();
        const double logP   = log(mGenerationMomentum);
        mN                  = 2.856-0.655*logP;
        if (mN < 0.1) mN = 0.1;

        const double sinTheta = sin(mGenerationTheta);
        mJPrime = 1600*exp(0.279*logP)*pow(cos(mGenerationTheta), mN)*(sinTheta*sin(mTheta0)*cos(mGenerationPhi)+cos(mGenerationTheta)*cos(mTheta0))*sinTheta;
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
    GenerateCharge();
  };
  ///////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////
  /// Generate a cosmic muon for the user-defined J
  ///////////////////////////////////////////////////////////////
  void GenerateFromCustomJ() {
    if (!mJ) {
      EMLogger(EMLog::ERROR, "GenerateFromCustomJ() called without a differential flux. Use SetDifferentialFlux() first.", EMLog::EcoMug);
      return;
    }

    if (mMaxCustomJ[mGenMethod] < 0) ComputeMaximumCustomJ();

    // Target-sphere generation, see Generate()
    if (mGenMethod == TargetSphere) {
      mAccepted = false;
      while (!mAccepted) {
        mRandAccRej      = mRandom.GenerateRandomDouble();
        mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        const double proposalDensity = GenerateMomentumCustomJ();
        mJPrime = mJ(mGenerationMomentum, mGenerationTheta)*sin(mGenerationTheta);
        if (mMaxCustomJ[mGenMethod]*mRandAccRej*proposalDensity < mJPrime) mAccepted = true;
      }
      mGenerationTheta = M_PI - mGenerationTheta;
      mGenerationPhi   = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      GeneratePositionTargetSphere();
      GenerateCharge();
      return;
    }

    // Sky or cylinder generation
    if (mGenMethod == Sky || mGenMethod == Cylinder) {
      // As in Generate(), an outward cylinder muon means the whole event is
      // resampled. This used to recurse into Generate(), i.e. into the BUILT-IN
      // flux, which silently mixed built-in muons into a custom-J sample.
      for (;;) {
        // Generation of the momentum and theta angle
        mAccepted = false;
        while (!mAccepted) {
          mRandAccRej  = mRandom.GenerateRandomDouble();
          mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
          const double proposalDensity = GenerateMomentumCustomJ();

          if (mGenMethod == Sky) {
            mJPrime = mJ(mGenerationMomentum, mGenerationTheta)*cos(mGenerationTheta)*sin(mGenerationTheta);
            if (mMaxCustomJ[mGenMethod]*mRandAccRej*proposalDensity < mJPrime) mAccepted = true;
          }

          if(mGenMethod == Cylinder)  {
            // No cos(phi) factor here: phi is not generated yet at this point and
            // the envelope (EMMaximization::CylFunc) does not contain it either.
            // Using the previous event's phi made this loop unable to terminate
            // whenever that cos was negative.
            const double sinTheta = sin(mGenerationTheta);
            mJPrime = mJ(mGenerationMomentum, mGenerationTheta)*sinTheta*sinTheta;
            if (mMaxCustomJ[mGenMethod]*mRandAccRej*proposalDensity < mJPrime) mAccepted = true;
          }
        }
        mGenerationTheta = M_PI - mGenerationTheta;

        // Generation of the position and phi angle
        if (mGenMethod == Sky) {
          GeneratePositionSky();
          mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
          break;
        }

        mAccepted = false;
        GeneratePositionCylinder();
        while (!mAccepted) {
          mRandAccRej  = mRandom.GenerateRandomDouble();
          mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
          if (mRandAccRej < fabs(cos(mGenerationPhi))) mAccepted = true;
        }
        mGenerationPhi = mGenerationPhi + mPhi0;
        if (mGenerationPhi >= 2.*M_PI) mGenerationPhi -= 2.*M_PI;

        // Keep the muon only if it points into the cylinder
        if (!IsOutwardCylinder()) break;
      }
    }

    // Half-sphere generation
    if (mGenMethod == HSphere) {
      mAccepted = false;
      // Generation point on the half-sphere
      mPhi0                 = mRandom.GenerateRandomDouble(mHSphereMinPositionPhi, mHSphereMaxPositionPhi);
      while (!mAccepted) {
        mRandAccRej         = mRandom.GenerateRandomDouble();
        mTheta0             = acos(mRandom.GenerateRandomDouble(mHSphereCosMaxPositionTheta, mHSphereCosMinPositionTheta));
        mGenerationTheta    = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationPhi      = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
        const double proposalDensity = GenerateMomentumCustomJ();

        const double sinTheta = sin(mGenerationTheta);
        mJPrime = mJ(mGenerationMomentum, mGenerationTheta)*(sin(mTheta0)*sinTheta*cos(mGenerationPhi) + cos(mTheta0)*cos(mGenerationTheta))*sinTheta;
        if (mMaxCustomJ[mGenMethod]*mRandAccRej*proposalDensity < mJPrime) mAccepted = true;
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
    GenerateCharge();
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
  mIndex(-1), mSigInstance(signal), mBckInstances{backgrounds}, mWeights(backgrounds.size()+1, 1.), 
  mPID(backgrounds.size()+1, 0) {
};

/// Set the seed for the internal PRNGs (if 0 a random seed is used).
/// Every generator instance gets its own independent stream, derived from seed.
void SetSeed(std::uint64_t seed) {
  if (seed == 0) return;
  std::uint64_t state = seed;
  mRandom.SetSeed(EMRandom::SplitMix64(state));
  mSigInstance.SetSeed(EMRandom::SplitMix64(state));
  for (auto& bck : mBckInstances) bck.SetSeed(EMRandom::SplitMix64(state));
};

/// Set the weights for all EcoMug background instance. The number of elements
/// must be equal to the background instances defined
void SetBckWeights(const std::vector<double>& weights) {
  if (mBckInstances.size() != weights.size()) {
    EMLogger(EMLog::ERROR, "Expected " + std::to_string(mBckInstances.size()) + " weights, but " + std::to_string(weights.size()) + " were provided. Setting them to 1.", EMLog::EMMultiGen);
    std::fill(mWeights.begin(), mWeights.end(), 1.);
  } else {
    for (std::size_t i = 0; i < weights.size(); ++i) mWeights[i+1] = weights[i];
  }
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
  return Selected().GetGenerationPosition();
};

/// Get the generation momentum
double GetGenerationMomentum() const {
  return Selected().GetGenerationMomentum();
};

/// Get the generation momentum
void GetGenerationMomentum(std::array<double, 3>& momentum) const {
  Selected().GetGenerationMomentum(momentum);
};

/// Get the generation theta
double GetGenerationTheta() const {
  return Selected().GetGenerationTheta();
};

/// Get the generation phi
double GetGenerationPhi() const {
  return Selected().GetGenerationPhi();
};

/// Get PID
int GetPID() const {
  // muon case
  if (mPID[Index()] == 0) {
    if (Selected().GetCharge() < 0) return 13;
    else return -13;
  }
  return mPID[Index()];
};

void Generate() {
  mIndex = SelectInstance();
  EcoMug& source = (mIndex == 0) ? mSigInstance : mBckInstances[mIndex-1];
  // Drive each source with the flux it was configured for. This used to call
  // Generate() unconditionally, so a component carrying a user-supplied flux
  // silently produced built-in muons instead.
  if (source.UsesCustomFlux()) source.GenerateFromCustomJ();
  else source.Generate();
};

private:
/// Index of the instance that generated the last muon: 0 is the signal,
/// i > 0 is background i-1. Before the first Generate() call it is the signal.
int Index() const {
  return (mIndex <= 0) ? 0 : mIndex;
};

/// The instance that generated the last muon
const EcoMug& Selected() const {
  if (Index() == 0) return mSigInstance;
  return mBckInstances[Index()-1];
};

/// Draw the instance generating the next muon, with probability proportional
/// to its weight
int SelectInstance() {
  double total = 0.;
  for (std::size_t i = 0; i < mWeights.size(); ++i) total += mWeights[i];
  if (total <= 0.) return 0;
  double x = mRandom.GenerateRandomDouble(0., total);
  double sum = 0.;
  for (std::size_t i = 0; i < mWeights.size(); ++i) {
    sum += mWeights[i];
    if (x < sum) return static_cast<int>(i);
  }
  return static_cast<int>(mWeights.size()) - 1;
};

int mIndex;
EcoMug mSigInstance;
std::vector<EcoMug> mBckInstances;
std::vector<double> mWeights;
std::vector<int> mPID;
EMRandom mRandom;
};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/// Library version as a C++ constant. The ECOMUG_VERSION macro is undefined
/// just below so that it does not leak into user code, which left the version
/// unreachable from outside this header.
inline constexpr const char* EcoMugVersion = ECOMUG_VERSION;

#ifdef ECOMUG_VERSION 
#undef ECOMUG_VERSION
#endif

#ifdef M_PI_NOT_DEFINED
#undef M_PI
#undef M_PI_NOT_DEFINED
#endif

#endif
