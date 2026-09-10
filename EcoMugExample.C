  /////////////////////////////////////////////////////////////////////////////////////
  // Worked examples for the EcoMug cosmic-ray muon generator                        //
  /////////////////////////////////////////////////////////////////////////////////////
  // EcoMug: Efficient COsmic MUon Generator                                         //
  // Copyright (C) 2022 Davide Pagano <davide.pagano@unibs.it>                       //
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
  //                                                                                 //
  // These are worked EXAMPLES, meant to be read and copied from. They print numbers //
  // but they do not check them: the assertion-based test suite is EcoMugTests.C.    //
  //                                                                                 //
  // Run one example with:                                                           //
  //     root -l -b -q 'EcoMugExample.C+(1)'          // default statistics          //
  //     root -l -b -q 'EcoMugExample.C+(2, 20000)'   // example 2, 20000 events     //
  //                                                                                 //
  //   1  a detector in a muon flux: counting rate and live time                     //
  //   2  a two-plane telescope, generated the slow way and the fast way             //
  //      (this one takes about a minute: measuring the slow way IS the point)       //
  //   3  handing muons to Geant4                                                    //
  //   4  supplying your own differential flux                                       //
  //   5  mixing several particle sources with EMMultiGen                            //
  //                                                                                 //
  // The trailing '+' compiles with ACLiC; interpreted it is far slower.             //
  /////////////////////////////////////////////////////////////////////////////////////

#include "EcoMug.h"

#include <TVector3.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <algorithm>

using namespace std;


///////////////////////////////////////////////////////////////////////////////
// Helpers used by the examples
///////////////////////////////////////////////////////////////////////////////

/// A user-supplied differential flux: the standard Gaisser-like parameterisation
/// of the sea-level muon spectrum. Momentum in GeV/c, theta in radians.
double J(double p, double theta) {
  double A = 1400*pow(p, -2.7);
  double B = 1. / (1. + 1.1*p*cos(theta)/115.);
  double C = 0.054 / (1. + 1.1*p*cos(theta)/850.);
  return A*(B+C);
};

class ErrorsUtility {
  public:
  /// Uncertainty on the ratio A/B given the uncertainties on A and B
  static double ErrorRatio(double A, double B, double errA, double errB, double cov = 0) {
    double f = A/B;
    return std::fabs(f)*std::sqrt(std::pow(errA/A, 2) + std::pow(errB/B, 2) - 2*cov/(A*B));
  }
};

/// A flat rectangular detector, defined by three of its corners.
///
///               / y
///              /
///             /         P3
///     -------------------
///    /                 /
///   /                 /--------------- x
///  /                 /
/// /                 /
/// ------------------
/// P1               P2
///
/// P1->P2 and P2->P3 are the two edges, so the rectangle may sit at any
/// orientation in space.
class PlaneDet {
public:
  PlaneDet(const TVector3 &p1, const TVector3 &p2, const TVector3 &p3) :
    mP1(p1), mU(p2 - p1), mV(p3 - p2), mNormal((p2 - p1).Cross(p3 - p2)) {};

  /// Area of the rectangle, correct for any orientation
  double GetArea() const {
    return mNormal.Mag();
  };

  /// True if the muon starting at Ro with momentum Po crosses the rectangle.
  bool IsCrossed(const TVector3 &Ro, const TVector3 &Po) const {
    const double denominator = Po.Dot(mNormal);
    if (std::fabs(denominator) < 1.e-12) return false;   // travelling parallel to the plane
    const double t = (mP1 - Ro).Dot(mNormal)/denominator;
    if (t <= 0.) return false;                           // the plane is behind the muon
    const TVector3 hit = Ro + Po*t - mP1;
    const double alpha = hit.Dot(mU)/mU.Mag2();
    const double beta  = hit.Dot(mV)/mV.Mag2();
    return alpha >= 0. && alpha <= 1. && beta >= 0. && beta <= 1.;
  };

private:
  TVector3 mP1;
  TVector3 mU;       ///< first edge, P1 -> P2
  TVector3 mV;       ///< second edge, P2 -> P3
  TVector3 mNormal;  ///< U x V, whose magnitude is the area
};

namespace {

TVector3 MuonMomentum(const EcoMug& gen) {
  const double p = gen.GetGenerationMomentum();
  const double t = gen.GetGenerationTheta();
  const double f = gen.GetGenerationPhi();
  return TVector3(p*sin(t)*cos(f), p*sin(t)*sin(f), p*cos(t));
};

TVector3 MuonOrigin(const EcoMug& gen) {
  const std::array<double, 3>& x = gen.GetGenerationPosition();
  return TVector3(x[0], x[1], x[2]);
};

double Seconds(std::chrono::steady_clock::time_point a,
               std::chrono::steady_clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
};

}


///////////////////////////////////////////////////////////////////////////////
// 1. A detector in a muon flux: counting rate and live time
///////////////////////////////////////////////////////////////////////////////
void ExampleDetectorRate(int number_of_events) {

  // A 1 m x 1 m horizontal detector lying in the z = 0 plane.
  PlaneDet detector(TVector3(-50.*EMUnits::cm, -50.*EMUnits::cm, 0.),
                    TVector3( 50.*EMUnits::cm, -50.*EMUnits::cm, 0.),
                    TVector3( 50.*EMUnits::cm,  50.*EMUnits::cm, 0.));

  const double sky_side = 2.*EMUnits::m;
  EcoMug gen;
  gen.SetUseSky();
  gen.SetSkySize({{sky_side, sky_side}});
  gen.SetSkyCenterPosition({0., 0., 1.*EMUnits::mm});
  gen.SetMinimumMomentum(500.*EMUnits::MeV);
  gen.SetMaximumMomentum(500.*EMUnits::GeV);
  gen.SetSeed(20250909);

  cout << "\n--- EcoMug v" << EcoMugVersion << ": counting rate of a 1 m x 1 m detector ---\n" << endl;

  long n_generated = 0, n_crossing = 0;
  while (n_crossing < number_of_events) {
    gen.Generate();
    ++n_generated;
    if (detector.IsCrossed(MuonOrigin(gen), MuonMomentum(gen))) ++n_crossing;
  }

  // GetEstimatedTime turns a number of GENERATED muons into the live time they
  // correspond to. It depends on the generator settings and the generation area.
  const double live_time = gen.GetEstimatedTime(n_generated);
  const double rate      = n_crossing/live_time;
  const double det_area  = detector.GetArea()/EMUnits::m2;
  const double sky_area  = gen.GetGenSurfaceArea()/EMUnits::m2;

  double flux = 0., flux_error = 0.;
  gen.GetAverageGenRateAndError(flux, flux_error, 1e7);
  const double rate_error = rate*sqrt(1./n_crossing + pow(flux_error/flux, 2));

  cout << "generation surface [m2]        = " << sky_area << endl;
  cout << "detector area [m2]             = " << det_area << endl;
  cout << "muons generated                = " << n_generated << endl;
  cout << "muons through the detector     = " << n_crossing << endl;
  cout << "live time of the sample [s]    = " << fixed << setprecision(2) << live_time << endl;
  cout << "\nCOUNTING RATE [Hz]             = " << setprecision(2) << rate
       << " +- " << rate_error << endl;
  cout << "              [Hz/m2]          = " << rate/det_area << endl;
  cout << "GetAverageGenRate()  [Hz/m2]   = " << flux << " +- " << flux_error
       << "   <- must agree" << endl;

  cout << "\ngeneration efficiency          = " << setprecision(2)
       << 100.*n_crossing/n_generated << " %" << endl;
  cout << "detector area / sky area       = " << 100.*det_area/sky_area
       << " %   <- must agree" << endl;
};


///////////////////////////////////////////////////////////////////////////////
// 2. A two-plane telescope, generated the slow way and the fast way
///////////////////////////////////////////////////////////////////////////////
void ExampleTargetSphere(int number_of_events) {

  // A muon telescope: two 20 cm x 20 cm planes, 50 cm apart. A muon counts only
  // if it crosses BOTH, which is the coincidence the electronics triggers on.
  const double side = 20.*EMUnits::cm;
  const double gap  = 50.*EMUnits::cm;
  const double h    = side/2.;

  PlaneDet bottom(TVector3(-h, -h, 0.), TVector3(h, -h, 0.), TVector3(h, h, 0.));
  PlaneDet top   (TVector3(-h, -h, gap), TVector3(h, -h, gap), TVector3(h, h, gap));

  // A sphere enclosing the whole telescope.
  const double radius = 0.5*sqrt(2.*side*side + gap*gap);

  cout << "\n--- EcoMug v" << EcoMugVersion << ": two-plane telescope, "
       << number_of_events << " coincidences ---\n" << endl;
  cout << "telescope : two " << side/EMUnits::cm << " x " << side/EMUnits::cm
       << " cm planes, " << gap/EMUnits::cm << " cm apart" << endl;
  cout << "enclosing sphere radius = " << setprecision(4) << radius/EMUnits::cm << " cm\n" << endl;

  cout << setw(16) << "generation" << setw(14) << "generated" << setw(12) << "efficiency"
       << setw(12) << "wall [s]" << setw(22) << "coincidence rate [Hz]" << endl;
  cout << string(76, '-') << endl;

  double rate[2] = {0., 0.}, error[2] = {0., 0.};

  for (int mode = 0; mode < 2; ++mode) {
    EcoMug gen;
    if (mode == 0) {
      // The traditional way: cover the telescope with a plane and throw away
      // every muon that misses. The plane has to be several metres wide even for
      // a 20 cm telescope.
      gen.SetUseSky();
      gen.SetSkySize({{4.*EMUnits::m, 4.*EMUnits::m}});
      gen.SetSkyCenterPosition({0., 0., gap + 1.*EMUnits::cm});
    } else {
      // Aim at the telescope instead. Every muon starts on the enclosing sphere
      // already pointing through it. The cost no longer depends on how large a
      // generation surface would have been needed.
      gen.SetUseTargetSphere();
      gen.SetTargetSphereRadius(radius);
      gen.SetTargetSphereCenterPosition({0., 0., gap/2.});
    }
    gen.SetMinimumMomentum(500.*EMUnits::MeV);
    gen.SetMaximumMomentum(500.*EMUnits::GeV);
    gen.SetSeed(20250909 + mode);

    long n_generated = 0, n_coincidences = 0;
    const auto t0 = chrono::steady_clock::now();
    while (n_coincidences < number_of_events) {
      gen.Generate();
      ++n_generated;
      const TVector3 origin = MuonOrigin(gen);
      const TVector3 p      = MuonMomentum(gen);
      if (top.IsCrossed(origin, p) && bottom.IsCrossed(origin, p)) ++n_coincidences;
    }
    const auto t1 = chrono::steady_clock::now();

    // The rate carries two independent uncertainties: the Poisson one on the
    // coincidence count, and the Monte Carlo one on the flux integral that
    // GetEstimatedTime is built from. Leaving the second out would make the two
    // geometries look inconsistent when they are not.
    double flux = 0., flux_error = 0.;
    gen.GetAverageGenRateAndError(flux, flux_error, 1e7);
    const double live_time = gen.GetEstimatedTime(n_generated);
    rate[mode]  = n_coincidences/live_time;
    error[mode] = rate[mode]*sqrt(1./n_coincidences + pow(flux_error/flux, 2));

    cout << setw(16) << (mode ? "target sphere" : "sky plane")
         << setw(14) << n_generated
         << setw(11) << fixed << setprecision(4) << 100.*n_coincidences/n_generated << "%"
         << setw(12) << setprecision(2) << Seconds(t0, t1)
         << setw(15) << setprecision(4) << rate[mode] << " +- " << error[mode] << endl;
  }
  cout << string(76, '-') << endl;

  const double ratio     = rate[1]/rate[0];
  const double ratio_err = ErrorsUtility::ErrorRatio(rate[1], rate[0], error[1], error[0]);
  cout << "\nrate(target sphere) / rate(sky) = " << setprecision(4) << ratio
       << " +- " << ratio_err << "   ("
       << setprecision(1) << std::fabs(ratio - 1.)/ratio_err << " sigma from 1)" << endl;
};


///////////////////////////////////////////////////////////////////////////////
// 3. Handing muons to Geant4
///////////////////////////////////////////////////////////////////////////////
void ExampleGeant4Handoff(int number_of_events) {

  // The generation surface is a sphere ENCLOSING the apparatus: every muon
  // starts on that sphere already pointing through it, so almost nothing is
  // thrown away, and it starts outside the setup ready to be transported.
  EcoMug gen;
  gen.SetUseTargetSphere();
  gen.SetTargetSphereRadius(80.*EMUnits::cm);          // must contain the whole setup
  gen.SetTargetSphereCenterPosition({0., 0., 50.*EMUnits::cm});
  gen.SetMinimumMomentum(500.*EMUnits::MeV);           // softer muons never reach the tracker
  gen.SetMaximumMomentum(500.*EMUnits::GeV);
  gen.SetSeed(1234);

  cout << "\n--- EcoMug v" << EcoMugVersion
       << ": muons into a G4VUserPrimaryGeneratorAction ---\n" << endl;
  cout << "target sphere R = " << gen.GetTargetSphereRadius()/EMUnits::cm
       << " cm centred at (0, 0, 50) cm" << endl;
  cout << "generation area = " << gen.GetGenSurfaceArea()/EMUnits::m2
       << " m2 (pi R^2, the same for every direction)\n" << endl;

  cout << setw(4) << "#" << setw(6) << "pdg" << setw(30) << "position [cm]"
       << setw(30) << "direction (unit)" << setw(12) << "|p| [GeV/c]"
       << setw(13) << "zenith [deg]" << endl;
  cout << string(95, '-') << endl;

  double max_norm_deviation = 0.;   // largest | |u| - 1 | over the whole sample
  int    n_upward           = 0;    // muons with u_z >= 0, i.e. flying the wrong way

  for (int i = 0; i < number_of_events; ++i) {
    gen.Generate();

    // (1) PDG code. EcoMug returns the electric CHARGE, and the PDG convention
    //     for charged leptons gives the NEGATIVE particle the POSITIVE code:
    //     mu- is 13, mu+ is -13. The code is therefore the opposite sign of the
    //     charge. This is also what EMMultiGen::GetPID does.
    const int pdg = (gen.GetCharge() < 0) ? 13 : -13;

    // (2) Starting point, on the target sphere. EcoMug's default length unit is
    //     the metre while Geant4 works in millimetres, so always attach the unit
    //     (x[0]*CLHEP::m) rather than passing the bare number.
    const std::array<double, 3>& x = gen.GetGenerationPosition();

    // (3) Direction. theta is ALREADY FLIPPED so the muon travels downward, so
    //     (sin t cos f, sin t sin f, cos t) is a unit vector with cos t < 0 and
    //     goes straight into SetParticleMomentumDirection, no sign fixing.
    const double theta = gen.GetGenerationTheta();
    const double phi   = gen.GetGenerationPhi();
    const TVector3 direction(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

    // (4) Momentum magnitude, GeV/c in EcoMug's default units.
    const double ptot = gen.GetGenerationMomentum();

    // The zenith angle from the vertical, the one worth histogramming, is pi - theta.
    const double zenith_deg = (M_PI - theta)/EMUnits::deg;

    max_norm_deviation = std::max(max_norm_deviation, std::fabs(direction.Mag() - 1.));
    if (direction.Z() >= 0.) ++n_upward;

    if (i < std::min(8, number_of_events)) {
      cout << fixed << setprecision(3) << setw(4) << i << setw(6) << pdg
           << setw(10) << x[0]/EMUnits::cm << setw(10) << x[1]/EMUnits::cm
           << setw(10) << x[2]/EMUnits::cm
           << setw(10) << direction.X() << setw(10) << direction.Y()
           << setw(10) << direction.Z()
           << setw(12) << ptot/EMUnits::GeV
           << setw(13) << setprecision(2) << zenith_deg << endl;
    }
  }

  cout << string(95, '-') << endl;
  cout << "muons checked                = " << number_of_events << endl;
  cout << "max | |direction| - 1 |      = " << scientific << setprecision(2) << max_norm_deviation
       << "   (Geant4 wants a unit vector)" << endl;
  cout << "muons with direction z >= 0  = " << n_upward << "   (expected 0: all downward)" << endl;
  cout << fixed << setprecision(3)
       << "live time of this sample [s] = " << gen.GetEstimatedTime(number_of_events)
       << "   (detector counts divided by this = rate in Hz)" << endl;

  // ---------------------------------------------------------------------------
  // The same four quantities inside a Geant4 primary generator action:
  //
  //   void MyPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  //     fEcoMug.Generate();
  //     const int pdg = (fEcoMug.GetCharge() < 0) ? 13 : -13;   // mu- = 13, mu+ = -13
  //     const std::array<double,3> x = fEcoMug.GetGenerationPosition();
  //     const double theta = fEcoMug.GetGenerationTheta();      // already downward
  //     const double phi   = fEcoMug.GetGenerationPhi();
  //     fParticleGun->SetParticleDefinition(
  //         G4ParticleTable::GetParticleTable()->FindParticle(pdg));
  //     fParticleGun->SetParticlePosition(
  //         G4ThreeVector(x[0]*CLHEP::m, x[1]*CLHEP::m, x[2]*CLHEP::m));
  //     fParticleGun->SetParticleMomentumDirection(
  //         G4ThreeVector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)));
  //     fParticleGun->SetParticleMomentum(fEcoMug.GetGenerationMomentum()*CLHEP::GeV);
  //     fParticleGun->GeneratePrimaryVertex(event);
  //   }
  // ---------------------------------------------------------------------------
};


///////////////////////////////////////////////////////////////////////////////
// 4. Supplying your own differential flux
///////////////////////////////////////////////////////////////////////////////
void ExampleCustomFlux(int number_of_events) {

  cout << "\n--- EcoMug v" << EcoMugVersion << ": built-in flux vs a user-supplied one ---\n" << endl;

  const double pmin = 1.*EMUnits::GeV, pmax = 100.*EMUnits::GeV;

  double mean_p[2] = {0., 0.}, mean_p2[2] = {0., 0.};
  double mean_t[2] = {0., 0.}, mean_t2[2] = {0., 0.};
  double seconds[2] = {0., 0.};

  for (int custom = 0; custom < 2; ++custom) {
    EcoMug gen;
    gen.SetUseSky();
    gen.SetSkySize({{2.*EMUnits::m, 2.*EMUnits::m}});
    gen.SetSkyCenterPosition({0., 0., 1.*EMUnits::m});
    gen.SetMinimumMomentum(pmin);
    gen.SetMaximumMomentum(pmax);
    // A user flux is any function of (momentum, theta). It is used by
    // GenerateFromCustomJ(), not by Generate().
    if (custom) gen.SetDifferentialFlux(&J);
    gen.SetSeed(4242);

    const auto t0 = chrono::steady_clock::now();
    for (int i = 0; i < number_of_events; ++i) {
      if (custom) gen.GenerateFromCustomJ(); else gen.Generate();
      const double p = gen.GetGenerationMomentum();
      // theta is returned already flipped, so the zenith angle is pi - theta
      const double t = M_PI - gen.GetGenerationTheta();
      mean_p[custom]  += p;   mean_p2[custom] += p*p;
      mean_t[custom]  += t;   mean_t2[custom] += t*t;
    }
    seconds[custom] = Seconds(t0, chrono::steady_clock::now());

    const double n = number_of_events;
    mean_p[custom] /= n;  mean_p2[custom] = sqrt((mean_p2[custom]/n - mean_p[custom]*mean_p[custom])/n);
    mean_t[custom] /= n;  mean_t2[custom] = sqrt((mean_t2[custom]/n - mean_t[custom]*mean_t[custom])/n);
  }

  cout << "momentum range " << pmin/EMUnits::GeV << " - " << pmax/EMUnits::GeV
       << " GeV/c, " << number_of_events << " muons each\n" << endl;
  cout << setw(22) << "" << setw(24) << "built-in flux" << setw(24) << "user-supplied flux" << endl;
  cout << string(70, '-') << endl;
  cout << setw(22) << "mean momentum [GeV/c]" << fixed << setprecision(4)
       << setw(17) << mean_p[0] << " +-" << setw(7) << mean_p2[0]
       << setw(17) << mean_p[1] << " +-" << setw(7) << mean_p2[1] << endl;
  cout << setw(22) << "mean zenith [rad]"
       << setw(17) << mean_t[0] << " +-" << setw(7) << mean_t2[0]
       << setw(17) << mean_t[1] << " +-" << setw(7) << mean_t2[1] << endl;
  cout << setw(22) << "microseconds per muon"
       << setw(24) << setprecision(3) << 1.e6*seconds[0]/number_of_events
       << setw(24) << 1.e6*seconds[1]/number_of_events << endl;
  cout << string(70, '-') << endl;

  // A user flux carries its own normalisation, which EcoMug cannot know. Two
  // consequences, both deliberate:
  //   - SetHorizontalRate has NO effect on the rate in custom-J mode;
  //   - GetEstimatedTime returns 0 and warns, because a live time cannot be
  //     derived. Normalise the flux yourself if you need absolute rates.
  cout << "\nA user flux carries its own normalisation, so SetHorizontalRate does not" << endl;
  cout << "affect it and GetEstimatedTime cannot work. Asking for it warns:" << endl;
  EcoMug warn;
  warn.SetUseSky();
  warn.SetSkySize({{2.*EMUnits::m, 2.*EMUnits::m}});
  warn.SetMinimumMomentum(pmin);
  warn.SetMaximumMomentum(pmax);
  warn.SetDifferentialFlux(&J);
  const double no_live_time = warn.GetEstimatedTime(10000);   // warns, then returns 0
  cout << "  GetEstimatedTime(10000) = " << no_live_time << " s" << endl;
};


///////////////////////////////////////////////////////////////////////////////
// 5. Mixing several particle sources with EMMultiGen
///////////////////////////////////////////////////////////////////////////////
void ExampleBackgroundMix(int number_of_events) {

  // EMMultiGen draws each event from one of several EcoMug instances.
  //  * SetBckWeights takes one RELATIVE weight per background. The signal is not
  //    in that list because it implicitly carries weight 1, so a component fires
  //    with probability w_i/(1 + sum w_j).
  //  * SetBckPID stamps a PDG code on each background. A code of 0, which the
  //    signal keeps, means "muon, charge chosen by the generator": GetPID then
  //    returns 13 or -13 following the built-in mu+/mu- = 128/100.
  //  * SetSeed seeds the component chooser AND every contained generator, each
  //    with its own independent stream derived from that one number.
  //  * Each source keeps its own configuration, including its own differential
  //    flux: a component set up with SetDifferentialFlux is driven with
  //    GenerateFromCustomJ, the others with Generate.
  const double w_electrons = 0.35, w_protons = 0.05;   // relative to the signal's 1
  const double w_total     = 1. + w_electrons + w_protons;

  // Built in a lambda so the identical mixture can be rebuilt from scratch when
  // reproducibility is checked at the end.
  auto build = [&](unsigned seed) {
    EcoMug muons;                                    // the signal
    muons.SetUseSky();
    muons.SetSkySize({{2.*EMUnits::m, 2.*EMUnits::m}});
    muons.SetSkyCenterPosition({0., 0., 3.*EMUnits::m});
    muons.SetMinimumMomentum(500.*EMUnits::MeV);
    muons.SetMaximumMomentum(500.*EMUnits::GeV);

    EcoMug electrons(muons);                         // a soft electromagnetic component
    electrons.SetMinimumMomentum(50.*EMUnits::MeV);
    electrons.SetMaximumMomentum(2.*EMUnits::GeV);

    EcoMug protons(muons);                           // harder and more vertical
    protons.SetMinimumMomentum(1.*EMUnits::GeV);
    protons.SetMaximumMomentum(20.*EMUnits::GeV);
    protons.SetMaximumTheta(45.*EMUnits::deg);

    EMMultiGen mix(muons, {electrons, protons});
    mix.SetBckWeights({w_electrons, w_protons});
    mix.SetBckPID({11, 2212});
    mix.SetSeed(seed);
    return mix;
  };

  cout << "\n--- EcoMug v" << EcoMugVersion << ": mixing three particle sources ---\n" << endl;
  cout << "relative weights  1 : " << w_electrons << " : " << w_protons
       << "   (signal : electrons : protons)" << endl;

  EMMultiGen mix = build(20250917);

  // The signal's weight is shared between the two muon charges by the 128/100
  // ratio, so all four expected fractions are known in closed form.
  const double mu_plus_share = 128./228.;
  const int    pid[4]      = {13, -13, 11, 2212};
  const char*  name[4]     = {"mu-", "mu+", "e-", "p"};
  const double expected[4] = {(1. - mu_plus_share)/w_total, mu_plus_share/w_total,
                              w_electrons/w_total, w_protons/w_total};

  long counts[4] = {0, 0, 0, 0};
  double sum_p[4] = {0., 0., 0., 0.};
  vector<string> first_events;

  for (int i = 0; i < number_of_events; ++i) {
    mix.Generate();
    const int    code = mix.GetPID();
    const double p    = mix.GetGenerationMomentum();
    for (int k = 0; k < 4; ++k) {
      if (code != pid[k]) continue;
      ++counts[k];
      sum_p[k] += p;
      if (first_events.size() < 6) {
        ostringstream os;
        os << setw(7) << code << setw(6) << name[k] << fixed << setprecision(4)
           << setw(14) << p << setw(12) << M_PI - mix.GetGenerationTheta()
           << setw(12) << mix.GetGenerationPhi();
        first_events.push_back(os.str());
      }
      break;
    }
  }

  cout << "\n" << setw(7) << "pdg" << setw(7) << "name" << setw(10) << "counts"
       << setw(12) << "observed" << setw(12) << "expected" << setw(9) << "pull"
       << setw(14) << "<p> [GeV/c]" << endl;
  cout << string(71, '-') << endl;
  for (int k = 0; k < 4; ++k) {
    const double observed = double(counts[k])/number_of_events;
    const double sigma    = sqrt(expected[k]*(1. - expected[k])/number_of_events);
    cout << setw(7) << pid[k] << setw(7) << name[k] << setw(10) << counts[k]
         << fixed << setprecision(4) << setw(12) << observed << setw(12) << expected[k]
         << setprecision(2) << setw(9) << (observed - expected[k])/sigma
         << setprecision(2) << setw(14) << (counts[k] ? sum_p[k]/counts[k] : 0.) << endl;
  }
  cout << string(71, '-') << endl;
  cout << "The mean momenta differ because each source keeps its own momentum range." << endl;
  cout << "muons of either charge = " << fixed << setprecision(4)
       << double(counts[0] + counts[1])/number_of_events
       << ",  expected " << 1./w_total << endl;
  if (counts[0] > 0) {
    cout << "mu+/mu- ratio          = " << double(counts[1])/counts[0]
         << ",  expected " << 128./100. << endl;
  }

  cout << "\nfirst events (theta is the zenith angle from the vertical):" << endl;
  cout << setw(7) << "pdg" << setw(6) << "name" << setw(14) << "p [GeV/c]"
       << setw(12) << "theta" << setw(12) << "phi" << endl;
  for (const string& line : first_events) cout << line << endl;

  // The whole mixture, selection included, follows the single seed.
  auto sequence = [&](unsigned seed) {
    EMMultiGen m = build(seed);
    vector<double> v;
    for (int i = 0; i < 200; ++i) {
      m.Generate();
      v.push_back(m.GetGenerationMomentum());
      v.push_back(double(m.GetPID()));
    }
    return v;
  };
  cout << "\nsame seed twice        : "
       << (sequence(99) == sequence(99) ? "identical sequence" : "DIFFERENT") << endl;
  cout << "a different seed       : "
       << (sequence(99) != sequence(100) ? "different sequence, as it should be" : "IDENTICAL") << endl;
};


///////////////////////////////////////////////////////////////////////////////
void EcoMugExample(int example_no = 1, int number_of_events = 20000) {
  switch (example_no) {
    case 1: ExampleDetectorRate(number_of_events);   break;
    case 2: ExampleTargetSphere(number_of_events);   break;
    case 3: ExampleGeant4Handoff(number_of_events);  break;
    case 4: ExampleCustomFlux(number_of_events);     break;
    case 5: ExampleBackgroundMix(number_of_events);  break;
    default:
      cout << "Unknown example. Valid values are 1 to 5:\n"
           << "  1  a detector in a muon flux: counting rate and live time\n"
           << "  2  a two-plane telescope, generated the slow way and the fast way\n"
           << "  3  handing muons to Geant4\n"
           << "  4  supplying your own differential flux\n"
           << "  5  mixing several particle sources with EMMultiGen" << endl;
  }
  cout << endl;
};
