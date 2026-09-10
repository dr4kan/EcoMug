  /////////////////////////////////////////////////////////////////////////////////////
  // Test suite for the EcoMug cosmic-ray muon generator                             //
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
  // Run with:                                                                       //
  //     root -l -b -q 'EcoMugTests.C+'              (default statistics)            //
  //     root -l -b -q 'EcoMugTests.C+(500000)'      (more statistics)               //
  //                                                                                 //
  // -b is optional: the macro switches to batch graphics on its own while it writes  //
  // the PDF, so it also runs safely as 'root -l EcoMugTests.C+'.                     //
  //                                                                                 //
  // The trailing '+' compiles the macro with ACLiC, which matters: interpreted it    //
  // is roughly an order of magnitude slower.                                        //
  //                                                                                 //
  // Every check has an explicit pass/fail criterion. The macro prints a summary and  //
  // writes EcoMugTests.pdf with the momentum, angular, position and charge           //
  // distributions overlaid on their analytic expectations.                          //
  /////////////////////////////////////////////////////////////////////////////////////

#include "EcoMug.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <TROOT.h>

#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <chrono>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Minimal assertion framework
///////////////////////////////////////////////////////////////////////////////
namespace {

int gPassed = 0, gFailed = 0;
string gSection;
vector<string> gFailures;

void Section(const string& name) {
  gSection = name;
  printf("\n\033[1m== %s ==\033[0m\n", name.c_str());
}

void Report(bool ok, const string& what, const string& detail) {
  if (ok) {
    ++gPassed;
    printf("  \033[32mPASS\033[0m  %-58s %s\n", what.c_str(), detail.c_str());
  } else {
    ++gFailed;
    gFailures.push_back(gSection + " / " + what + "  [" + detail + "]");
    printf("  \033[31mFAIL\033[0m  %-58s %s\n", what.c_str(), detail.c_str());
  }
}

void Check(const string& what, bool ok, const string& detail = "") {
  Report(ok, what, detail);
}

/// |a - b| <= tol
void CheckClose(const string& what, double a, double b, double tol) {
  char buf[160];
  snprintf(buf, sizeof(buf), "got %.10g, expected %.10g (tol %.3g)", a, b, tol);
  Report(fabs(a - b) <= tol, what, buf);
}

/// relative agreement
void CheckRel(const string& what, double a, double b, double relTol) {
  char buf[160];
  snprintf(buf, sizeof(buf), "got %.6g, expected %.6g (%.2f%%, tol %.2f%%)",
           a, b, b != 0. ? 100.*fabs(a-b)/fabs(b) : 0., 100.*relTol);
  Report(fabs(a - b) <= relTol*fabs(b), what, buf);
}

void CheckExact(const string& what, double a, double b) {
  char buf[160];
  snprintf(buf, sizeof(buf), "%.17g vs %.17g", a, b);
  Report(a == b, what, buf);
}

///////////////////////////////////////////////////////////////////////////////
// Physics reference: the built-in EcoMug differential flux
//
//   J(p,theta) = 1600 * (p + 2.68)^-3.175 * p^0.279 * cos^n(theta)
//   n(p)       = max(0.1, 2.856 - 0.655 ln p)
//
// The (p+2.68)^-3.175 factor is the density of the F1 distribution that
// Generate() samples the momentum from; it is part of the flux, not merely a
// sampling device. The same product appears explicitly in the rate integrands.
///////////////////////////////////////////////////////////////////////////////
double JBuiltin(double p, double theta) {
  const double n = std::max(0.1, 2.856 - 0.655*log(p));
  return 1600.*pow(p + 2.68, -3.175)*pow(p, 0.279)*pow(cos(theta), n);
}

/// A realistic user-supplied flux, used to exercise the custom-J code paths.
double JCustom(double p, double theta) {
  const double A = 1400.*pow(p, -2.7);
  const double B = 1./(1. + 1.1*p*cos(theta)/115.);
  const double C = 0.054/(1. + 1.1*p*cos(theta)/850.);
  return A*(B + C);
}

/// Average over the dome of max(0, cos(psi)), the projection of the muon
/// direction on the inward normal. Depends on theta only, so it is tabulated.
struct HSphereFactor {
  vector<double> tab;
  double tmin, tmax;
  HSphereFactor(double a, double b, int n = 2000) : tab(n), tmin(a), tmax(b) {
    const int N0 = 160, NPHI = 320;
    for (int k = 0; k < n; ++k) {
      const double t = tmin + (tmax-tmin)*(k+0.5)/n;
      double acc = 0.;
      for (int i = 0; i < N0; ++i) {
        const double t0 = acos((i+0.5)/N0);          // cos(theta0) uniform => area average
        for (int j = 0; j < NPHI; ++j) {
          const double ph = 2.*M_PI*(j+0.5)/NPHI;
          const double c = sin(t0)*sin(t)*cos(ph) + cos(t0)*cos(t);
          if (c > 0.) acc += c;
        }
      }
      tab[k] = acc/(N0*NPHI);
    }
  }
  double operator()(double t) const {
    int k = int((t - tmin)/(tmax - tmin)*tab.size());
    return tab[std::min(std::max(k, 0), int(tab.size()) - 1)];
  }
};

/// Target density in (p, theta), up to normalisation, for each geometry.
double TargetDensity(int geom, const function<double(double,double)>& J,
                     double p, double t, const HSphereFactor& hsf) {
  if (geom == EcoMug::Sky)      return J(p,t)*cos(t)*sin(t);
  if (geom == EcoMug::Cylinder) return J(p,t)*sin(t)*sin(t);
  // the target-sphere generation disc is perpendicular to the muon, so there is
  // no projection cosine: only the solid-angle Jacobian
  if (geom == EcoMug::TargetSphere) return J(p,t)*sin(t);
  return J(p,t)*sin(t)*hsf(t);
}

/// Predicted bin contents of a 1-D marginal of the target density.
/// which = 0 -> marginal in momentum, which = 1 -> marginal in theta.
vector<double> PredictMarginal(int which, int geom,
                               const function<double(double,double)>& J,
                               const HSphereFactor& hsf, int nbins,
                               double pmin, double pmax, double tmin, double tmax) {
  vector<double> pred(nbins, 0.);
  const int NI = 260;
  for (int b = 0; b < nbins; ++b) {
    const double lo = (which == 0) ? pmin + (pmax-pmin)*b/nbins : tmin + (tmax-tmin)*b/nbins;
    const double hi = (which == 0) ? pmin + (pmax-pmin)*(b+1)/nbins : tmin + (tmax-tmin)*(b+1)/nbins;
    double s = 0.;
    for (int i = 0; i < NI; ++i) {
      const double u = lo + (hi-lo)*(i+0.5)/NI;
      for (int j = 0; j < NI; ++j) {
        if (which == 0) {
          const double t = tmin + (tmax-tmin)*(j+0.5)/NI;
          s += TargetDensity(geom, J, u, t, hsf);
        } else {
          const double p = pmin + (pmax-pmin)*(j+0.5)/NI;
          s += TargetDensity(geom, J, p, u, hsf);
        }
      }
    }
    pred[b] = s;
  }
  return pred;
}

/// chi2 of a histogram against a predicted shape, normalised to the same total.
/// Returns chi2/ndf and fills the prediction histogram for plotting.
double Chi2Shape(TH1D* h, const vector<double>& pred, TH1D* hpred, int& ndf) {
  double tp = 0.;
  for (double v : pred) tp += v;
  const double th = h->Integral();
  double chi2 = 0.;
  ndf = 0;
  for (int b = 0; b < h->GetNbinsX(); ++b) {
    const double e = pred[b]/tp*th;
    if (hpred) hpred->SetBinContent(b+1, e);
    if (e < 25.) continue;                       // keep the Gaussian approximation valid
    const double d = h->GetBinContent(b+1) - e;
    chi2 += d*d/e;
    ++ndf;
  }
  ndf = std::max(ndf - 1, 1);
  return chi2/ndf;
}

/// Configure a generator for one of the three geometries, on a common scale.
void Configure(EcoMug& g, int geom, double pmin, double pmax) {
  if (geom == EcoMug::Sky) {
    g.SetUseSky();
    g.SetSkySize({{4., 4.}});
    g.SetSkyCenterPosition({0., 0., 3.});
  } else if (geom == EcoMug::Cylinder) {
    g.SetUseCylinder();
    g.SetCylinderRadius(2.);
    g.SetCylinderHeight(4.);
    g.SetCylinderCenterPosition({0., 0., 0.});
  } else if (geom == EcoMug::HSphere) {
    g.SetUseHSphere();
    g.SetHSphereRadius(3.);
    g.SetHSphereCenterPosition({0., 0., 0.});
  } else {
    g.SetUseTargetSphere();
    g.SetTargetSphereRadius(1.);
    g.SetTargetSphereCenterPosition({0., 0., 0.});
  }
  g.SetMinimumMomentum(pmin);
  g.SetMaximumMomentum(pmax);
}

const char* GeomName(int g) {
  return (g == EcoMug::Sky)      ? "Sky"
       : (g == EcoMug::Cylinder) ? "Cylinder"
       : (g == EcoMug::HSphere)  ? "HSphere" : "TargetSphere";
}

/// Horizontal square detector of side L centred on the origin, in the z = 0 plane.
bool CrossesDetector(const array<double,3>& x0, double th, double ph, double L) {
  const double dz = cos(th);
  if (fabs(dz) < 1e-12) return false;
  const double t = -x0[2]/dz;
  if (t <= 0.) return false;
  const double x = x0[0] + t*sin(th)*cos(ph);
  const double y = x0[1] + t*sin(th)*sin(ph);
  return fabs(x) <= L/2. && fabs(y) <= L/2.;
}

} // anonymous namespace


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void EcoMugTests(long nevents = 200000) {

  printf("\n\033[1mEcoMug v%s test suite\033[0m   (%ld events per distribution)\n",
         EcoMugVersion, nevents);

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPaperSize(28., 18.);        // match the canvas aspect, else the PDF
  gStyle->SetPadLeftMargin(0.15);        // pages come out mostly empty
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetTitleSize(0.055, "t");
  gStyle->SetLabelSize(0.045, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gStyle->SetTitleOffset(1.45, "y");

  const double pmin = 1., pmax = 50.;
  HSphereFactor hsf(0., M_PI/2.);

  /////////////////////////////////////////////////////////////////////////////
  Section("EMRandom");
  /////////////////////////////////////////////////////////////////////////////
  {
    EMRandom r;
    r.SetSeed(12345);
    double lo = 2., hi = -1.;
    bool everOne = false, everNeg = false;
    for (long i = 0; i < 2000000; ++i) {
      const double v = r.GenerateRandomDouble();
      lo = std::min(lo, v); hi = std::max(hi, v);
      if (v >= 1.) everOne = true;
      if (v < 0.)  everNeg = true;
    }
    Check("GenerateRandomDouble() stays in [0,1)", !everOne && !everNeg,
          Form("min %.3g max %.17g", lo, hi));

    // never exactly 1.0: EMMaximization indexes with floor(rand*size)
    EMRandom q; q.SetSeed(1);
    double top = 0.;
    for (long i = 0; i < 2000000; ++i) top = std::max(top, q.GenerateRandomDouble());
    Check("never returns exactly 1.0 (index safety)", top < 1.0, Form("max %.17g", top));

    // reproducibility
    EMRandom a, b;
    a.SetSeed(999); b.SetSeed(999);
    bool same = true;
    for (int i = 0; i < 1000; ++i) same = same && (a.GenerateRandomDouble() == b.GenerateRandomDouble());
    Check("same seed gives the same stream", same, "");

    EMRandom c, d;
    c.SetSeed(1); d.SetSeed(2);
    bool diff = false;
    for (int i = 0; i < 100; ++i) diff = diff || (c.GenerateRandomDouble() != d.GenerateRandomDouble());
    Check("different seeds give different streams", diff, "");

    // The seeding must not start from a degenerate state. With s[0] = s[1] = seed
    // the first draw was 2*seed/2^64 -- exactly 0 for any seed below 2048 -- and
    // the second was a plain bit-rotation of the seed. The right check is that the
    // first draw of a run is uniformly distributed ACROSS seeds, which is what a
    // batch job doing one seed per job depends on.
    {
      const int NB = 20, NSEED = 20000;
      vector<double> first(NB, 0.), second(NB, 0.);
      double meanFirst = 0.;
      for (int s = 1; s <= NSEED; ++s) {
        EMRandom e; e.SetSeed(s);
        const double v1 = e.GenerateRandomDouble();
        const double v2 = e.GenerateRandomDouble();
        first[std::min(NB-1, int(v1*NB))] += 1.;
        second[std::min(NB-1, int(v2*NB))] += 1.;
        meanFirst += v1;
      }
      const double e0 = double(NSEED)/NB;
      double c1 = 0., c2 = 0.;
      for (int b = 0; b < NB; ++b) {
        c1 += (first[b]-e0)*(first[b]-e0)/e0;
        c2 += (second[b]-e0)*(second[b]-e0)/e0;
      }
      Check("first draw of a run is uniform across seeds", c1/(NB-1) < 2.,
            Form("chi2/ndf = %.3f over seeds 1..%d", c1/(NB-1), NSEED));
      Check("second draw of a run is uniform across seeds", c2/(NB-1) < 2.,
            Form("chi2/ndf = %.3f", c2/(NB-1)));
      CheckClose("first draw of a run averages 0.5 across seeds",
                 meanFirst/NSEED, 0.5, 0.01);
    }

    // seed 0 must not lock the generator in the all-zero absorbing state
    EMRandom z; z.SetSeed(0);
    const double z1 = z.GenerateRandomDouble(), z2 = z.GenerateRandomDouble();
    Check("SetSeed(0) does not produce a frozen state", z1 != 0. && z2 != 0. && z1 != z2,
          Form("%.6f %.6f", z1, z2));

    // uniformity
    EMRandom u; u.SetSeed(4242);
    const int NB = 100;
    vector<double> bins(NB, 0.);
    const long NU = 2000000;
    for (long i = 0; i < NU; ++i) bins[std::min(NB-1, int(u.GenerateRandomDouble()*NB))] += 1.;
    double chi2 = 0.;
    for (double v : bins) { const double e = double(NU)/NB; chi2 += (v-e)*(v-e)/e; }
    Check("uniform over [0,1)", chi2/(NB-1) < 1.6, Form("chi2/ndf = %.3f", chi2/(NB-1)));

    // ranged variant
    EMRandom w; w.SetSeed(7);
    double wlo = 1e9, whi = -1e9, sum = 0.;
    const long NW = 500000;
    for (long i = 0; i < NW; ++i) {
      const double v = w.GenerateRandomDouble(-3., 5.);
      wlo = std::min(wlo, v); whi = std::max(whi, v); sum += v;
    }
    Check("GenerateRandomDouble(a,b) respects its range", wlo >= -3. && whi < 5.,
          Form("[%.4f, %.4f]", wlo, whi));
    CheckClose("GenerateRandomDouble(a,b) has the right mean", sum/NW, 1., 0.02);

    // SplitMix64 is deterministic and mixes
    uint64_t st1 = 12345, st2 = 12345;
    Check("SplitMix64 is deterministic",
          EMRandom::SplitMix64(st1) == EMRandom::SplitMix64(st2), "");
    uint64_t st3 = 0;
    Check("SplitMix64(0) is non-zero", EMRandom::SplitMix64(st3) != 0, "");
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Parameter setters and getters");
  /////////////////////////////////////////////////////////////////////////////
  {
    EcoMug g;
    g.SetMinimumMomentum(0.3);   CheckExact("SetMinimumMomentum / GetMinimumMomentum", g.GetMinimumMomentum(), 0.3);
    g.SetMaximumMomentum(123.);  CheckExact("SetMaximumMomentum / GetMaximumMomentum", g.GetMaximumMomentum(), 123.);
    g.SetMinimumTheta(0.1);      CheckExact("SetMinimumTheta / GetMinimumTheta", g.GetMinimumTheta(), 0.1);
    g.SetMaximumTheta(1.2);      CheckExact("SetMaximumTheta / GetMaximumTheta", g.GetMaximumTheta(), 1.2);
    g.SetMinimumPhi(0.5);        CheckExact("SetMinimumPhi / GetMinimumPhi", g.GetMinimumPhi(), 0.5);
    g.SetMaximumPhi(4.0);        CheckExact("SetMaximumPhi / GetMaximumPhi", g.GetMaximumPhi(), 4.0);
    g.SetHorizontalRate(150.);   CheckExact("SetHorizontalRate / GetHorizontalRate", g.GetHorizontalRate(), 150.);

    g.SetSkySize({{7., 9.}});
    CheckExact("SetSkySize / GetSkySize(0)", g.GetSkySize(0), 7.);
    CheckExact("SetSkySize / GetSkySize(1)", g.GetSkySize(1), 9.);

    g.SetCylinderRadius(2.5);    CheckExact("SetCylinderRadius / GetCylinderRadius", g.GetCylinderRadius(), 2.5);
    g.SetCylinderHeight(6.5);    CheckExact("SetCylinderHeight / GetCylinderHeight", g.GetCylinderHeight(), 6.5);
    g.SetCylinderCenterPosition({1., 2., 3.});
    Check("SetCylinderCenterPosition / Get...", g.GetCylinderCenterPosition() == array<double,3>{1.,2.,3.}, "");

    g.SetHSphereRadius(4.5);     CheckExact("SetHSphereRadius / GetHSphereRadius", g.GetHSphereRadius(), 4.5);
    g.SetHSphereCenterPosition({-1., 0., 2.});
    Check("SetHSphereCenterPosition / Get...", g.GetHSphereCenterPosition() == array<double,3>{-1.,0.,2.}, "");

    EcoMug m;
    m.SetUseSky();      Check("SetUseSky selects the plane",        m.GetGenSurfaceArea() == 0., "");
    m.SetGenerationMethod(EcoMug::Cylinder);
    m.SetCylinderRadius(1.); m.SetCylinderHeight(2.);
    CheckRel("SetGenerationMethod(Cylinder) takes effect", m.GetGenSurfaceArea(), 2.*M_PI*1.*2., 1e-12);
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Generation surface area (closed form)");
  /////////////////////////////////////////////////////////////////////////////
  {
    EcoMug g;
    g.SetUseSky(); g.SetSkySize({{3., 5.}});
    CheckRel("Sky: A = a*b", g.GetGenSurfaceArea(), 15., 1e-12);

    g.SetUseCylinder(); g.SetCylinderRadius(2.); g.SetCylinderHeight(7.);
    CheckRel("Cylinder: A = 2*pi*r*h", g.GetGenSurfaceArea(), 2.*M_PI*2.*7., 1e-12);
    g.SetCylinderMinPositionPhi(0.); g.SetCylinderMaxPositionPhi(M_PI);
    CheckRel("Cylinder with restricted phi: A = dphi*r*h", g.GetGenSurfaceArea(), M_PI*2.*7., 1e-12);

    EcoMug h;
    h.SetUseHSphere(); h.SetHSphereRadius(3.);
    CheckRel("HSphere: A = 2*pi*r^2", h.GetGenSurfaceArea(), 2.*M_PI*9., 1e-12);
    h.SetHSphereMinPositionTheta(0.); h.SetHSphereMaxPositionTheta(M_PI/3.);
    CheckRel("HSphere cap: A = 2*pi*r^2*(1-cos(theta_max))",
             h.GetGenSurfaceArea(), 2.*M_PI*9.*(1. - cos(M_PI/3.)), 1e-12);
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Seeding and reproducibility of complete events");
  /////////////////////////////////////////////////////////////////////////////
  {
    // This is the property a Geant4 user depends on: two runs with the same seed
    // must agree event by event, INCLUDING the charge. The charge used to come
    // from a separate engine that SetSeed() never touched.
    for (int geom = 0; geom <= 3; ++geom) {
      auto run = [&](uint64_t seed) {
        EcoMug g; Configure(g, geom, pmin, pmax); g.SetSeed(seed);
        vector<double> v;
        for (int i = 0; i < 500; ++i) {
          g.Generate();
          const array<double,3>& x = g.GetGenerationPosition();
          v.insert(v.end(), {g.GetGenerationMomentum(), g.GetGenerationTheta(),
                             g.GetGenerationPhi(), x[0], x[1], x[2],
                             double(g.GetCharge())});
        }
        return v;
      };
      const vector<double> a = run(12345), b = run(12345), c = run(54321);
      Check(Form("%s: same seed reproduces every event exactly", GeomName(geom)), a == b, "");
      Check(Form("%s: a different seed gives a different sample", GeomName(geom)), a != c, "");
    }

    // charge specifically, since that was the reported defect
    auto charges = [&](uint64_t seed) {
      EcoMug g; Configure(g, EcoMug::Sky, pmin, pmax); g.SetSeed(seed);
      vector<int> v;
      for (int i = 0; i < 200; ++i) { g.Generate(); v.push_back(g.GetCharge()); }
      return v;
    };
    Check("charge sequence follows SetSeed", charges(777) == charges(777), "");
    Check("charge sequence differs for a different seed", charges(777) != charges(778), "");
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Generated values respect their limits");
  /////////////////////////////////////////////////////////////////////////////
  {
    for (int custom = 0; custom < 2; ++custom) {
      for (int geom = 0; geom <= 3; ++geom) {
        EcoMug g; Configure(g, geom, 2., 20.);
        g.SetMinimumTheta(0.2); g.SetMaximumTheta(1.0);
        g.SetMinimumPhi(0.5);   g.SetMaximumPhi(3.0);
        if (custom) g.SetDifferentialFlux(&JCustom);
        g.SetSeed(31337);
        bool okP = true, okT = true, finite = true;
        for (int i = 0; i < 20000; ++i) {
          if (custom) g.GenerateFromCustomJ(); else g.Generate();
          const double p = g.GetGenerationMomentum();
          const double t = M_PI - g.GetGenerationTheta();   // undo the downward flip
          okP = okP && (p >= 2. && p <= 20.);
          okT = okT && (t >= 0.2 - 1e-12 && t <= 1.0 + 1e-12);
          const array<double,3>& x = g.GetGenerationPosition();
          finite = finite && std::isfinite(p) && std::isfinite(t) &&
                   std::isfinite(g.GetGenerationPhi()) &&
                   std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
        }
        const char* tag = custom ? "customJ" : "builtin";
        Check(Form("%s %s: momentum within [pmin,pmax]", tag, GeomName(geom)), okP, "");
        Check(Form("%s %s: theta within [theta_min,theta_max]", tag, GeomName(geom)), okT, "");
        Check(Form("%s %s: all quantities finite", tag, GeomName(geom)), finite, "");
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Muon charge ratio");
  /////////////////////////////////////////////////////////////////////////////
  {
    EcoMug g; Configure(g, EcoMug::Sky, pmin, pmax); g.SetSeed(2024);
    const long N = std::max(200000L, nevents);
    long npos = 0;
    for (long i = 0; i < N; ++i) { g.Generate(); if (g.GetCharge() > 0) ++npos; }
    const double frac = double(npos)/N, expected = 128./228.;
    const double sigma = sqrt(expected*(1.-expected)/N);
    Check("mu+ fraction = 128/228 (charge ratio 1.28)", fabs(frac-expected) < 5.*sigma,
          Form("%.5f vs %.5f (%.1f sigma)", frac, expected, fabs(frac-expected)/sigma));
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Rate estimation");
  /////////////////////////////////////////////////////////////////////////////
  {
    // reference values obtained by direct numerical integration of the built-in
    // flux, for p in [0.5, 500] GeV/c and the default 129 Hz/m2 normalisation
    const double refSky = 111.221, refHSphere = 64.992, refCylinder = 27.832;

    struct { int geom; double ref; } cases[] = {
      {EcoMug::Sky, refSky}, {EcoMug::Cylinder, refCylinder}, {EcoMug::HSphere, refHSphere}
    };
    for (auto& c : cases) {
      EcoMug g;
      if (c.geom == EcoMug::Sky)      { g.SetUseSky(); g.SetSkySize({{6.,6.}}); g.SetSkyCenterPosition({0.,0.,0.05}); }
      else if (c.geom == EcoMug::Cylinder) { g.SetUseCylinder(); g.SetCylinderRadius(1.); g.SetCylinderHeight(2.); g.SetCylinderCenterPosition({0.,0.,0.}); }
      else                            { g.SetUseHSphere(); g.SetHSphereRadius(1.); g.SetHSphereCenterPosition({0.,0.,0.}); }
      g.SetMinimumMomentum(0.5); g.SetMaximumMomentum(500.); g.SetSeed(11);
      double rate = 0., err = 0.;
      g.GetAverageGenRateAndError(rate, err, 4000000);
      Check(Form("%s rate matches the reference integral", GeomName(c.geom)),
            fabs(rate - c.ref) < 5.*err,
            Form("%.4f +- %.4f vs %.4f (%.1f sigma)", rate, err, c.ref,
                 err > 0. ? fabs(rate-c.ref)/err : 0.));
    }

    // the MC error must shrink as 1/sqrt(npoints)
    EcoMug g; g.SetUseSky(); g.SetSkySize({{6.,6.}}); g.SetSkyCenterPosition({0.,0.,0.05});
    g.SetMinimumMomentum(0.5); g.SetMaximumMomentum(500.); g.SetSeed(11);
    double r1, e1, r2, e2;
    g.GetAverageGenRateAndError(r1, e1, 250000);
    g.GetAverageGenRateAndError(r2, e2, 4000000);
    CheckRel("MC error scales as 1/sqrt(npoints)", e1/e2, 4., 0.25);

    // a rate query must not disturb a seeded generation sequence, and must not
    // overwrite the last generated muon
    auto sequence = [&](bool query) {
      EcoMug h; Configure(h, EcoMug::Sky, pmin, pmax); h.SetSeed(2024);
      vector<double> v;
      for (int i = 0; i < 20; ++i) {
        h.Generate();
        v.push_back(h.GetGenerationMomentum());
        if (query && i == 5) h.GetAverageGenRate(20000);
      }
      return v;
    };
    Check("GetAverageGenRate does not consume the generator's stream",
          sequence(false) == sequence(true), "");

    EcoMug k; Configure(k, EcoMug::HSphere, pmin, pmax); k.SetSeed(5);
    k.Generate();
    const double kp = k.GetGenerationMomentum(), kt = k.GetGenerationTheta();
    const array<double,3> kx = k.GetGenerationPosition();
    k.GetAverageGenRate(20000);
    Check("GetAverageGenRate does not overwrite the last muon",
          kp == k.GetGenerationMomentum() && kt == k.GetGenerationTheta() &&
          kx == k.GetGenerationPosition(), "");

    // The rate must scale linearly with the sea-level normalisation...
    auto scaled = [&](double hr, bool custom) {
      EcoMug h; h.SetUseSky(); h.SetSkySize({{2., 2.}});
      h.SetMinimumMomentum(100.); h.SetMaximumMomentum(1000.);
      if (custom) h.SetDifferentialFlux(&JCustom);
      h.SetHorizontalRate(hr);
      h.SetSeed(3);
      return h.GetAverageGenRate(2000000);
    };
    CheckRel("rate scales linearly with SetHorizontalRate",
             scaled(150., false)/scaled(129., false), 150./129., 1e-9);
    // ...but not for a user-supplied flux, which carries its own normalisation.
    // GetAverageGenRateAndError returns before applying the factor in that case.
    CheckRel("a user-supplied flux ignores SetHorizontalRate (by design)",
             scaled(200., true)/scaled(129., true), 1., 1e-9);

    // GetEstimatedTime is consistent with GetAverageGenRate and the surface area
    EcoMug t; Configure(t, EcoMug::Sky, 0.5, 500.); t.SetSeed(3);
    const double expected = 100000./(t.GetGenSurfaceArea()*t.GetAverageGenRate(2000000));
    CheckRel("GetEstimatedTime = N / (area * rate)", t.GetEstimatedTime(100000), expected, 0.05);
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Geometry consistency (physics invariant)");
  /////////////////////////////////////////////////////////////////////////////
  {
    // The rate through ONE fixed detector must not depend on the generation
    // surface used. The half-sphere is a dome over the plane, so Sky and HSphere
    // must agree. (The cylinder generates on its lateral surface only, so it
    // legitimately misses muons entering through the top and is excluded here.)
    const double L = 1.0;
    const long NCROSS = std::max(20000L, nevents/5);
    double rate[2] = {0., 0.}, sig[2] = {0., 0.};
    for (int pass = 0; pass < 2; ++pass) {
      EcoMug g;
      if (pass == 0) { g.SetUseSky(); g.SetSkySize({{6.,6.}}); g.SetSkyCenterPosition({0.,0.,0.05}); }
      else           { g.SetUseHSphere(); g.SetHSphereRadius(1.); g.SetHSphereCenterPosition({0.,0.,0.}); }
      g.SetMinimumMomentum(0.5); g.SetMaximumMomentum(500.);
      g.SetSeed(20250909u + pass);
      long ngen = 0, ncross = 0;
      while (ncross < NCROSS) {
        g.Generate(); ++ngen;
        if (CrossesDetector(g.GetGenerationPosition(), g.GetGenerationTheta(),
                            g.GetGenerationPhi(), L)) ++ncross;
      }
      // the rate carries two uncertainties: the binomial one on the crossing
      // count, and the Monte Carlo one on the generation rate itself
      double r = 0., e = 0.;
      g.GetAverageGenRateAndError(r, e, 4000000);
      const double t = g.GetEstimatedTime(ngen);
      rate[pass] = ncross/t;
      const double relCross = 1./sqrt(double(ncross));
      const double relRate  = (r > 0.) ? e/r : 0.;
      sig[pass] = rate[pass]*sqrt(relCross*relCross + relRate*relRate);
    }
    const double diff = fabs(rate[0] - rate[1]);
    const double sigd = sqrt(sig[0]*sig[0] + sig[1]*sig[1]);
    Check("Sky and HSphere give the same detector rate", diff < 4.*sigd,
          Form("%.3f +- %.3f vs %.3f +- %.3f  (%.1f sigma)",
               rate[0], sig[0], rate[1], sig[1], sigd > 0. ? diff/sigd : 0.));
    Check("detector rate is physically sensible (50-200 Hz/m2)",
          rate[0] > 50. && rate[0] < 200., Form("%.2f Hz/m2", rate[0]));
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Target-sphere generation");
  /////////////////////////////////////////////////////////////////////////////
  {
    const double R = 0.75;
    const array<double,3> C = {0.3, -0.2, 1.5};

    EcoMug g;
    g.SetUseTargetSphere();
    g.SetTargetSphereRadius(R);
    g.SetTargetSphereCenterPosition(C);
    CheckExact("SetTargetSphereRadius / Get...", g.GetTargetSphereRadius(), R);
    Check("SetTargetSphereCenterPosition / Get...", g.GetTargetSphereCenterPosition() == C, "");
    CheckRel("generation area = pi*R^2", g.GetGenSurfaceArea(), M_PI*R*R, 1e-12);

    // a non-positive radius must be refused, leaving the previous value intact
    g.SetTargetSphereRadius(-1.);
    CheckExact("a non-positive radius is refused", g.GetTargetSphereRadius(), R);

    // every muon must start exactly on the sphere, and head into it
    g.SetMinimumMomentum(0.5); g.SetMaximumMomentum(200.); g.SetSeed(31337);
    double worstOnSphere = 0.;
    long misses = 0;
    const long N = std::max(50000L, nevents/4);
    for (long i = 0; i < N; ++i) {
      g.Generate();
      const array<double,3>& x = g.GetGenerationPosition();
      const double dx = x[0]-C[0], dy = x[1]-C[1], dz = x[2]-C[2];
      worstOnSphere = std::max(worstOnSphere, fabs(sqrt(dx*dx+dy*dy+dz*dz) - R));
      // distance of the sphere centre from the muon line must be below R
      const double th = g.GetGenerationTheta(), ph = g.GetGenerationPhi();
      const double ux = sin(th)*cos(ph), uy = sin(th)*sin(ph), uz = cos(th);
      const double t = -(dx*ux + dy*uy + dz*uz);          // closest approach
      const double cx = dx + t*ux, cy = dy + t*uy, cz = dz + t*uz;
      if (t <= 0. || sqrt(cx*cx+cy*cy+cz*cz) > R*(1.+1e-9)) ++misses;
    }
    Check("every muon starts on the target sphere", worstOnSphere < 1e-9,
          Form("max |r - R| = %.3g m", worstOnSphere));
    Check("every muon is aimed through the target sphere", misses == 0,
          Form("%ld of %ld miss", misses, N));

    // rate normalisation, against the reference integral of J*sin(theta)
    EcoMug r; r.SetUseTargetSphere(); r.SetTargetSphereRadius(1.);
    r.SetTargetSphereCenterPosition({0.,0.,0.});
    r.SetMinimumMomentum(0.5); r.SetMaximumMomentum(500.); r.SetSeed(11);
    double rate = 0., err = 0.;
    r.GetAverageGenRateAndError(rate, err, 4000000);
    Check("target-sphere rate matches the reference integral",
          fabs(rate - 148.749) < 5.*err,
          Form("%.4f +- %.4f vs %.4f (%.1f sigma)", rate, err, 148.749,
               err > 0. ? fabs(rate-148.749)/err : 0.));

    // The decisive check: aiming at a sphere must not change the physics. The
    // rate through a fixed detector has to come out the same as with the plane,
    // it just costs far fewer generated muons to get there.
    const double L = 1.0;
    const long NCROSS = std::max(20000L, nevents/5);
    double dRate[2] = {0., 0.}, dSig[2] = {0., 0.};
    long dGen[2] = {0, 0};
    for (int pass = 0; pass < 2; ++pass) {
      EcoMug h;
      if (pass == 0) { h.SetUseSky(); h.SetSkySize({{6.,6.}}); h.SetSkyCenterPosition({0.,0.,0.05}); }
      else {
        h.SetUseTargetSphere();
        h.SetTargetSphereRadius(L/sqrt(2.));            // encloses the square
        h.SetTargetSphereCenterPosition({0.,0.,0.});
      }
      h.SetMinimumMomentum(0.5); h.SetMaximumMomentum(500.);
      h.SetSeed(20250909u + pass);
      long ngen = 0, ncross = 0;
      while (ncross < NCROSS) {
        h.Generate(); ++ngen;
        if (CrossesDetector(h.GetGenerationPosition(), h.GetGenerationTheta(),
                            h.GetGenerationPhi(), L)) ++ncross;
      }
      double rr = 0., ee = 0.;
      h.GetAverageGenRateAndError(rr, ee, 4000000);
      dRate[pass] = ncross/h.GetEstimatedTime(ngen);
      const double relCross = 1./sqrt(double(ncross));
      const double relRate  = (rr > 0.) ? ee/rr : 0.;
      dSig[pass] = dRate[pass]*sqrt(relCross*relCross + relRate*relRate);
      dGen[pass] = ngen;
    }
    const double dd = fabs(dRate[0]-dRate[1]);
    const double ds = sqrt(dSig[0]*dSig[0] + dSig[1]*dSig[1]);
    Check("target sphere gives the same detector rate as the plane", dd < 4.*ds,
          Form("%.3f +- %.3f vs %.3f +- %.3f  (%.1f sigma)",
               dRate[0], dSig[0], dRate[1], dSig[1], ds > 0. ? dd/ds : 0.));
    Check("target sphere needs far fewer generated muons", dGen[1] < dGen[0]/3,
          Form("%ld generated vs %ld from the plane (%.1fx fewer)",
               dGen[1], dGen[0], double(dGen[0])/dGen[1]));
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Accept-reject envelope");
  /////////////////////////////////////////////////////////////////////////////
  {
    // The envelope is cached per geometry. If a setter that changes the sampled
    // domain does not invalidate it, a stale (too small) envelope silently clips
    // the generated distribution.
    EcoMug a; Configure(a, EcoMug::Sky, 1., 5.); a.SetSeed(1);
    for (int i = 0; i < 200; ++i) a.Generate();       // force the envelope
    a.SetMaximumMomentum(200.);                        // widen after the fact
    double maxP = 0.;
    for (int i = 0; i < 200000; ++i) { a.Generate(); maxP = std::max(maxP, a.GetGenerationMomentum()); }
    Check("widening the momentum range invalidates the cached envelope", maxP > 100.,
          Form("largest momentum after widening = %.2f GeV/c", maxP));

    // the same generator built directly with the wide range must agree
    EcoMug b; Configure(b, EcoMug::Sky, 1., 200.); b.SetSeed(1);
    TH1D hA("hA", "", 30, 1., 60.), hB("hB", "", 30, 1., 60.);
    EcoMug a2; Configure(a2, EcoMug::Sky, 1., 5.); a2.SetSeed(9);
    for (int i = 0; i < 200; ++i) a2.Generate();
    a2.SetMaximumMomentum(200.);
    for (long i = 0; i < 200000; ++i) { a2.Generate(); hA.Fill(a2.GetGenerationMomentum()); }
    for (long i = 0; i < 200000; ++i) { b.Generate();  hB.Fill(b.GetGenerationMomentum()); }
    Check("post-widening spectrum matches a freshly configured generator",
          hA.Chi2Test(&hB, "UU NORM") > 0.001,
          Form("Chi2Test p-value = %.4f", hA.Chi2Test(&hB, "UU NORM")));
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Custom differential flux");
  /////////////////////////////////////////////////////////////////////////////
  {
    // GenerateFromCustomJ on a cylinder used to spin forever whenever the
    // previous event left cos(phi) < 0, because the accept-reject test used a
    // phi that had not been generated yet.
    for (int geom = 0; geom <= 3; ++geom) {
      EcoMug g; Configure(g, geom, 0.5, 200.);
      g.SetDifferentialFlux(&JCustom);
      g.SetSeed(7);
      const auto t0 = chrono::steady_clock::now();
      bool finished = true;
      for (int i = 0; i < 3000; ++i) {
        g.GenerateFromCustomJ();
        if (chrono::duration<double>(chrono::steady_clock::now()-t0).count() > 60.) {
          finished = false; break;
        }
      }
      Check(Form("customJ %s: 3000 events complete without hanging", GeomName(geom)),
            finished, "");
    }

    // calling it without a flux must be reported, not crash
    EcoMug n; Configure(n, EcoMug::Sky, 1., 10.); n.SetSeed(1);
    n.GenerateFromCustomJ();
    Check("GenerateFromCustomJ without SetDifferentialFlux does not crash", true,
          "an error is logged");

    // a custom flux equal to the built-in one must reproduce the built-in spectrum
    EcoMug u; Configure(u, EcoMug::Sky, pmin, pmax); u.SetSeed(4242);
    EcoMug v; Configure(v, EcoMug::Sky, pmin, pmax);
    v.SetDifferentialFlux(&JBuiltin); v.SetSeed(4242);
    TH1D hu("hu", "", 30, pmin, pmax), hv("hv", "", 30, pmin, pmax);
    TH1D tu("tu", "", 30, 0., M_PI/2.), tv("tv", "", 30, 0., M_PI/2.);
    const long NC = std::max(100000L, nevents);
    for (long i = 0; i < NC; ++i) {
      u.Generate();              hu.Fill(u.GetGenerationMomentum()); tu.Fill(M_PI-u.GetGenerationTheta());
      v.GenerateFromCustomJ();   hv.Fill(v.GetGenerationMomentum()); tv.Fill(M_PI-v.GetGenerationTheta());
    }
    Check("custom J = built-in J reproduces the built-in momentum spectrum",
          hu.Chi2Test(&hv, "UU NORM") > 0.001,
          Form("Chi2Test p-value = %.4f", hu.Chi2Test(&hv, "UU NORM")));
    Check("custom J = built-in J reproduces the built-in zenith distribution",
          tu.Chi2Test(&tv, "UU NORM") > 0.001,
          Form("Chi2Test p-value = %.4f", tu.Chi2Test(&tv, "UU NORM")));
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Copy semantics");
  /////////////////////////////////////////////////////////////////////////////
  {
    EcoMug a; Configure(a, EcoMug::Cylinder, 2., 40.);
    a.SetMinimumTheta(0.1); a.SetMaximumTheta(1.1);
    a.SetHorizontalRate(140.); a.SetSeed(31337);
    EcoMug b(a);
    Check("copy preserves the momentum range",
          b.GetMinimumMomentum() == 2. && b.GetMaximumMomentum() == 40., "");
    Check("copy preserves the angular range",
          b.GetMinimumTheta() == 0.1 && b.GetMaximumTheta() == 1.1, "");
    Check("copy preserves the horizontal rate", b.GetHorizontalRate() == 140., "");
    CheckExact("copy preserves the generation surface", b.GetGenSurfaceArea(), a.GetGenSurfaceArea());
    vector<double> va, vb;
    for (int i = 0; i < 100; ++i) { a.Generate(); va.push_back(a.GetGenerationMomentum()); }
    for (int i = 0; i < 100; ++i) { b.Generate(); vb.push_back(b.GetGenerationMomentum()); }
    Check("copy continues the same random stream (documented behaviour)", va == vb, "");
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("EMMultiGen");
  /////////////////////////////////////////////////////////////////////////////
  {
    // Three sources with disjoint momentum bands, so the reported momentum
    // identifies which instance actually generated the muon. The getters used to
    // read the wrong instance, and ran off the end of the vector for the last
    // background.
    auto make = [&](double lo, double hi) {
      EcoMug g; g.SetUseSky(); g.SetSkySize({{4.,4.}}); g.SetSkyCenterPosition({0.,0.,3.});
      g.SetMinimumMomentum(lo); g.SetMaximumMomentum(hi);
      return g;
    };
    EMMultiGen mg(make(1., 2.), {make(100., 200.), make(500., 600.)});
    mg.SetBckWeights({0.5, 0.25});
    mg.SetBckPID({11, -11});
    mg.SetSeed(4242);

    const long N = std::max(60000L, nevents/4);
    long nsig = 0, nb1 = 0, nb2 = 0, outside = 0;
    long pid13 = 0, pid11 = 0, pidm11 = 0, pidOther = 0;
    for (long i = 0; i < N; ++i) {
      mg.Generate();
      const double p = mg.GetGenerationMomentum();
      if      (p >= 1.   && p <= 2.)   ++nsig;
      else if (p >= 100. && p <= 200.) ++nb1;
      else if (p >= 500. && p <= 600.) ++nb2;
      else ++outside;
      const int pid = mg.GetPID();
      if      (abs(pid) == 13) ++pid13;
      else if (pid == 11)      ++pid11;
      else if (pid == -11)     ++pidm11;
      else                     ++pidOther;
    }
    Check("every muon comes from one of the configured sources", outside == 0,
          Form("%ld outside all bands", outside));

    const double tot = 1. + 0.5 + 0.25;
    CheckRel("signal fraction follows its weight",     double(nsig)/N, 1./tot,   0.03);
    CheckRel("background 1 fraction follows its weight", double(nb1)/N, 0.5/tot, 0.04);
    CheckRel("background 2 fraction follows its weight", double(nb2)/N, 0.25/tot, 0.05);

    Check("PID mapping is consistent with the source", pidOther == 0, "");
    CheckRel("PID 11 count matches background 1", double(pid11), double(nb1), 1e-9);
    CheckRel("PID -11 count matches background 2", double(pidm11), double(nb2), 1e-9);
    CheckRel("muon PID count matches the signal", double(pid13), double(nsig), 1e-9);

    // reproducibility
    auto multiRun = [&](uint64_t seed) {
      EMMultiGen m(make(1., 2.), {make(100., 200.), make(500., 600.)});
      m.SetBckWeights({0.5, 0.25}); m.SetBckPID({11, -11}); m.SetSeed(seed);
      vector<double> v;
      for (int i = 0; i < 300; ++i) {
        m.Generate();
        v.insert(v.end(), {m.GetGenerationMomentum(), m.GetGenerationTheta(),
                           m.GetGenerationPhi(), double(m.GetPID())});
      }
      return v;
    };
    Check("EMMultiGen::SetSeed makes runs reproducible", multiRun(99) == multiRun(99), "");
    Check("EMMultiGen with a different seed differs", multiRun(99) != multiRun(100), "");

    // A component carrying a user-supplied flux must be driven with that flux.
    // EMMultiGen used to call Generate() unconditionally, silently substituting
    // the built-in muon flux for the user's.
    {
      auto spike = [](double p, double) { return pow(p, -8.); };   // extreme on purpose
      auto plain = [&]() {
        EcoMug g; g.SetUseSky(); g.SetSkySize({{2.,2.}});
        g.SetSkyCenterPosition({0.,0.,1.});
        g.SetMinimumMomentum(1.); g.SetMaximumMomentum(50.);
        return g;
      };
      EcoMug sig = plain();
      sig.SetDifferentialFlux(spike);
      EMMultiGen m(sig, {plain()});
      m.SetBckWeights({1e-4});           // essentially always the signal
      m.SetSeed(1);
      double viaMix = 0.;
      const int NM = 20000;
      for (int i = 0; i < NM; ++i) { m.Generate(); viaMix += m.GetGenerationMomentum(); }
      viaMix /= NM;

      EcoMug direct = plain();
      direct.SetDifferentialFlux(spike);
      direct.SetSeed(1);
      double viaDirect = 0.;
      for (int i = 0; i < NM; ++i) { direct.GenerateFromCustomJ(); viaDirect += direct.GetGenerationMomentum(); }
      viaDirect /= NM;

      EcoMug builtin = plain();
      builtin.SetSeed(1);
      double viaBuiltin = 0.;
      for (int i = 0; i < NM; ++i) { builtin.Generate(); viaBuiltin += builtin.GetGenerationMomentum(); }
      viaBuiltin /= NM;

      CheckRel("EMMultiGen honours a component's custom flux", viaMix, viaDirect, 0.05);
      Check("EMMultiGen does not fall back to the built-in flux",
            fabs(viaMix - viaBuiltin) > 0.5*fabs(viaBuiltin - viaDirect),
            Form("mix %.3f, custom %.3f, built-in %.3f GeV/c", viaMix, viaDirect, viaBuiltin));
    }

    // wrong-sized weight/PID vectors must be handled, not crash
    EMMultiGen bad(make(1., 2.), {make(3., 4.)});
    bad.SetBckWeights({1., 2., 3.});
    bad.SetBckPID({1, 2, 3});
    bad.SetSeed(1);
    bad.Generate();
    Check("mismatched weight/PID vectors are rejected without crashing",
          std::isfinite(bad.GetGenerationMomentum()), "");
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Distributions vs analytic expectation");
  /////////////////////////////////////////////////////////////////////////////

  vector<TH1D*> momHists, thetaHists, momPred, thetaPred;
  vector<TH1D*> phiHists;
  vector<string> caseNames;

  for (int custom = 0; custom < 2; ++custom) {
    for (int geom = 0; geom <= 3; ++geom) {
      const string name = string(custom ? "customJ " : "builtin ") + GeomName(geom);
      caseNames.push_back(name);

      EcoMug g; Configure(g, geom, pmin, pmax);
      function<double(double,double)> J = custom ? function<double(double,double)>(&JCustom)
                                                 : function<double(double,double)>(&JBuiltin);
      if (custom) g.SetDifferentialFlux(&JCustom);
      g.SetSeed(4242);

      const int NB = 25;
      TH1D* hp = new TH1D(Form("hp%d%d", custom, geom),
                          Form("%s;p [GeV/c];muons / bin", name.c_str()), NB, pmin, pmax);
      TH1D* ht = new TH1D(Form("ht%d%d", custom, geom),
                          Form("%s;#theta [rad];muons / bin", name.c_str()), NB, 0., M_PI/2.);
      TH1D* hf = new TH1D(Form("hf%d%d", custom, geom),
                          Form("%s;#varphi [rad];muons / bin", name.c_str()), 36, 0., 2.*M_PI);
      hp->Sumw2(); ht->Sumw2(); hf->Sumw2();

      for (long i = 0; i < nevents; ++i) {
        if (custom) g.GenerateFromCustomJ(); else g.Generate();
        hp->Fill(g.GetGenerationMomentum());
        ht->Fill(M_PI - g.GetGenerationTheta());
        hf->Fill(g.GetGenerationPhi());
      }

      TH1D* pp = new TH1D(Form("pp%d%d", custom, geom), "", NB, pmin, pmax);
      TH1D* pt = new TH1D(Form("pt%d%d", custom, geom), "", NB, 0., M_PI/2.);
      int ndfP = 0, ndfT = 0;
      const double c2p = Chi2Shape(hp, PredictMarginal(0, geom, J, hsf, NB, pmin, pmax, 0., M_PI/2.), pp, ndfP);
      const double c2t = Chi2Shape(ht, PredictMarginal(1, geom, J, hsf, NB, pmin, pmax, 0., M_PI/2.), pt, ndfT);

      Check(Form("%-18s momentum spectrum matches the analytic target", name.c_str()),
            c2p < 3., Form("chi2/ndf = %.3f (%d dof)", c2p, ndfP));
      Check(Form("%-18s zenith distribution matches the analytic target", name.c_str()),
            c2t < 3., Form("chi2/ndf = %.3f (%d dof)", c2t, ndfT));

      momHists.push_back(hp);   momPred.push_back(pp);
      thetaHists.push_back(ht); thetaPred.push_back(pt);
      phiHists.push_back(hf);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("High-momentum regime (clamped zenith exponent)");
  /////////////////////////////////////////////////////////////////////////////
  {
    // n(p) = 2.856 - 0.655 ln p goes negative above about 67 GeV/c and is clamped
    // at 0.1. The distribution checks above run at 1-50 GeV/c and never reach it,
    // so this exercises the clamped branch explicitly.
    const double hpmin = 100., hpmax = 1000.;
    Check("the clamp is actually active over this range",
          2.856 - 0.655*log(hpmin) < 0.1,
          Form("n(%.0f) = %+.4f, n(%.0f) = %+.4f",
               hpmin, 2.856-0.655*log(hpmin), hpmax, 2.856-0.655*log(hpmax)));

    for (int geom = 0; geom <= 3; ++geom) {
      EcoMug g; Configure(g, geom, hpmin, hpmax); g.SetSeed(808);
      const int NB = 20;
      TH1D hp("hhp", "", NB, hpmin, hpmax), ht("hht", "", NB, 0., M_PI/2.);
      hp.Sumw2(); ht.Sumw2();
      const long N = std::max(100000L, nevents);
      for (long i = 0; i < N; ++i) {
        g.Generate();
        hp.Fill(g.GetGenerationMomentum());
        ht.Fill(M_PI - g.GetGenerationTheta());
      }
      function<double(double,double)> J = &JBuiltin;
      int ndfP = 0, ndfT = 0;
      const double c2p = Chi2Shape(&hp, PredictMarginal(0, geom, J, hsf, NB, hpmin, hpmax, 0., M_PI/2.), nullptr, ndfP);
      const double c2t = Chi2Shape(&ht, PredictMarginal(1, geom, J, hsf, NB, hpmin, hpmax, 0., M_PI/2.), nullptr, ndfT);
      Check(Form("%-13s momentum spectrum at 100-1000 GeV/c", GeomName(geom)),
            c2p < 3., Form("chi2/ndf = %.3f (%d dof)", c2p, ndfP));
      Check(Form("%-13s zenith distribution at 100-1000 GeV/c", GeomName(geom)),
            c2t < 3., Form("chi2/ndf = %.3f (%d dof)", c2t, ndfT));
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  Section("Plots");
  /////////////////////////////////////////////////////////////////////////////
  {
    // The plots are only ever written to a PDF, so no on-screen canvas is
    // wanted. Forcing batch mode also avoids a ROOT/macOS crash: opening a real
    // canvas window goes through TGCocoa::SetApplicationIcon(), which faults
    // while loading the application icon when the macro is run without -b.
    const Bool_t wasBatch = gROOT->IsBatch();
    gROOT->SetBatch(kTRUE);

    const char* out = "EcoMugTests.pdf";
    TCanvas c("c", "EcoMug tests", 1400, 900);
    c.Print(Form("%s[", out));

    // --- momentum spectra
    c.Clear(); c.Divide(4, 2);
    for (size_t i = 0; i < momHists.size(); ++i) {
      c.cd(i+1); gPad->SetLogy();
      momHists[i]->SetLineColor(kBlack); momHists[i]->SetMarkerStyle(20); momHists[i]->SetMarkerSize(0.6);
      momHists[i]->Draw("E");
      momPred[i]->SetLineColor(kRed+1); momPred[i]->SetLineWidth(2); momPred[i]->Draw("HIST SAME");
      if (i == 0) {
        TLegend* l = new TLegend(0.45, 0.72, 0.88, 0.88);
        l->AddEntry(momHists[i], "generated", "lep");
        l->AddEntry(momPred[i], "analytic target", "l");
        l->SetBorderSize(0); l->Draw();
      }
    }
    c.cd(0);
    TLatex tt; tt.SetNDC(); tt.SetTextSize(0.022);
    tt.DrawLatex(0.01, 0.985, "EcoMug momentum spectra: generated sample vs analytic target");
    c.Print(out);

    // --- zenith angle
    c.Clear(); c.Divide(4, 2);
    for (size_t i = 0; i < thetaHists.size(); ++i) {
      c.cd(i+1);
      thetaHists[i]->SetLineColor(kBlack); thetaHists[i]->SetMarkerStyle(20); thetaHists[i]->SetMarkerSize(0.6);
      thetaHists[i]->Draw("E");
      thetaPred[i]->SetLineColor(kRed+1); thetaPred[i]->SetLineWidth(2); thetaPred[i]->Draw("HIST SAME");
    }
    c.cd(0);
    tt.DrawLatex(0.01, 0.985, "EcoMug zenith-angle distributions: generated sample vs analytic target");
    c.Print(out);

    // --- azimuth
    c.Clear(); c.Divide(4, 2);
    for (size_t i = 0; i < phiHists.size(); ++i) {
      c.cd(i+1);
      phiHists[i]->SetMinimum(0.);
      phiHists[i]->SetLineColor(kBlue+1);
      phiHists[i]->Draw("E");
    }
    c.cd(0);
    tt.DrawLatex(0.01, 0.985, "EcoMug azimuth distributions");
    c.Print(out);

    // --- the cos^n(theta) law, in momentum slices
    c.Clear(); c.Divide(2, 2);
    const double slices[4][2] = {{1., 2.}, {2., 5.}, {5., 15.}, {15., 50.}};
    for (int s = 0; s < 4; ++s) {
      EcoMug g; Configure(g, EcoMug::Sky, slices[s][0], slices[s][1]); g.SetSeed(101 + s);
      TH1D* h = new TH1D(Form("cs%d", s),
                         Form("p in [%.0f, %.0f] GeV/c;cos#theta;muons / bin",
                              slices[s][0], slices[s][1]), 24, 0., 1.);
      h->Sumw2();
      for (long i = 0; i < nevents; ++i) { g.Generate(); h->Fill(cos(M_PI - g.GetGenerationTheta())); }
      // expected shape: cos^(n+1) from the flux times the cos from the horizontal
      // projection, integrated over the momentum slice
      TH1D* pr = new TH1D(Form("csp%d", s), "", 24, 0., 1.);
      for (int b = 1; b <= 24; ++b) {
        const double ct = pr->GetBinCenter(b);
        const double th = acos(ct);
        double acc = 0.;
        for (int i = 0; i < 400; ++i) {
          const double p = slices[s][0] + (slices[s][1]-slices[s][0])*(i+0.5)/400;
          acc += JBuiltin(p, th)*ct;
        }
        pr->SetBinContent(b, acc);
      }
      pr->Scale(h->Integral()/pr->Integral());
      c.cd(s+1);
      h->SetMarkerStyle(20); h->SetMarkerSize(0.6); h->Draw("E");
      pr->SetLineColor(kRed+1); pr->SetLineWidth(2); pr->Draw("HIST SAME");
    }
    c.cd(0);
    tt.DrawLatex(0.01, 0.985, "Zenith dependence in momentum slices: cos^{n(p)}#theta law");
    c.Print(out);

    // --- generation positions
    c.Clear(); c.Divide(3, 4);
    for (int geom = 0; geom <= 3; ++geom) {
      EcoMug g; Configure(g, geom, pmin, pmax); g.SetSeed(55);
      TH1D* hx = new TH1D(Form("px%d", geom), Form("%s;x [m];muons", GeomName(geom)), 40, -4., 4.);
      TH1D* hy = new TH1D(Form("py%d", geom), Form("%s;y [m];muons", GeomName(geom)), 40, -4., 4.);
      TH1D* hz = new TH1D(Form("pz%d", geom), Form("%s;z [m];muons", GeomName(geom)), 40, -4., 4.);
      for (long i = 0; i < nevents; ++i) {
        g.Generate();
        const array<double,3>& x = g.GetGenerationPosition();
        hx->Fill(x[0]); hy->Fill(x[1]); hz->Fill(x[2]);
      }
      c.cd(3*geom+1); hx->SetLineColor(kAzure+2); hx->Draw("HIST");
      c.cd(3*geom+2); hy->SetLineColor(kAzure+2); hy->Draw("HIST");
      c.cd(3*geom+3); hz->SetLineColor(kAzure+2); hz->Draw("HIST");
    }
    c.cd(0);
    tt.DrawLatex(0.01, 0.985, "Generation positions on the four generation surfaces");
    c.Print(out);

    // --- charge
    c.Clear();
    {
      EcoMug g; Configure(g, EcoMug::Sky, pmin, pmax); g.SetSeed(2024);
      TH1D* h = new TH1D("chg", "muon charge;charge;muons", 2, -2., 2.);
      for (long i = 0; i < nevents; ++i) { g.Generate(); h->Fill(g.GetCharge()); }
      h->SetFillColor(kAzure-9); h->Draw("HIST");
      const double ratio = h->GetBinContent(2)/h->GetBinContent(1);
      TLatex l; l.SetNDC(); l.SetTextSize(0.035);
      l.DrawLatex(0.15, 0.85, Form("#mu^{+}/#mu^{-} = %.4f  (expected 1.2800)", ratio));
      c.Print(out);
    }

    c.Print(Form("%s]", out));
    gROOT->SetBatch(wasBatch);
    Check("plots written to EcoMugTests.pdf", true, "6 pages");
  }

  /////////////////////////////////////////////////////////////////////////////
  // Summary
  /////////////////////////////////////////////////////////////////////////////
  printf("\n\033[1m================ SUMMARY ================\033[0m\n");
  printf("  passed : \033[32m%d\033[0m\n", gPassed);
  printf("  failed : %s%d\033[0m\n", gFailed ? "\033[31m" : "\033[32m", gFailed);
  if (gFailed) {
    printf("\n  failing checks:\n");
    for (const string& f : gFailures) printf("    - %s\n", f.c_str());
  }
  printf("\033[1m=========================================\033[0m\n\n");
}
