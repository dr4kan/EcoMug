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

#include "EcoMug.h"
#include <ctime>
#include <chrono>

using namespace std; 

class ErrorsUtility {
  public:
  static double ErrorRatio(double A, double B, double errA, double errB, double cov = 0) {
    double f = A/B;
    return std::fabs(f)*std::sqrt(std::pow(errA/A, 2) + std::pow(errB/B, 2) - 2*cov/(A*B));
  }
};

double J(double p, double theta) {
  double A = 1400*pow(p, -2.7);
  double B = 1. / (1. + 1.1*p*cos(theta)/115.);
  double C = 0.054 / (1. + 1.1*p*cos(theta)/850.);
  return A*(B+C);
};

class PlaneDet {
  //               / y
  //              /
  //             /         P3
  //     -------------------
  //    /                 /
  //   /                 /--------------- x
  //  /                 /
  // /                 /
  // ------------------
  // P1               P2
public:

  PlaneDet(const TVector3 &p1, const TVector3 &p2, const TVector3 &p3) {
    mP1 = p1; 
    mP2 = p2; 
    mP3 = p3; 
  };

  double GetArea() const {
    double dx = std::max({mP1.X(), mP2.X(), mP3.X()}) - std::min({mP1.X(), mP2.X(), mP3.X()});
    double dy = std::max({mP1.Y(), mP2.Y(), mP3.Y()}) - std::min({mP1.Y(), mP2.Y(), mP3.Y()});
    return dx*dy;
  };

  bool IsCrossed(TVector3 Ro, TVector3 Po) {
    TMatrixD* num = new TMatrixD(4,4);
    (*num)(0,0) = 1.;        (*num)(0,1) = 1.;        (*num)(0,2) = 1.;        (*num)(0,3) = 1.;
    (*num)(1,0) = mP1.X();   (*num)(1,1) = mP2.X();   (*num)(1,2) = mP3.X();   (*num)(1,3) = Ro.X();
    (*num)(2,0) = mP1.Y();   (*num)(2,1) = mP2.Y();   (*num)(2,2) = mP3.Y();   (*num)(2,3) = Ro.Y();
    (*num)(3,0) = mP1.Z();   (*num)(3,1) = mP2.Z();   (*num)(3,2) = mP3.Z();   (*num)(3,3) = Ro.Z();

    TMatrixD* den = new TMatrixD(4,4);
    (*den)(0,0) = 1.;        (*den)(0,1) = 1.;        (*den)(0,2) = 1.;        (*den)(0,3) = 0.;
    (*den)(1,0) = mP1.X();   (*den)(1,1) = mP2.X();   (*den)(1,2) = mP3.X();   (*den)(1,3) = Po.X();
    (*den)(2,0) = mP1.Y();   (*den)(2,1) = mP2.Y();   (*den)(2,2) = mP3.Y();   (*den)(2,3) = Po.Y();
    (*den)(3,0) = mP1.Z();   (*den)(3,1) = mP2.Z();   (*den)(3,2) = mP3.Z();   (*den)(3,3) = Po.Z();

    if(std::fabs(den->Determinant())<1.e-9) return false;

    Double_t t = -num->Determinant()/den->Determinant();
    TVector3 intercept = Ro + Po*t;

    if ((intercept.Y() <= mP3.Y() && intercept.Y() >= mP2.Y()) &&
        (intercept.X() <= mP2.X() && intercept.X() >= mP1.X())) return true;

    return false;
  };

private:
  TVector3 mP1;
  TVector3 mP2;
  TVector3 mP3;
};



void SuiteNo1(int number_of_events) {
  // ----------------------------------------
  //               Suite No. 1
  // ----------------------------------------
  
  EcoMug genPlane;
  genPlane.SetUseSky();
  genPlane.SetSkySize({{200.*EMUnits::cm, 200.*EMUnits::cm}});
  genPlane.SetSkyCenterPosition({0., 0., 1.*EMUnits::mm});

  EcoMug genHSphere;
  genHSphere.SetUseHSphere();
  genHSphere.SetHSphereRadius(200*EMUnits::cm);
  genHSphere.SetHSphereCenterPosition({0., 0., 0.});

  TVector3 P1 = {-50.*EMUnits::cm, -50.*EMUnits::cm, 0.};
  TVector3 P2 = { 50.*EMUnits::cm, -50.*EMUnits::cm, 0.};
  TVector3 P3 = { 50.*EMUnits::cm,  50.*EMUnits::cm, 0.};
  PlaneDet detector(P1, P2, P3);

  auto n_gen_events  = 0;
  auto n_good_events = 0;
  while (n_good_events < number_of_events) {
    genPlane.Generate();
    n_gen_events++;
    std::array<double, 3> muon_position = genPlane.GetGenerationPosition();
    double muon_ptot = genPlane.GetGenerationMomentum();
    double muon_theta = genPlane.GetGenerationTheta();
    double muon_phi = genPlane.GetGenerationPhi();

    TVector3 muon_origin = {muon_position[0], muon_position[1], muon_position[2]};
    TVector3 muon_p  = {
      muon_ptot*sin(muon_theta)*cos(muon_phi),
      muon_ptot*sin(muon_theta)*sin(muon_phi),
      muon_ptot*cos(muon_theta)
    };

    if (!detector.IsCrossed(muon_origin, muon_p)) continue;
    n_good_events++;
  }
  cout << "\n--- Generation from horizontal plane ---" << endl;
  cout << "number of generated muons                   = " << n_gen_events << endl;
  cout << "number of muons through the detector        = " << n_good_events << endl;
  cout << "number of gen muons/generation surface [m2] = " << n_gen_events/(genPlane.GetGenSurfaceArea()/EMUnits::m2) << endl;
  cout << "Estimated time [s]                          = " << genPlane.GetEstimatedTime(n_gen_events) << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  n_gen_events  = 0;
  n_good_events = 0;
  while (n_good_events < number_of_events) {
    genHSphere.Generate();
    n_gen_events++;
    std::array<double, 3> muon_position = genHSphere.GetGenerationPosition();
    double muon_ptot = genHSphere.GetGenerationMomentum();
    double muon_theta = genHSphere.GetGenerationTheta();
    double muon_phi = genHSphere.GetGenerationPhi();

    TVector3 muon_origin = {muon_position[0], muon_position[1], muon_position[2]};
    TVector3 muon_p  = {
      muon_ptot*sin(muon_theta)*cos(muon_phi),
      muon_ptot*sin(muon_theta)*sin(muon_phi),
      muon_ptot*cos(muon_theta)
    };

    if (!detector.IsCrossed(muon_origin, muon_p)) continue;
    n_good_events++;
  }

  cout << "\n--- Generation from half-sphere ---" << endl;
  cout << "number of generated muons                   = " << n_gen_events << endl;
  cout << "number of muons through the detector        = " << n_good_events << endl;
  cout << "number of gen muons/generation surface [m2] = " << n_gen_events/(genHSphere.GetGenSurfaceArea()/EMUnits::m2) << endl;
  cout << "Estimated time [s]                          = " << genHSphere.GetEstimatedTime(n_gen_events) << endl;

  cout << "\nhorizontal to half-spherical rate           = " << (n_good_events/detector.GetArea())/(n_gen_events/genHSphere.GetGenSurfaceArea()) << endl;

  cout << endl;
  gApplication->Terminate();
};



void SuiteNo2(int number_of_events) {
  // ----------------------------------------
  //               Suite No. 2
  // ----------------------------------------

  EcoMug genPlane;
  genPlane.SetUseSky();
  genPlane.SetSkySize({{200.*EMUnits::cm, 200.*EMUnits::cm}});
  genPlane.SetMinimumMomentum(100.*EMUnits::GeV);
  genPlane.SetMaximumMomentum(1000.*EMUnits::GeV);

  EcoMug genCylinder(genPlane);
  genCylinder.SetUseCylinder();
  genCylinder.SetCylinderRadius(100.*EMUnits::cm);
  genCylinder.SetCylinderHeight(10.*EMUnits::m);

  EcoMug genHSphere(genPlane);
  genHSphere.SetUseHSphere();
  genHSphere.SetHSphereRadius(300*EMUnits::cm);

  EcoMug genCustomSky(genPlane);
  genCustomSky.SetDifferentialFlux(&J);

  EcoMug genCustomCylinder(genCylinder);
  genCustomCylinder.SetDifferentialFlux(&J);

  EcoMug genCustomHSphere(genHSphere);
  genCustomHSphere.SetDifferentialFlux(&J);


  double rateSky, rateCyl, rateHS, rateCustomSky, rateCustomCylinder, rateCustomHSphere;
  double errorSky, errorCyl, errorHS, errorCustomSky, errorCustomCylinder, errorCustomHSphere;
  genPlane.GetAverageGenRateAndError(rateSky, errorSky, 1e7);
  genCylinder.GetAverageGenRateAndError(rateCyl, errorCyl, 1e7);
  genHSphere.GetAverageGenRateAndError(rateHS, errorHS, 1e7);
  genCustomSky.GetAverageGenRateAndError(rateCustomSky, errorCustomSky, 1e7);
  genCustomCylinder.GetAverageGenRateAndError(rateCustomCylinder, errorCustomCylinder, 1e7);
  genCustomHSphere.GetAverageGenRateAndError(rateCustomHSphere, errorCustomHSphere, 1e7);

  cout << "rate sky                  = " << rateSky << " +- " << errorSky << endl;
  cout << "rate cylinder             = " << rateCyl << " +- " << errorCyl << endl;
  cout << "rate half-sphere          = " << rateHS << " +- " << errorHS << endl;
  cout << "rate custom J sky         = " << rateCustomSky << " +- " << errorCustomSky << endl;
  cout << "rate custom J cylinder    = " << rateCustomCylinder << " +- " << errorCustomCylinder << endl;
  cout << "rate custom J half-sphere = " << rateCustomHSphere << " +- " << errorCustomHSphere << endl;

  gApplication->Terminate();
};



void SuiteNo3(int number_of_events) {
  // ----------------------------------------
  //               Suite No. 3
  // ----------------------------------------

  EcoMug genSky;
  genSky.SetUseSky();
  genSky.SetSkySize({{200.*EMUnits::cm, 200.*EMUnits::cm}});
  genSky.SetMinimumMomentum(100.*EMUnits::GeV);

  EcoMug genSky200(genSky);
  genSky200.SetHorizontalRate(150.*EMUnits::hertz/EMUnits::m2);

  EcoMug genCustomSky(genSky);
  genCustomSky.SetDifferentialFlux(&J);

  EcoMug genCustomSky200(genSky);
  genCustomSky200.SetDifferentialFlux(&J);
  genCustomSky200.SetHorizontalRate(200.*EMUnits::hertz/EMUnits::m2);

  double rateSky, rateSky200, rateCustomSky, rateCustomSky200;
  double errorSky, errorSky200, errorCustomSky, errorCustomSky200;
  genSky.GetAverageGenRateAndError(rateSky, errorSky, 1e7);
  genSky200.GetAverageGenRateAndError(rateSky200, errorSky200, 1e7);
  genCustomSky.GetAverageGenRateAndError(rateCustomSky, errorCustomSky, 1e7);
  genCustomSky200.GetAverageGenRateAndError(rateCustomSky200, errorCustomSky200, 1e7);

  cout << "rate sky (129 Hz/m2)        = " << rateSky << " +- " << errorSky << endl;
  cout << "rate sky (150 Hz/m2)        = " << rateSky200 << " +- " << errorSky200 << endl;
  cout << "rate custom sky (129 Hz/m2) = " << rateCustomSky << " +- " << errorCustomSky << endl;
  cout << "rate custom sky (150 Hz/m2) = " << rateCustomSky200 << " +- " << errorCustomSky200 << endl;

  gApplication->Terminate();
};



void SuiteNo4(int number_of_events) {
  // ----------------------------------------
  //               Suite No. 4
  // ----------------------------------------

  EcoMug muonGen;
  muonGen.SetUseSky();
  muonGen.SetSkySize({{200.*EMUnits::cm, 200.*EMUnits::cm}});
  muonGen.SetSkyCenterPosition({0., 0., 1.*EMUnits::mm});

  EcoMug electronGen(muonGen);
  electronGen.SetDifferentialFlux(&J);

  EcoMug positronsGen(muonGen);
  electronGen.SetDifferentialFlux(&J);

  EMMultiGen genSuite(muonGen, {electronGen, positronsGen});
  genSuite.SetBckWeights({0.2, 0.1});
  genSuite.SetBckPID({11, -11});

  map<int, int> counts;
  for (auto i = 0; i < number_of_events; ++i) {
    genSuite.Generate();
    counts[genSuite.GetPID()]++;
  }
  cout << std::right << std::setw(5) << "PID" << std::right << std::setw(8) << " counts" 
       << "     ratio" << endl; 
  for (auto& [key, value]: counts) {
    cout << std::right << std::setw(5) << key << std::right << std::setw(8) << value 
         << "   (" << std::setprecision(3) << (double) value/(counts[13]+counts[-13]) << ")" << endl;
  }

  gApplication->Terminate();
};



void TestSuite(int suite_no, int number_of_events) {
  if (suite_no == 1) {
    return SuiteNo1(number_of_events);
  } else if (suite_no == 2) {
    return SuiteNo2(number_of_events);
  } else if (suite_no == 3) {
    return SuiteNo3(number_of_events);
  } else if (suite_no == 4) {
    return SuiteNo4(number_of_events);
  } else {
    cout << "Unknown suite number! Valid values are numbers from 1 to 4" << endl;
    gApplication->Terminate();
  }
};