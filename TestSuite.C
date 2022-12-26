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
using namespace std; 
using namespace EMUnits; 

Bool_t IsDetCrossed(TVector3 P1, TVector3 P2, TVector3 P3, TVector3 Ro, TVector3 Po, TVector3 &intercept) {
  TMatrixD* num = new TMatrixD(4,4);
  (*num)(0,0) = 1.;       (*num)(0,1) = 1.;       (*num)(0,2) = 1.;       (*num)(0,3) = 1.;
  (*num)(1,0) = P1.X();   (*num)(1,1) = P2.X();   (*num)(1,2) = P3.X();   (*num)(1,3) = Ro.X();
  (*num)(2,0) = P1.Y();   (*num)(2,1) = P2.Y();   (*num)(2,2) = P3.Y();   (*num)(2,3) = Ro.Y();
  (*num)(3,0) = P1.Z();   (*num)(3,1) = P2.Z();   (*num)(3,2) = P3.Z();   (*num)(3,3) = Ro.Z();

  TMatrixD* den = new TMatrixD(4,4);
  (*den)(0,0) = 1.;       (*den)(0,1) = 1.;       (*den)(0,2) = 1.;       (*den)(0,3) = 0.;
  (*den)(1,0) = P1.X();   (*den)(1,1) = P2.X();   (*den)(1,2) = P3.X();   (*den)(1,3) = Po.X();
  (*den)(2,0) = P1.Y();   (*den)(2,1) = P2.Y();   (*den)(2,2) = P3.Y();   (*den)(2,3) = Po.Y();
  (*den)(3,0) = P1.Z();   (*den)(3,1) = P2.Z();   (*den)(3,2) = P3.Z();   (*den)(3,3) = Po.Z();

  if(std::fabs(den->Determinant())<1.e-9) return false;

  Double_t t = -num->Determinant()/den->Determinant();
  intercept = Ro + Po*t;
  if (std::fabs(intercept.Y()) <= P3.Y() && intercept.Z() <= P3.Z() && intercept.Z() >= 0) return true;
  return false;
};

void SuiteNo1(int number_of_events) {
  // ----------------------------------------
  //               Suite No. 1
  // ----------------------------------------
  // 1 m2 horizontal plane detector
  //
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
  //
  // Example: root -l TestSuite.C\(1\,10000\)
  // ----------------------------------------
  
  EcoMug genPlane;
  genPlane.SetUseSky();
  genPlane.SetSkySize({{200.*cm, 200.*cm}});
  genPlane.SetSkyCenterPosition({0., 0., 1.*mm});

  EcoMug genHSphere;
  genHSphere.SetUseHSphere();
  genHSphere.SetHSphereRadius(300*cm);
  genHSphere.SetHSphereCenterPosition({0., 0., 0.});

  TVector3 P1 = {-50.*cm, -50.*cm, 0.};
  TVector3 P2 = { 50.*cm, -50.*cm, 0.};
  TVector3 P3 = { 50.*cm,  50.*cm, 0.};

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

    TVector3 layer_1_v, layer_2_v, layer_3_v;

    if (!IsDetCrossed(P1, P2, P3, muon_origin, muon_p, layer_1_v)) continue;
    n_good_events++;
  }
  cout << "\n--- Generation from horizontal plane ---" << endl;
  cout << "number of generated muons            = " << n_gen_events << endl;
  cout << "number of muons through the detector = " << n_good_events << endl;
  cout << "generation surface are [m2]          = " << genPlane.GetGenSurfaceArea()/m2 << endl;
  cout << "generation rate [Hz/m2]              = " << genPlane.GetAverageGenRate()/hertz*m2 << endl;
  cout << "estimated time of data taking [s]    = " << genPlane.GetEstimatedTime(n_gen_events)/s << endl;

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

    TVector3 layer_1_v, layer_2_v, layer_3_v;

    if (!IsDetCrossed(P1, P2, P3, muon_origin, muon_p, layer_1_v)) continue;
    n_good_events++;
  }
  cout << "\n--- Generation from half-sphere ---" << endl;
  cout << "number of generated muons            = " << n_gen_events << endl;
  cout << "number of muons through the detector = " << n_good_events << endl;
  cout << "generation surface are [m2]          = " << genHSphere.GetGenSurfaceArea()/m2 << endl;
  cout << "(average) generation rate [Hz/m2]    = " << genHSphere.GetAverageGenRate()/hertz*m2 << endl;
  cout << "estimated time of data taking [s]    = " << genHSphere.GetEstimatedTime(n_gen_events)/s << endl;

  cout << endl;
  gApplication->Terminate();
};

void SuiteNo2(int number_of_events) {
  EcoMug genPlane;
  genPlane.SetUseSky();
  genPlane.SetSkySize({{200.*cm, 200.*cm}});
  genPlane.SetSkyCenterPosition({0., 0., 1.*mm});
  genPlane.SetMaximumMomentum(1000*GeV);

  double rate, error;
  genPlane.MCJprimeSkyIntegration(rate, error, number_of_events);
  cout << "rate = " << rate << " +- " << error << endl;

  genPlane.MCJprimeSkyIntegrationStrat(rate, error, number_of_events, 10);
  cout << "rate = " << rate << " +- " << error << endl;

  gApplication->Terminate();
};

void TestSuite(int suite_no, int number_of_events) {

  if (suite_no == 1) {
    return SuiteNo1(number_of_events);
  } else if (suite_no == 2) {
    return SuiteNo2(number_of_events);
  } else {
    cout << "Unknown suite number! Valid values: 1, 2" << endl;
    gApplication->Terminate();
  }
};