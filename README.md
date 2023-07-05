# EcoMug: Efficient COsmic MUon Generator

EcoMug is a header-only C++11 library for the generation of cosmic ray (CR) muons, based on a parametrization of experimental data. Unlike other tools, EcoMug gives the possibility of generating from different surfaces (plane, cylinder and half-sphere), while keeping the correct angular and momentum distribution of generated tracks. EcoMug also allows the generation of CR muons according to user-defined parametrizations of their differential flux.

If you use, or want to refer to, EcoMug please cite the following paper:

> Pagano, D., Bonomi, G., Donzella, A., Zenoni, A., Zumerle, G., & Zurlo, N. (2021). EcoMug: An Efficient COsmic MUon Generator for cosmic-ray muon applications. Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, 1014, 165732.

Latest release: [EcoMug v2.1](https://github.com/dr4kan/EcoMug/releases/tag/v2.1)



# Basic Usage

The use of the library requires the initialization of the `EcoMug` class, the choice of the generation method, and the definition of the size and position of the generation surface. Once the setup of the instance of the `EcoMug` class is done, the generation of a cosmic-ray muon can be invoked with the method `Generate()`, which will compute its position, direction, momentum, and charge. All these quantities can be accessed with the methods `GetGenerationPosition()`, `GetGenerationTheta()`, `GetGenerationPhi()`, `GetGenerationMomentum()`, and `GetCharge()`, as shown in the examples below. The charge for generated muons takes into account the excess of positive muons over negative ones, assuming a constant charge ratio (see the above mentioned paper for more details). Angles are in radians, momentum is in GeV/c, whereas the unit of measure of the position is arbitrary and depends on the choice done in the simulation code where EcoMug is used.



### Plane-based generation

```
EcoMug gen; // initialization of the class
gen.SetUseSky(); // plane surface generation
gen.SetSkySize({{10., 10.}}); // x and y size of the plane
// (x,y,z) position of the center of the plane
gen.SetSkyCenterPosition({{0., 0., 20.}});

// The array storing muon generation position
std::array<double, 3> muon_position;

for (auto event = 0; event < number_of_events; ++event) {
  gen.Generate();
  muon_position = gen.GetGenerationPosition();
  double muon_p = gen.GetGenerationMomentum();
  double muon_theta = gen.GetGenerationTheta();
  double muon_phi = gen.GetGenerationPhi();
  double muon_charge = gen.GetCharge();
  ...
}
```



### Cylinder-based generation

```
EcoMug gen; // initialization of the class
gen.SetUseCylinder(); // cylindrical surface generation
gen.SetCylinderRadius(10.); // cylinder radius
gen.SetCylinderHeight(30.); // cylinder height
// (x,y,z) position of the center of the cylinder
gen.SetCylinderCenterPosition({{0., 0., 15.}});

// The array storing muon generation position
std::array<double, 3> muon_position;

for (auto event = 0; event < number_of_events; ++event) {
  gen.Generate();
  muon_position = gen.GetGenerationPosition();
  double muon_p = gen.GetGenerationMomentum();
  double muon_theta = gen.GetGenerationTheta();
  double muon_phi = gen.GetGenerationPhi();
  double muon_charge = gen.GetCharge();
  ...
}
```



### Hsphere-based generation

```
EcoMug gen; // initialization of the class
gen.SetUseHSphere(); // half-spherical surface generation
gen.SetHSphereRadius(30.); // half-sphere radius
// (x,y,z) position of the center of the half-sphere
gen.SetHSphereCenterPosition({{0., 0., 0.}});

// The array storing muon generation position
std::array<double, 3> muon_position;

for (auto event = 0; event < number_of_events; ++event) {
  gen.Generate();
  muon_position = gen.GetGenerationPosition();
  double muon_p = gen.GetGenerationMomentum();
  double muon_theta = gen.GetGenerationTheta();
  double muon_phi = gen.GetGenerationPhi();
  double muon_charge = gen.GetCharge();
  ...
}
```



# More Advanced Usage

It is possible to set the seed in EcoMug, for reproducible generations. This can be done with the method `SetSeed`, as shown in the example below. If the seed is set to 0 (or the method is not invoked at all), a random seed is used.

```
EcoMug gen;
gen.SetUseSky();
gen.SetSkySize({{10., 10.}});
gen.SetSkyCenterPosition({{0., 0., 20.}});

// set the seed (only positive integers are accepted)
gen.SetSeed(1234);
```

In several scenarios, one could be interested in generating tracks from a subset of these parameters, saving space and computation time. EcoMug allows this by exposing to the user the following methods:

- `SetMinimumMomentum` - Set the minimum momentum for generated cosmic-ray muons;
- `SetMaximumMomentum` - Set the maximum momentum for generated cosmic-ray muons;
- `SetMinimumTheta` - Set the minimum zenith angle 𝜃 for generated cosmic-ray muons;
- `SetMaximumTheta` - Set the maximum zenith angle 𝜃 for generated cosmic-ray muons;
- `SetMinimumPhi` - Set the minimum azimuthal angle 𝜙 for generated cosmic-ray muons;
- `SetMaximumPhi` - Set the maximum azimuthal angle 𝜙 for generated cosmic-ray muons.

```
#include <math.h> // necessary for the M_PI constant

EcoMug gen;
gen.SetUseSky();
gen.SetSkySize({{10., 10.}});
gen.SetSkyCenterPosition({{0., 0., 20.}});

gen.SetMinimumMomentum(80.);
gen.SetMaximumMomentum(800.);
gen.SetMinimumTheta(0.);
gen.SetMaximumTheta(M_PI/4);
gen.SetMinimumPhi(0.);
gen.SetMaximumPhi(M_PI);
```

In those cases where the proposed parametrization of the differential flux *J* of CR muons does not fit the user needs, EcoMug gives the possibility to use a custom function for *J*, as shown in the example below.

```
double J(double p, double theta) {
  double A = 0.14*pow(p, -2.7);
  double B = 1. / (1. + 1.1*p*cos(theta)/115.);
  double C = 0.054 / (1. + 1.1*p*cos(theta)/850.);
  return A*(B+C);
}

EcoMug gen;
gen.SetUseSky();
gen.SetSkySize({{x, y}});
gen.SetSkyCenterPosition({0., 0., z});
gen.SetMinimumMomentum(150);
gen.SetDifferentialFlux(&J);

for (auto event = 0; event < nevents; ++event) {
  gen.GenerateFromCustomJ(); // generate from user-defined J
  ... // retrieve and use muon data
  gen.Generate(); // generate from J as in equation 2
  ... // retrieve and use muon data
}
```


# Rate and time estimation

EcoMug allows to estimate the rate and time to collect a given number of muons, also in those cases where the user has constrained the generation (for example by cutting on the momentum or angles). The user can specify the average expected rate (via `SetHorizontalRate`) to take into account, for example, the effect of altitude. Default value is 129 $Hz/m^2$.

While the rate and time estimation also works with custom definitions of the flux, it is up to the user to define a properly normalized J: `SetHorizontalRate` does not work in this case.

```
EcoMug genPlane;
genPlane.SetUseSky();
genPlane.SetSkySize({{200.*EMUnits::cm, 200.*EMUnits::cm}});
genPlane.SetSkyCenterPosition({0., 0., 1.*EMUnits::mm});

cout << "Estimated time [s] = " << genPlane.GetEstimatedTime(10000) << endl;
```



# Deal with background

The class `EMMultiGen` allows to handle the generation of the background as well as the signal. It requires a `EcoMug` instance for the signal and one or more instances for the background. Additionally the user has to specify the differential flux (even unnormalized), the PID ([Monte Carlo particle numbering scheme](https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)) and the relative weight (w.r.t. signal) for all backgrounds.

```
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
```

