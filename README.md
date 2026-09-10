# EcoMug: Efficient COsmic MUon Generator

EcoMug is a header-only C++17 library for the generation of cosmic ray (CR) muons, based on a parametrization of experimental data. Unlike other tools, EcoMug gives the possibility of generating from different surfaces (plane, cylinder, half-sphere, and a sphere enclosing the detector), while keeping the correct angular and momentum distribution of generated tracks. EcoMug also allows the generation of CR muons according to user-defined parametrizations of their differential flux.

If you use, or want to refer to, EcoMug please cite the following paper:

> Pagano, D., Bonomi, G., Donzella, A., Zenoni, A., Zumerle, G., & Zurlo, N. (2021). EcoMug: An Efficient COsmic MUon Generator for cosmic-ray muon applications. Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, 1014, 165732.

Latest release: [EcoMug v3.0](https://github.com/dr4kan/EcoMug/releases/tag/v3.0)



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


### Target-sphere generation

The three generation surfaces above cover the apparatus and then rely on the user
to discard the muons that miss it. For a typical detector that is almost all of
them: a 20x20 cm telescope with 50 cm between the planes keeps fewer than one
generated muon in a thousand, and the waste grows with the square of the
generation surface, so enlarging the sky to avoid losing inclined muons makes it
worse.

`SetUseTargetSphere()` instead aims the generation at a sphere enclosing the
detector, so every muon is produced already pointing through it:

```
EcoMug gen; // initialization of the class
gen.SetUseTargetSphere(); // generation aimed at the detector
gen.SetTargetSphereRadius(0.29); // radius of a sphere enclosing the detector
// (x,y,z) position of the center of that sphere
gen.SetTargetSphereCenterPosition({{0., 0., 0.25}});

for (auto event = 0; event < number_of_events; ++event) {
  gen.Generate();
  std::array<double, 3> muon_position = gen.GetGenerationPosition();
  double muon_p = gen.GetGenerationMomentum();
  ...
}
```

The muon starts on the surface of the sphere, so it is already outside the
detector volume and can be handed straight to a transport code.

This is an exact construction, not an approximation: for a given direction the
muons crossing a sphere of radius R are exactly those crossing the disc of radius
R perpendicular to that direction. Rates and `GetEstimatedTime()` are normalised
to the disc area pi*R^2 and stay directly comparable with the other geometries --
the same detector sees the same rate, it just takes far fewer generated muons to
measure it. For the telescope above, against a plane covering the same solid
angle:

| generation surface | generated muons | time |
| ------------------ | --------------- | ---- |
| sky 2 x 2 m        | 4.3 M           | 3.4 s |
| sky 8 x 8 m        | 70.5 M          | 54.5 s |
| sky 16 x 16 m      | 278.8 M         | 216.7 s |
| target sphere, R = 0.29 m | 0.4 M    | 0.4 s |

Note that the gain depends on how well a sphere encloses the apparatus: a
compact detector benefits most, a long thin one least.


# More Advanced Usage

It is possible to set the seed in EcoMug, for reproducible generations. This can be done with the method `SetSeed`, as shown in the example below. If the seed is set to 0 (or the method is not invoked at all), a random seed is used.

One seed drives everything a generator produces, the muon charge included, so two runs
with the same seed agree event by event. This matters when EcoMug feeds a transport code:
a muon of the opposite charge makes Geant4 consume a different number of random numbers,
and everything after it diverges. (Before v3.0 the charge came from a separate engine
that `SetSeed` did not reach.) `EMMultiGen` has its own `SetSeed`, which seeds the source
selection and every contained generator, each with an independent stream.

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

With a user-supplied flux the normalisation is not known to EcoMug, so rates cannot be
put on an absolute scale: `SetHorizontalRate` has no effect, and `GetEstimatedTime`
returns 0 and prints a warning. `GetAverageGenRate` still returns the integral of the
flux you supplied, in whatever units that flux is expressed in, so normalise J yourself
if you need absolute rates.

```
EcoMug genPlane;
genPlane.SetUseSky();
genPlane.SetSkySize({{200.*EMUnits::cm, 200.*EMUnits::cm}});
genPlane.SetSkyCenterPosition({0., 0., 1.*EMUnits::mm});

cout << "Estimated time [s] = " << genPlane.GetEstimatedTime(10000) << endl;
```



# Deal with background

The class `EMMultiGen` allows to handle the generation of the background as well as the signal. It requires a `EcoMug` instance for the signal and one or more instances for the background. Additionally the user has to specify the differential flux (even unnormalized), the PID ([Monte Carlo particle numbering scheme](https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)) and the relative weight (w.r.t. signal) for all backgrounds.

The weights are relative and the signal implicitly carries weight 1, so a component fires
with probability `w_i/(1 + sum_j w_j)`. A PID of 0, which the signal keeps, means "muon,
charge chosen by the generator", and `GetPID()` then returns 13 or -13. Each source keeps
its own configuration, its differential flux included: a source set up with
`SetDifferentialFlux` is driven with `GenerateFromCustomJ`, the others with `Generate`.
`EMMultiGen::SetSeed` makes the whole mixture, source selection included, reproducible.

```
EcoMug muonGen;
muonGen.SetUseSky();
muonGen.SetSkySize({{200.*EMUnits::cm, 200.*EMUnits::cm}});
muonGen.SetSkyCenterPosition({0., 0., 1.*EMUnits::mm});

EcoMug electronGen(muonGen);
electronGen.SetDifferentialFlux(&J);

EcoMug positronsGen(muonGen);
positronsGen.SetDifferentialFlux(&J);

EMMultiGen genSuite(muonGen, {electronGen, positronsGen});
genSuite.SetBckWeights({0.2, 0.1});
genSuite.SetBckPID({11, -11});

map<int, int> counts;
for (auto i = 0; i < number_of_events; ++i) {
  genSuite.Generate();
  counts[genSuite.GetPID()]++;
}
```

In case you want to compile the previous code, please take a look at the following example, which should be compiled with the -std=c++11 flag.

```
#include <iostream>
#include <map>
#include <iomanip>
#include "EcoMug.h"

double J(double p, double theta) {
  double A = 1400*pow(p, -2.7);
  double B = 1. / (1. + 1.1*p*cos(theta)/115.);
  double C = 0.054 / (1. + 1.1*p*cos(theta)/850.);
  return A*(B+C);
};

int main() {

    // EcoMug logs at WARNING by default; raise it to silence warnings, or
    // lower it to see more. Assign, do not redefine: the header already
    // defines this member.
    EMLog::ReportingLevel = EMLog::WARNING;


    EcoMug muonGen;
    muonGen.SetUseSky();
    muonGen.SetSkySize({{200.*EMUnits::cm, 200.*EMUnits::cm}});
    muonGen.SetSkyCenterPosition({0., 0., 1.*EMUnits::mm});

    EcoMug electronGen(muonGen);
    electronGen.SetDifferentialFlux(&J);

    EcoMug positronsGen(muonGen);
    positronsGen.SetDifferentialFlux(&J);

    EMMultiGen genSuite(muonGen, {electronGen, positronsGen});
    genSuite.SetBckWeights({0.2, 0.1});
    genSuite.SetBckPID({11, -11});

    std::map<int, int> counts;
    for (auto i = 0; i < 10000; ++i) {
        genSuite.Generate();
        counts[genSuite.GetPID()]++;
    }
    std::cout << std::right << std::setw(5) << "PID" << std::right << std::setw(8) << " counts" 
        << "     ratio" << std::endl; 

    for (auto const& x : counts) {
        std::cout << std::right << std::setw(5) << x.first << std::right << std::setw(8) << x.second 
            << "   (" << std::setprecision(3) << (double) x.second/(counts[13]+counts[-13]) << ")" << std::endl;
    }
}
```


# Tests and examples

Two ROOT macros ship with the library.

`EcoMugExample.C` is a set of worked examples, meant to be read and copied from:

```
root -l -b -q 'EcoMugExample.C+(1)'        # a detector in a muon flux: rate and live time
root -l -b -q 'EcoMugExample.C+(2)'        # a two-plane telescope, the slow way and the fast way
root -l -b -q 'EcoMugExample.C+(3)'        # handing muons to Geant4
root -l -b -q 'EcoMugExample.C+(4)'        # supplying your own differential flux
root -l -b -q 'EcoMugExample.C+(5)'        # mixing several particle sources with EMMultiGen
```

`EcoMugTests.C` is the test suite. Every check has an explicit pass/fail criterion,
and the macro writes `EcoMugTests.pdf` with the momentum, angular, position and
charge distributions overlaid on their analytic expectations:

```
root -l -b -q 'EcoMugTests.C+'             # default statistics
root -l -b -q 'EcoMugTests.C+(500000)'     # more statistics
```

The trailing `+` compiles with ACLiC; interpreted, both are far slower.
