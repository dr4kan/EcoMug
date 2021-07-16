# EcoMug: Efficient COsmic MUon Generator

EcoMug is a header-only C++11 library for the generation of cosmic muons.

Latest release: [EcoMug v1.3](https://github.com/dr4kan/EcoMug/releases/tag/v1.3)

# Examples
## Plane-based generation

```
EcoMug gen;
gen.SetUseSky();
gen.SetSkySize({{10., 10.}});

for (auto event = 0; event < number_of_events; ++event) {
  gen.Generate();
  std::array<double, 3> muon_position = gen.GetGenerationPosition();
  double muon_p = gen.GetGenerationMomentum();
  double muon_theta = gen.GetGenerationTheta();
  double muon_phi = gen.GetGenerationPhi();
  double muon_charge = gen.GetCharge();
  ...
}
```

## Cylinder-based generation

```
EcoMug gen;
gen.SetUseCylinder();
gen.SetCylinderRadius(10.);
gen.SetCylinderHeight(30.);

for (auto event = 0; event < number_of_events; ++event) {
  gen.Generate();
  std::array<double, 3> muon_position = gen.GetGenerationPosition();
  double muon_p = gen.GetGenerationMomentum();
  double muon_theta = gen.GetGenerationTheta();
  double muon_phi = gen.GetGenerationPhi();
  double muon_charge = gen.GetCharge();
  ...
}
```

## Hsphere-based generation

```
EcoMug gen;
gen.SetUseHSphere();
gen.SetHSphereRadius(30.);

for (auto event = 0; event < number_of_events; ++event) {
  gen.Generate();
  std::array<double, 3> muon_position = gen.GetGenerationPosition();
  double muon_p = gen.GetGenerationMomentum();
  double muon_theta = gen.GetGenerationTheta();
  double muon_phi = gen.GetGenerationPhi();
  double muon_charge = gen.GetCharge();
  ...
}
```
