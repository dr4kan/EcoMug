# EcoMug: Efficient COsmic MUon Generator

EcoMug is a header-only C++11 library for the generation of cosmic muons.

# Examples
## Plane-based generation

```
EcoMug gen;
gen.SetUseSky();
gen.SetSkySize({{10., 10., 10.}});

for (auto event = 0; event < number_of_events; ++event) {
  gen.Generate();
  std::array<double, 3> muon_position = gen.GetGenerationPosition();
  double muon_p = gen.GetGenerationMomentum();
  double muon_theta = gen.GetGenerationTheta();
  double muon_phi = gen.GetGenerationPhi();
  ...
}
```
