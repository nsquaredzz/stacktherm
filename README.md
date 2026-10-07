# stacktherm

A 3D thermal simulator for chip-stack unit cells that computes its own material
inputs. You describe a slice of a stacked chip (transistor tier, wiring levels,
bond layer); it returns the 3D temperature field and the effective anisotropic
conductivity and thermal resistance of every layer, in formats a package-scale
tool can load.

Three things that bulk-property models leave out are built in:

- **Size effects.** A 16 nm copper line, a 7 nm silicon fin and a 2 µm bond pad
  get different conductivities, computed from the feature's own dimensions.
- **Interface resistance.** Every material boundary carries a thermal boundary
  resistance, with its provenance (measured, user, assumed) written to the output.
- **Glass conductivity from atoms.** The conductivity of amorphous SiO2 is
  computed from an atomistic model with the Wigner formulation of heat
  transport, rather than looked up, and follows the run temperature.

The numerical core is C++17: finite-volume solver with algebraic multigrid,
stack voxeliser, material models, glass molecular dynamics and the Wigner
engine, with a command-line program. Python is used only for the browser studio
and for the bridge to kALDo.

## Build

Requirements: a C++17 compiler, CMake 3.20+, BLAS/LAPACK (Accelerate on macOS,
OpenBLAS elsewhere) and, for the studio, Python 3.10+. AMGCL and yaml-cpp are
fetched by CMake.

```bash
python3 -m venv .venv && .venv/bin/pip install pyyaml numpy pytest pybind11 cmake ninja
```

```bash
PATH="$PWD/.venv/bin:$PATH" cmake -S . -B build -G Ninja -DPython_EXECUTABLE="$PWD/.venv/bin/python"
```

```bash
PATH="$PWD/.venv/bin:$PATH" cmake --build build
```

This produces `build/stacktherm` (the command line), the test programs, and the
Python extension `stacktherm/_core*.so`. Pass `-DSTACKTHERM_PYTHON=OFF` to build
the C++ parts alone.

```bash
PATH="$PWD/.venv/bin:$PATH" ctest --test-dir build --output-on-failure
```

```bash
.venv/bin/python -m pytest
```

## Quick start

```bash
build/stacktherm init my_stack.yaml
```

```bash
build/stacktherm run my_stack.yaml -o out --report --budget
```

```bash
build/stacktherm validate
```

```bash
.venv/bin/python -m stacktherm serve
```

| command | what it does |
|---|---|
| `stacktherm init` | writes a commented starter stack file |
| `stacktherm run` | solves a stack file; `-o DIR`, `--report`, `--budget`, `--converge`, `--temps`, `--refine` |
| `stacktherm validate` | closed-form and published-data checks |
| `stacktherm wigner` | glass conductivity from atoms with the native engine; `--recompute`, `--model` |
| `stacktherm glass` | melts, quenches and relaxes a new amorphous SiO2 sample |
| `python -m stacktherm serve` | interactive studio in the browser |
| `python -m stacktherm kaldo` | the same glass through kALDo, or `--import` external force constants |

`run -o DIR` writes:

| file | contents |
|---|---|
| `layers.csv` | kx, ky, kz, resistance and heat capacity of every layer |
| `layers_apdl.mac` | Ansys APDL orthotropic material cards (`MP,KXX/KYY/KZZ`), temperature tables with `--temps` |
| `result.json` | everything: layers, computed material inputs, every interface resistance and where it came from, solver record |
| `field.vtr` | temperature, conductivity, heat source and material id on the mesh (ParaView / VisIt) |
| `report.html` | resistance budget, interactive temperature slices, tables (`--report`) |

## Interactive studio

`python -m stacktherm serve` opens a local page (http://127.0.0.1:8770) backed by
the C++ core:

- **3D cell.** Four parametric cells: wiring levels, hybrid bond, microbump,
  transistor tier. Moving a slider (linewidth, metal, via placement, pad
  misalignment, bond conductance, underfill, ...) re-solves the 3D field in a
  fraction of a second. Rotate it, cut it open, colour by temperature, heat flux
  or conductivity, and switch on heat-flow paths. Every solve is paired with a
  bulk-property solve of the same geometry, so the cost of ignoring size effects
  and interfaces is always on screen.
- **Your stack.** Paste or edit a stack file, run it through the full solver
  (sub-stacks included), view the result in 3D and download the layer table,
  APDL cards, JSON and report. Tick "error budget" to rank the inputs the answer
  depends on.
- **Materials engine.** Conductivity against feature size for Cu, Ru, Co, W and
  silicon, recomputed by the core on every change, with the published silicon
  measurements overlaid.
- **Wigner engine.** The atomistic calculation described below, native and kALDo.

The page loads three.js from a CDN, so the 3D tab needs a network connection.

## Describing a stack

A stack is a bottom-to-top list of layers in a laterally periodic unit cell
(see `examples/`). Lengths accept units (`42 nm`, `1 um`), power `100 W/cm2`,
interfaces either as a conductance `150 MW/m2K` or a resistance `20 m2K/GW`.

| type | meaning | key fields |
|---|---|---|
| `film` | uniform material | `material`, `thickness` |
| `lines` | parallel lines in a matrix (wiring levels, fins) | `material`, `matrix`, `width`, `pitch`, `direction`, `barrier` |
| `posts` | 2D array of square or round posts (vias, bond pads, bumps, TSVs) | `material`, `matrix`, `size` or `diameter`, `pitch`, `offset`, `barrier` |
| `interface` | extra resistance on a plane, per material pair (bond interfaces) | `pairs`, `R` / `tbc` |
| `substack` | a group of fine-pitch layers solved in its own small cell | `cell`, `layers`, `flip`, `lump` |
| `aniso` | a block with given kx, ky, kz | `kx`, `ky`, `kz` |

Any layer can carry `heat: {flux: 50 W/cm2, region: [x0, x1, y0, y1]}`.

**Sub-stacks are how the scales are bridged.** A bond pad is a micrometre wide
and an M1 line 16 nm; resolving both in one mesh would need billions of cells.
A `substack` is characterised in its own nanometre-scale periodic cell and
enters the parent cell as anisotropic layers that preserve its through-plane
resistance and lateral conductivities exactly. Identical definitions (both
tiers, via YAML anchors) are solved once.

## How it works

**Material inputs** (`cpp/src/materials.cpp`)

- *Metals* (Cu, W, Co, Ru, Al, Mo): Fuchs-Sondheimer surface scattering plus
  Mayadas-Shatzkes grain-boundary scattering give the resistivity of a wire or
  film; Wiedemann-Franz converts it to thermal conductivity with a Lorenz number
  pinned to the bulk metal. Barriers narrow the conducting core and add series
  resistance on the walls.
- *Silicon*: kinetic-theory integral over a Born-von Karman dispersion with
  isotope and Umklapp scattering. The isotope term is computed, the two Umklapp
  parameters are fitted to the bulk k(T) curve (within 2 % from 150 to 600 K),
  and films, fins and nanosheets are then predictions: exact Fuchs-Sondheimer
  suppression in-plane, a ballistic-diffusive interpolation cross-plane, a
  Casimir length for wires.
- *Dielectrics, solder, underfill*: tabulated in `data/library.yaml`, each
  with its source. Porous variants scale as (1 - porosity)^1.5.
- *Interfaces*: measured boundary conductances where the library has one, a
  stated default otherwise. A diffuse-mismatch estimate is available but off by
  default (see validation).

**Solver** (`cpp/src/fvm.cpp`, `cpp/src/solver.cpp`): cell-centred finite volumes on a
rectilinear mesh that conforms to every feature edge. Each face conductance is
the series sum of two half-cells and the interface resistance on that face, so
material jumps and boundary resistance are exact for axis-aligned geometry.
Conjugate gradients with classical algebraic multigrid (AMGCL); cost is linear in
the number of cells.

**What is solved** (`cpp/src/characterize.cpp`): one unit-temperature-drop solve splits
the total through-plane resistance into a term per layer and per interface,
using plane-averaged temperatures (so a 1D stack of the exported layers
reproduces the average temperature of every plane). Two periodic solves give the
in-situ lateral conductivity of every layer. A fourth solve applies your
boundary conditions and heat sources for the temperature field.

The table also lists `kz mix`, the rule-of-mixtures value of the layer on its
own. Where it is far above the in-situ kz (bond pads on top of low-conductivity
wiring), heat cannot reach the conductor without spreading through its
neighbours; that penalty is invisible to a per-layer mixing rule.

## Error budget

`stacktherm run --budget` (or the tick box in the studio) raises each class of
input by 10 % in turn, re-solves the stack and reports how much the total
resistance moves:

```
  input class                                            elasticity assumed +/-  effect
  conductivities marked as estimates in the library           -0.20        30 %   6.1 %
  silicon phonon model                                        -0.27        10 %   2.7 %
  bond interfaces and barriers given in the stack file        +0.08        30 %   2.4 %
```

(microbump example). The elasticities are computed; the uncertainty assigned to
each class is a stated assumption, so read the output as a ranking of what to
measure or pin down first. For the microbump stack that is the underfill and
solder conductivities; for the hybrid-bond stack it is the bond-interface
values in the stack file.

## Wigner engine: dielectric conductivity from atoms

The wiring levels of a stack are limited by their amorphous dielectrics, and a
glass is exactly where the phonon-gas picture behind most thermal models fails:
its vibrations do not travel, they hand energy to one another. The Wigner
formulation of heat transport (Simoncelli, Marzari and Mauri, 2019) covers this
regime and the crystalline one with a single expression. It is used in
first-principles materials codes for bulk conductivity; here it supplies a
material input to a device-scale solver.

```bash
build/stacktherm wigner
```

What `cpp/src/glass.cpp` and `cpp/src/wigner.cpp` do:

1. Melt and quench a 648-atom SiO2 cell by molecular dynamics (BKS pair
   potential with truncated Coulomb), then relax it until the forces vanish.
2. Build the exact Hessian, diagonalise it, and form the velocity operator
   between every pair of the 1941 vibrational modes.
3. Get cubic force constants by finite differences of the Hessian along sampled
   modes, and from them three-phonon linewidths at each temperature.
4. Evaluate the Wigner conductivity, in which heat is carried by tunnelling
   between modes whose frequencies overlap within their linewidths.

Result for the first stored model (W/m/K):

| T (K) | native Wigner | kALDo | harmonic limit | reference |
|---|---|---|---|---|
| 100 | 0.47 | 0.42 | 0.43 | 0.69 (handbook, recalled) |
| 200 | 0.99 | 0.96 | 0.88 | 1.14 (handbook, recalled) |
| 300 | **1.34** | **1.32** | 1.17 | **1.38** (measured) |
| 400 | 1.55 | 1.53 | 1.34 | 1.51 (handbook, recalled) |
| 600 | 1.77 | 1.74 | 1.51 | 1.75 (handbook, recalled) |

No transport parameter is fitted. At 300 K the native result is 3 % below the
measured value and kALDo's is 5 % below. The native result moves by under 3 %
when the numerical regularisation is varied sixteen-fold, and not at all with
the number of sampled modes or the width used for energy conservation. The
harmonic (Allen-Feldman) limit, without linewidths, does depend on the
regularisation at this cell size (1.00 to 1.33 over the same range); the
anharmonic linewidths are what make the calculation well defined.

Use it in a stack as `material: SiO2_wigner`. `stacktherm wigner --recompute`
redoes steps 2 to 4 in about 25 seconds; `stacktherm glass --name NAME --seed N`
builds a new sample.

**Two backends, one result.** The same force constants are also run through
[kALDo](https://github.com/nanotheorygroup/kaldo), an established open-source
package for thermal transport in crystalline and disordered solids, using its
quasi-harmonic Green-Kubo method (`python -m stacktherm kaldo`, about ten
minutes; needs `pip install kaldo`). The two implementations share nothing but the force
constants:

| at 300 K | W/m/K |
|---|---|
| measured, bulk vitreous silica | 1.38 |
| native Wigner engine, sample 1 / sample 2 | 1.338 / 1.326 |
| kALDo, sample 1 / sample 2 | 1.316 / 1.301 |
| `SiO2_wigner` as used by the solver (kALDo, two-sample mean) | 1.309 |

Agreement between the two codes within 3 % from 200 K upward is the check on
the implementation; agreement between two independently quenched glasses within
about 1 % is the check on the sample. Below about 150 K the codes part (0.42
against 0.47 at 100 K) because they treat the discrete spectrum of the small
cell differently, and neither is reliable there. The library averages the
stored samples and prefers a kALDo table where one exists.

One trap worth knowing if you script kALDo yourself: with in-memory storage it
keeps the linewidths of the first temperature when the temperature is changed.
The backend clears that cache per temperature and refuses to write a table
whose linewidths do not vary; a first run without this gave a wrong table,
which the comparison with the native engine exposed.

**Bring your own material.** kALDo reads force constants written by LAMMPS,
DFT codes and machine-learned potentials. A material that stacktherm cannot
model itself becomes a solver input with

```bash
.venv/bin/python -m stacktherm kaldo --import path/to/force_constants --format lammps --name my-lowk
```

and then, in a stack file, `materials: {my-lowk: {table: my-lowk}}`. The native
engine remains for pair-potential glasses, where it is about ten times faster
because it samples the linewidths instead of computing every mode.

What this result does and does not show:

- Only the 300 K reference was checked against a source for this work. The
  other reference values are recalled from a standard handbook table and are
  shown for the trend only.
- The model is too low at 100 K (native -32 %, kALDo -41 %). The cell has no
  vibrations below 33 cm^-1, so long-wavelength sound waves, which carry more
  of the heat when it is cold, are missing.
- BKS is a classical potential, and the glass was quenched at about 1e14 K/s,
  far faster than any real glass. Agreement within 3 % at one temperature for
  one potential should not be read as 3 % accuracy in general.
- It is bulk glass. Thickness dependence, porosity and the low-k SiCOH family,
  where a tabulated number is least trustworthy, are not covered yet. Those
  need structures that are process dependent, which is where validation gets
  hard again.

## How this relates to existing tools

Each link in the chain exists somewhere in a stronger form. imec's in-house
BTE-FEM solves the Boltzmann equation directly for wiring stacks and is
validated on real nodes; RPI and IBM have published a homogenisation workflow
up to die scale; kALDo and phono3py compute glass and crystal conductivity from
atoms; Ansys, Cadence and Siemens cover package and die scale with effective
properties as inputs. What I did not find in a search is a public tool that
connects them: atoms to material inputs, size effects and interface resistance
with provenance, a 3D unit cell, and a layer table a package tool can load.
That chain, with its validation and error budget, is what this package is.

## Validation

`stacktherm validate` reruns all of this. Current results:

**Against closed-form solutions**

| check | result |
|---|---|
| layered stack with interface resistance, laminate (Wiener) bounds | exact to solver tolerance |
| manufactured solution, anisotropic, non-uniform mesh | second order (observed 1.99) |
| square array of cylinders vs Rayleigh series, conductor in insulator | 0.15 % at 128 cells |
| same, insulator in conductor | 0.8 % at 128 cells |
| cylinder with boundary resistance vs the exact equivalent-inclusion result | within 1.5 % |

**Against published measurements** (`data/published.yaml` records
each number and whether it was stated in the text or read from a plot)

| data | result |
|---|---|
| In-plane k of Si membranes, 15 nm to 1.5 µm, 11 points (Cuffe et al., PRB 2015) | mean -7 %, worst -15 %, no fitted parameters |
| 9 nm Si membrane (Chavez-Angel et al., 2014): 9 ± 2 W/m/K | model 13.8, outside the error bar |
| SOI films 74-240 nm (Ju & Goodson, 1999) | model 24-39 % low; these early data sit above later membrane measurements |
| Cross-plane k of M1-M5 in a 22 nm node chip, TDTR (Huang et al., arXiv:2609.20379) | 4 of 5 layers inside the error bars; a bulk-property model without size effects or boundary resistance lands 27-73 % above the measurements |
| Interface conductance: diffuse mismatch model vs measurement | agrees with the Si/SiO2 lower limit, overestimates Cu/SiO2 about tenfold |
| Amorphous SiO2 at 300 K, native Wigner engine: 1.34 W/m/K | 3 % below the measured 1.38, nothing fitted |
| Same glass through kALDo (independent code): 1.32 W/m/K | within 2 % of the native engine, 5 % below measured |
| Second, independently quenched glass: 1.33 native, 1.30 kALDo | within about 1 % of the first |

Read these with their limits in mind:

- The wiring-layer comparison is a consistency check, not a blind prediction. It
  uses the dielectric conductivity and Cu/dielectric boundary resistance that the
  authors extracted from the same measurements, and the line height, level
  height and linewidth are my estimates from their TEM image (the paper gives
  only Cu fractions). The paper's larger set of 40+ points is published as a
  plot without the per-point geometry needed to simulate each one.
- The copper size-effect model is checked against its limits and against the
  generally reported resistivity range, not against a specific measured dataset.
  Its grain-boundary reflection coefficient (R = 0.3) is an assumed typical value.
- The imported-force-constant route (`kaldo --import`) goes through kALDo's
  own readers and has not been exercised here with a real LAMMPS or DFT data
  set; only the built-in glass has been run end to end.
- No measured hybrid-bond versus microbump stack resistance with a fully
  specified geometry was found; the two examples are a like-for-like comparison
  by this tool, not a reproduction of an experiment.
- Bonded SiO2-SiO2 uses the published lower limit of 150 MW/m²K (Zhang et al.,
  arXiv:2601.03106), so it is a conservative bound. The Cu-Cu bond value in the
  examples is assumed.
- SiN, SiCN, underfill and solder conductivities are labelled estimates in the
  library and should be replaced with process data.

## Performance

`examples/hybrid_bond.yaml` on an Apple M4, single thread, including six
sub-stack characterisations and four solves on the main cell:

| cells | time | peak memory | total resistance |
|---|---|---|---|
| 278,400 | 2.8 s | 0.35 GB | 0.8430 mm²K/W |
| 1,171,456 | 10.8 s | 1.5 GB | 0.8402 mm²K/W |
| 3,796,992 | 38 s | 3.8 GB | 0.8390 mm²K/W |

The Wigner calculation for a 648-atom glass (modes, 160 sampled linewidths, 11
temperatures) takes about 25 s.

The core was first written and validated in Python and then ported. The C++
version reproduces that reference to solver tolerance (layer tables to 1e-9,
Wigner conductivities to 1e-10) and is a little over twice as fast; the
reference values are frozen into `cpp/tests`.

## Layout

```
cpp/include/stacktherm/   public header
cpp/src/                  solver, materials, stack, glass, Wigner engine, exports, validation
cpp/cli/                  the stacktherm command
cpp/python/               pybind11 bindings
cpp/tests/                C++ tests (ctest)
stacktherm/               Python shell: studio server and page, kALDo bridge
data/                     material library, published measurements, glass models and tables
examples/                 stack files
tests/                    tests of the Python shell
```

## Limits of version one

- Heat conduction only, steady state. Heat capacity is exported for transient
  package models, but there is no transient solve here.
- Geometry is extruded per layer and axis-aligned. Round posts are a staircase
  with an area-exact axial conductivity and a perimeter-corrected wall
  resistance.
- Size effects enter through each feature's conductivity; transport between
  features is diffusive. Fully ballistic effects between neighbouring
  nanostructures are not resolved.
- A homogenised sub-stack meets its neighbours with no extra boundary
  resistance, and its layers are characterised between isothermal planes.
- Materials are evaluated at one temperature per run (`temperature:`); use
  `--temps` for tables. There is no temperature feedback within a run.

## Licence

MIT, see [LICENSE](LICENSE). The dependencies fetched at build time keep their own
licences: AMGCL and yaml-cpp are both MIT.

## References

- M. Simoncelli, N. Marzari, F. Mauri, *Unified theory of thermal transport in crystals and glasses*, Nat. Phys. 15, 809 (2019); *Wigner formulation of thermal transport in solids*, Phys. Rev. X 12, 041011 (2022); M. Simoncelli, F. Mauri, N. Marzari, npj Comput. Mater. 9, 106 (2023).
- B. W. H. van Beest, G. J. Kramer, R. A. van Santen, Phys. Rev. Lett. 64, 1955 (1990); A. Carre et al., J. Chem. Phys. 127, 114512 (2007). P. B. Allen and J. L. Feldman, Phys. Rev. B 48, 12581 (1993).
- B. P. Barua, M. R. I. Udoy, A. Aziz, *A Review of Multiscale Thermal Modeling in Heterogeneous 3D ICs*, arXiv:2604.03290 (2026), which describes size effects, interfacial resistance and anisotropic conductivity in stacked chips as "first-order terms rather than small corrections."
- Z. Huang, Y. Sun, T. Jia, R. Wang, Z. Cheng, *Predictive Structure-to-Thermal Conductivity Modeling Framework for BEOL Interconnect Stacks…*, arXiv:2609.20379 (2026).
- X. Zhang, L. Chang, L. Li, Z. Cheng, *Thermal conductance across bonded SiOx-SiOx interfaces in hybrid bonding process*, arXiv:2601.03106 (2026).
- J. Kimling et al., Phys. Rev. B 95, 184305 (2017). J. Cuffe et al., Phys. Rev. B 91, 245423 (2015). E. Chavez-Angel et al., APL Mater. 2, 012113 (2014). Y. S. Ju and K. E. Goodson, Appl. Phys. Lett. 74, 3005 (1999).
- R. L. Graham et al., Appl. Phys. Lett. 96, 042116 (2010). D. Gall, J. Appl. Phys. 119, 085101 (2016). J. Wakil, E. G. Colgan, S. Chen, IEEE Trans. CPMT 1, 1007 (2011).
- Model references are in the module docstrings.
