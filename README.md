# ezQMMM2.0 — Easy QM/MM Input File Generator

Generates QM/MM single-point calculation input files from CHARMM/NAMD MD
trajectories. Supports **ORCA**, **Q-Chem**, and **Psi4**.

---

## Requirements

```bash
pip install MDAnalysis numpy pyyaml
```

---

## Quick Start

```bash
# Generate an example config
python ezQMMM2.py --example

# Edit config_example.yaml, then run
python ezQMMM2.py config.yaml
```

---

## Configuration Reference

```yaml
# ── Trajectory inputs ────────────────────────────────────────────────────────
psf_file: system.psf          # CHARMM PSF topology file
dcd_file: trajectory.dcd      # DCD trajectory file

# ── QM region ────────────────────────────────────────────────────────────────
qm_selection: "resid 100"     # MDAnalysis atom selection string
charge: 0                     # Total QM charge
multiplicity: 1               # Spin multiplicity

# ── QM method ────────────────────────────────────────────────────────────────
method: B3LYP                 # DFT functional or wavefunction method
basis: 6-31G*                 # Basis set (use 'gen' for custom, see *_blocks)
program: qchem                # Target program: orca | qchem | psi4

# ── MM environment ───────────────────────────────────────────────────────────
mm_cutoff: 40.0               # Cutoff radius (Angstrom) for MM point charges
mm_switchdist: 35.0           # Start of switching zone (Angstrom)
                              # Default: mm_cutoff - 5.0

# ── Boundary scheme ──────────────────────────────────────────────────────────
boundary_scheme: RCD          # RCD | CS | Z1 | NONE

# ── Frame selection ──────────────────────────────────────────────────────────
first_frame: 0                # First frame index (0-based)
last_frame: 100               # Last frame index (-1 = last frame)
stride: 10                    # Step between frames

# ── Output ───────────────────────────────────────────────────────────────────
output_dir: ./qmmm_calcs      # Output directory (created if absent)
output_prefix: qmmm           # Prefix for all output files

# ── Periodic images (supercell) ──────────────────────────────────────────────
supercell_axes: []            # Axes to tile: any combo of x/y/z or a/b/c
                              # e.g. [x, y] for membrane (z = normal)
                              #      [x, y, z] for fully periodic system
                              # Leave empty to disable

# ── Structure output ─────────────────────────────────────────────────────────
pdb_stride: null              # Write PSF+PDB every N frames
                              # Options: all | half | tenth | integer | null

# ── Program-specific keywords ────────────────────────────────────────────────
qchem_keywords: |
  scf_convergence    8
  mem_total          8000

qchem_blocks: |
  $basis
  C 0
  cc-pVDZ
  ****
  $end
```

---

## Output Files

For each frame the following are generated under `output_dir`:

| File | Description |
|---|---|
| `<prefix>_frame<N>_qchem.in` | Q-Chem input (or `_orca.inp` / `_psi4.dat`) |
| `<prefix>_frame<N>_charges.pc` | ORCA only — external point charge file |
| `<prefix>_frame<N>_struct.pdb` | Full-system PDB with beta-column labels |
| `<prefix>_frame<N>_struct.psf` | Copy of original PSF (paired topology) |
| `<prefix>_boundary.log` | Per-frame boundary charge modification detail |
| `<prefix>_switching.log` | Per-frame switching-zone charge detail |

### Beta column in PDB

| Value | Region |
|---|---|
| 9 | QM atom |
| 8 | MM point charge (within cutoff) |
| 0 | Outside cutoff |

In VMD: **Graphics → Representations → Coloring Method → Beta**

---

## MDAnalysis Selection Syntax

The `qm_selection` field uses MDAnalysis selection language, which is similar
to but **not identical** to VMD's Tcl-based atom selection. The key
differences and common patterns are shown below.

### A Few Selection Examples

```yaml
# Single residue by number
qm_selection: "resid 100"

# Multiple residues
qm_selection: "resid 100 101 102"

# Residue range
qm_selection: "resid 100:105"

# Specific segment and residue (CHARMM segid)
qm_selection: "segid PROA and resid 100"

# Residue by name (e.g. all HEM residues)
qm_selection: "resname HEM"

# Specific atoms within a residue
qm_selection: "resid 100 and name FE"

# Exclude backbone — keep only sidechain + Cα
qm_selection: "resid 100 and not backbone"

# Cofactor plus coordinating residues
qm_selection: "resname FAD or (resid 45 112 and not backbone)"

# All atoms within a distance of a residue (sphere around ligand)
qm_selection: "resname LIG or (around 3.5 resname LIG)"

# Metal ion plus its coordination shell
qm_selection: "resname ZN or (sphzone 2.5 resname ZN)"
```

### Key Differences from VMD

| Feature | MDAnalysis | VMD |
|---|---|---|
| Residue number | `resid 100` | `resid 100` |
| Residue range | `resid 100:105` | `resid 100 to 105` |
| Segment ID | `segid PROA` | `segid PROA` |
| Residue name | `resname HEM` | `resname HEM` |
| Atom name | `name CA` | `name CA` |
| Logical AND | `and` | `and` |
| Logical OR | `or` | `or` |
| Logical NOT | `not` | `not` |
| Distance selection | `around 3.5 resname LIG` | `within 3.5 of resname LIG` |
| Sphere zone | `sphzone 2.5 resname ZN` | `within 2.5 of resname ZN` |
| Backbone atoms | `backbone` | `backbone` |
| Protein atoms | `protein` | `protein` |
| By index (0-based) | `index 0 1 2` | `index 1 2 3` (1-based) |
| Atom type | `type CT1` | `type CT1` |
| Chain | `segid` (CHARMM) | `chain` (PDB) |

The most common pitfall when translating a VMD selection to MDAnalysis is:

- **`within` → `around`**: VMD uses `within X of ...`, MDAnalysis uses `around X ...`
- **Index offset**: VMD atom indices are **1-based**; MDAnalysis indices are **0-based**
- **Chain vs segid**: VMD uses `chain A` for PDB chain IDs; CHARMM PSF files use `segid` (e.g. `PROA`, `MEMB`). For PSF-based systems always use `segid` in MDAnalysis.
- **Residue ranges**: VMD uses `resid 100 to 105`; MDAnalysis uses `resid 100:105`

### Verifying Your Selection

It is strongly recommended to verify the atom count before running a full
trajectory. Use the following Python snippet:

```python
import MDAnalysis as mda
u = mda.Universe("system.psf", "trajectory.dcd")
sel = u.select_atoms("resid 100 and not backbone")
print(f"Selected {len(sel)} atoms")
for atom in sel:
    print(f"  {atom.segid} {atom.resid} {atom.resname} {atom.name} "
          f"mass={atom.mass:.3f} charge={atom.charge:.4f}")
```

Cross-check the atom count and identity against VMD before running ezQMMM2.0.


### MM Point Charge Selection

MM atoms are selected by whole-residue inclusion: if any atom of a residue
falls within `mm_cutoff` of any QM atom, the entire residue is included.
Distance is measured as the minimum distance to any QM atom (not center of
mass), using the minimum image convention (PBC-aware).

### Switching Function

A NAMD-style quintic switching function smoothly attenuates charges between
`mm_switchdist` and `mm_cutoff`:

$S(r) = 1 - 10t^3 + 15t^4 - 6t^5, \qquad t = \frac{r - r_\text{sw}}{r_\text{cut} - r_\text{sw}}$

$S(r) = 1 \quad \text{for } r \leq r_\text{sw}, \qquad S(r) = 0 \quad \text{for } r \geq r_\text{cut}$

Distance $r$ is the minimum distance from the charge to any QM atom. All
switching events are recorded in `<prefix>_switching.log`.

### Boundary Schemes

The QM/MM boundary is where a covalent bond is cut between the QM and MM
regions. A link hydrogen atom is placed along each cut bond at 1.09 Å from
the QM atom. The MM1 atom (directly bonded to QM) is treated according to
the chosen scheme:

| Scheme | MM1 treatment | Purpose |
|---|---|---|
| `RCD` | Removed; charge redistributed to midpoint virtual charges and MM2 atoms adjusted | Preserves bond dipole and charge neutrality |
| `CS` | Removed; charge split into a dipole pair around MM2 | Charge shift scheme |
| `Z1` | Removed (charge zeroed) | Simplest; may break neutrality |
| `NONE` | No modification | Use only if no QM/MM covalent bonds exist |

**RCD virtual charge factor of 2**: the virtual charge placed at the
MM1–MM2 midpoint is `2 * q_MM1 / n` (not `q_MM1 / n`). This preserves the
dipole moment of the original MM1–MM2 bond: `q·d = 2q·(d/2)`. The MM2
charge is simultaneously adjusted by `-q_MM1 / n` to maintain overall
charge neutrality.

All modifications are recorded in `<prefix>_boundary.log`.

### Coordinate Remapping

`distance_array(..., box=dimensions)` selects MM atoms using the minimum
image convention but retains raw trajectory coordinates. Before writing to
any input file, all MM charge positions are remapped to their minimum image
equivalent relative to the QM centroid:

$x_\text{remapped} = x - \text{round}\!\left(\frac{x - x_\text{QM}}{L_x}\right) L_x$

and equivalently for $y$ and $z$. This ensures charge coordinates in the
input files are physically near the QM region and that periodic image
tiling does not double-count atoms.

### Supercell / Periodic Images

In a periodic MD simulation, the QM region interacts not only with MM atoms
in the primary simulation box but also with their periodic images in
neighboring cells. For most large solvated systems the box is sufficiently
large that no image falls within the cutoff, and the standard non-periodic
treatment is adequate. However, for smaller boxes or systems where the QM
region is close to the box boundary, neglecting periodic images introduces
an asymmetry in the electrostatic environment that can affect the QM
wavefunction and energetics.

When `supercell_axes` is set, ezQMMM2.0 tiles the primary MM charge set along
the specified axes to generate explicit image copies. The number of shells
along each axis is derived automatically:

$n_x = \left\lceil \frac{r_\text{cut}}{L_x} \right\rceil, \qquad
  n_y = \left\lceil \frac{r_\text{cut}}{L_y} \right\rceil, \qquad
  n_z = \left\lceil \frac{r_\text{cut}}{L_z} \right\rceil$

This guarantees that every image within $r_\text{cut}$ of the QM region is
evaluated, regardless of box shape. For each shell combination
$(i_x, i_y, i_z)$ excluding $(0, 0, 0)$, the entire primary charge set is
translated by $(i_x L_x,\, i_y L_y,\, i_z L_z)$ and each candidate image
is distance-filtered against all QM atoms. Only images within $r_\text{cut}$
are retained.

**Axis selection** is intentionally flexible. Common use cases are:

| System | Recommended setting |
|---|---|
| Large solvated protein | `supercell_axes: []` — images never within cutoff |
| Membrane protein (normal along $z$) | `supercell_axes: [x, y]` — periodic in bilayer plane only |
| Small solvated active site model | `supercell_axes: [x, y, z]` — full 3D tiling |
| Elongated system (e.g. DNA fiber, fibril) | `supercell_axes: [z]` — periodic along fiber axis only |

**Coordinate remapping before tiling** is essential for correctness. The
MDAnalysis `distance_array` call that selects primary MM atoms uses the
minimum image convention internally, which finds atoms correctly across PBC
but retains their raw trajectory coordinates. An atom physically close to
QM via PBC may have a raw coordinate a full box length away. Before tiling,
all primary charges are remapped to their minimum image position relative to
the QM centroid:

$x_\text{remapped} = x - \text{round}\!\left(\frac{x - x_\text{QM}}{L_x}\right) L_x$

Without this remapping, adding $\pm L_x$ to an already-displaced coordinate
would place image copies at double the box length rather than at the correct
neighboring cell position, producing ghost charge shells in the input files.

**Switching is applied after tiling**, so image charges are attenuated by
the same quintic switching function as primary charges, using their explicit
Cartesian distance to the nearest QM atom. Images that survive the cutoff
filter but fall in the switching zone will have their charges smoothly
reduced to zero at $r_\text{cut}$, ensuring no discontinuity at the
boundary of the image shell.

**For large boxes `img=0` is expected and correct.** After remapping, every
primary charge is within $r_\text{cut}$ of QM by construction. The nearest
possible image of any such charge is at approximately $L - r_\text{cut}$
from QM. For a box of $L = 141\,\mathring{\text{A}}$ and $r_\text{cut} = 60\,\mathring{\text{A}}$ this gives a
minimum image distance of ${\sim}81\,\mathring{\text{A}}$ — well outside the cutoff. Images
will only contribute when $L \lesssim 2\, r_\text{cut}$.

---

## Caveats and Limitations

### QM Region
- **Must not straddle a periodic boundary.** The remapping formula uses the
  mean QM position as a reference. If QM atoms are split across the box, the
  centroid will be meaningless and all charge remapping will be wrong. Wrap
  the QM region before use.
- **Element assignment is mass-based.** Elements are inferred from atomic
  masses because PSF files store CHARMM atom types, not element symbols. The
  lookup table covers H, D, C, N, O, Na, Mg, P, S, Cl, K, Ca, Fe, Cu, Zn.
  Unusual elements will be assigned `X` and the QM input will be wrong.
  Check `QM=N` in the console output against expectation.

### MM Environment
- **Whole-residue inclusion** means a single atom near the cutoff edge pulls
  in its entire residue. For large residues (lipids, polymers) this can
  significantly increase the MM region size.
- **Charges are PSF partial charges.** No polarization is applied. The MM
  environment is purely electrostatic.
- **No MM geometry optimization** is performed. Input files are for single-
  point energy calculations only.

### Boundary
- **Bonds must be in the PSF.** `_find_boundary_bonds` relies on
  `atom.bonds` from MDAnalysis. If the PSF was built without explicit bonds
  (rare but possible), no boundary bonds will be detected and no link atoms
  or charge corrections will be applied regardless of the chosen scheme.
- **Only C–C type cuts are well-tested.** Cutting across polar bonds (C–N,
  C–O) with a link hydrogen introduces additional errors that none of the
  boundary schemes fully correct. Avoid such cuts where possible.

### Supercell
- **Orthorhombic boxes only.** The image tiling uses simple `n * L` offsets
  along Cartesian axes. Triclinic boxes require lattice vector arithmetic and
  are not supported.
- **For large boxes `img=0` is expected and correct.** Images are only
  non-zero when `L < mm_cutoff + max_QM_to_edge_distance`. For boxes much
  larger than the cutoff, no periodic neighbours will ever be within range.
- **`supercell_axes` disables PBC for the primary selection.** When tiling
  is active, primary charges are selected using raw Cartesian distances (no
  PBC wrap). PBC is handled exclusively through explicit image copies. This
  is intentional and correct but means the cutoff sphere may appear flat-
  faced for very asymmetric QM placements. If in doubt, visualize the PDB.

### PDB / PSF Output
- **Long bonds in VMD are a visualization artifact.** Atoms near the box
  edge that are connected via PBC will appear connected by a line crossing
  the entire box. This does not affect the QM/MM input files.
- **The copied PSF contains all system atoms.** The paired PDB also contains
  all atoms. This is intentional — VMD uses the PSF for connectivity and the
  PDB for coordinates. The beta column identifies which atoms are active.
- **tempfactors are not in PSF files.** MDAnalysis initializes them to zero
  at load time. The note printed at startup is informational only.
- **Enabling `pdb_stride` significantly slows down the run.** Writing a
  full-system PDB requires MDAnalysis to format and write coordinates for
  every atom in the system (e.g. 436k atoms) for each requested frame. For
  large systems this can take several seconds per frame and dominate the
  total runtime. Use `pdb_stride: tenth` or a large integer during
  production runs and reserve `pdb_stride: all` for quick validation of a
  small number of frames. If you only need to verify the QM/MM partitioning,
  running a single frame with `pdb_stride: all` is sufficient.

### Program-Specific
- **ORCA**: point charges are written to a separate `_charges.pc` file.
  Keep the `.inp` and `.pc` files in the same directory when running ORCA.
- **Psi4**: coordinates in the `molecule` block are in Angstrom (Psi4
  default); point charges are converted to Bohr (`× 1.88973`).
- **Q-Chem**: charges appear inline in the `$external_charges` block in
  `x y z q` order (Angstrom).
- **Custom basis sets**: set `basis: gen` and provide the basis block in
  `qchem_blocks` / `orca_blocks` / `psi4_blocks`. For ORCA, the custom
  block goes in the `%basis` section.

---

## References

- NAMD QM/MM: https://www.ks.uiuc.edu/Research/qmmm/
- Melo et al., *Nature Methods* **15**, 351–354 (2018)
- RCD and RC boundary schemes: Lin & Truhlar, *J. Phys. Chem. A* **109**, 3991–4004 (2005). DOI: 10.1021/jp0446332
- MDAnalysis: Michaud-Agrawal et al., *J. Comput. Chem.* **32**, 2319 (2011)
