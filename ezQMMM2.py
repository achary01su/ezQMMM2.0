#!/usr/bin/env python3
"""
ezQMMM 2.0 — Easy QM/MM Input File Generator
Generates QM/MM single point calculation input files from MD trajectories
Supports: ORCA, Q-Chem, and Psi4

Usage:
    python ezQMMM2.py --example
    python ezQMMM2.py config.yaml

Reference:
    NAMD QM/MM: https://www.ks.uiuc.edu/Research/qmmm/
    Paper: Melo et al., Nature Methods 15, 351-354 (2018)
"""

import shutil
import warnings
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import yaml
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple


# ---------------------------------------------------------------------------
# Data records
# ---------------------------------------------------------------------------

class ChargeMod:
    """Immutable record of a single charge modification at the QM/MM boundary."""

    __slots__ = (
        'frame', 'mod_type',
        'atom_index', 'segid', 'resid', 'resname', 'name',
        'psf_charge', 'applied_charge', 'delta',
        'position', 'reason',
    )

    def __init__(self, frame: int, mod_type: str, reason: str,
                 psf_charge: float, applied_charge: float,
                 position: np.ndarray,
                 atom_index: Optional[int] = None,
                 segid: str = '', resid: int = 0,
                 resname: str = '', name: str = ''):
        self.frame          = frame
        self.mod_type       = mod_type
        self.atom_index     = atom_index
        self.segid          = segid
        self.resid          = resid
        self.resname        = resname
        self.name           = name
        self.psf_charge     = psf_charge
        self.applied_charge = applied_charge
        self.delta          = applied_charge - psf_charge
        self.position       = np.array(position)
        self.reason         = reason


class SwitchRecord:
    """Record of a charge scaled by the switching function."""

    __slots__ = ('frame', 'psf_charge', 'scaled_charge', 'scale',
                 'dist', 'position', 'is_image')

    def __init__(self, frame: int, psf_charge: float, scaled_charge: float,
                 scale: float, dist: float, position: np.ndarray,
                 is_image: bool = False):
        self.frame         = frame
        self.psf_charge    = psf_charge
        self.scaled_charge = scaled_charge
        self.scale         = scale
        self.dist          = dist
        self.position      = np.array(position)
        self.is_image      = is_image


# ---------------------------------------------------------------------------
# Main generator
# ---------------------------------------------------------------------------

class QMMMGenerator:
    """Generate QM/MM input files from MD trajectories."""

    MASS_TO_ELEMENT = {
        (0.5,  1.5):  'H',  (1.5,  2.5):  'D',  (11.5, 12.5): 'C',
        (13.5, 14.5): 'N',  (15.5, 16.5): 'O',  (22.5, 23.5): 'Na',
        (24.0, 25.0): 'Mg', (30.5, 31.5): 'P',  (31.5, 32.5): 'S',
        (34.5, 36.0): 'Cl', (38.5, 40.0): 'K',  (39.5, 41.0): 'Ca',
        (54.5, 56.5): 'Fe', (63.0, 64.0): 'Cu', (64.0, 66.0): 'Zn',
    }

    def __init__(self, psf_file: str, dcd_file: str):
        print("Loading trajectory...")
        print(f"  PSF: {psf_file}")
        print(f"  DCD: {dcd_file}")
        self.universe = mda.Universe(psf_file, dcd_file)
        print(f"  Atoms: {len(self.universe.atoms)}")
        print(f"  Frames: {len(self.universe.trajectory)}")

        # Cache PSF charges (topology reference, frame-independent)
        self._psf_charges: Dict[int, float] = {
            atom.index: float(atom.charge) for atom in self.universe.atoms
        }

        # Ensure tempfactors attribute exists — PSF files do not carry it
        try:
            _ = self.universe.atoms.tempfactors
        except AttributeError:
            self.universe.add_TopologyAttr(
                'tempfactors', np.zeros(len(self.universe.atoms))
            )
            print("  Note: tempfactors not found in PSF — initialized to 0")

    # ------------------------------------------------------------------
    # Element helpers
    # ------------------------------------------------------------------

    def _get_element_from_mass(self, mass: float) -> str:
        for (lo, hi), elem in self.MASS_TO_ELEMENT.items():
            if lo <= mass <= hi:
                return elem
        if mass < 2.0:             return 'H'
        elif 11.0 <= mass <= 13.0: return 'C'
        elif 13.5 <= mass <= 15.0: return 'N'
        elif 15.5 <= mass <= 17.0: return 'O'
        elif 31.0 <= mass <= 33.0: return 'S'
        else:                      return 'X'

    # ------------------------------------------------------------------
    # Config parsers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_axes(axes_config) -> Tuple[bool, bool, bool]:
        """
        Parse supercell_axes from config into (expand_x, expand_y, expand_z).
        Accepts a list or comma-separated string: x/a, y/b, z/c (case-insensitive).
        Examples: ['x','y'], 'x,y,z', ['a','b','c'], 'z'
        """
        if not axes_config:
            return False, False, False
        if isinstance(axes_config, str):
            tokens = [t.strip().lower() for t in axes_config.split(',')]
        else:
            tokens = [str(t).strip().lower() for t in axes_config]
        return (
            any(t in ('x', 'a') for t in tokens),
            any(t in ('y', 'b') for t in tokens),
            any(t in ('z', 'c') for t in tokens),
        )

    @staticmethod
    def _parse_pdb_stride(value) -> Optional[int]:
        """
        Parse pdb_stride config value.
          'all'   -> 1  (every frame)
          'half'  -> 2  (every other frame)
          'tenth' -> 10 (every 10th frame)
          integer -> that integer
          null / absent -> None (disabled)
        """
        if value is None:
            return None
        if isinstance(value, int):
            return value
        mapping = {'all': 1, 'half': 2, 'tenth': 10}
        v = str(value).strip().lower()
        if v in mapping:
            return mapping[v]
        try:
            return int(v)
        except ValueError:
            raise ValueError(
                f"pdb_stride '{value}' not recognised. "
                f"Use: all, half, tenth, or an integer."
            )

    # ------------------------------------------------------------------
    # Coordinate / link-atom extraction
    # ------------------------------------------------------------------

    def extract_coordinates(self, qm_selection: str, frame: int):
        self.universe.trajectory[frame]
        qm_atoms = self.universe.select_atoms(qm_selection)
        coords = [
            (self._get_element_from_mass(m), p[0], p[1], p[2])
            for m, p in zip(qm_atoms.masses, qm_atoms.positions)
        ]
        for qm_idx, mm_idx in self._find_boundary_bonds(qm_atoms):
            try:
                lp = self._place_link_atom(qm_idx, mm_idx, frame)
                coords.append(('H', lp[0], lp[1], lp[2]))
            except ValueError:
                pass
        return coords

    # ------------------------------------------------------------------
    # Minimum image remapping
    # ------------------------------------------------------------------

    def _remap_position(self, pos: np.ndarray, qm_center: np.ndarray,
                        box: np.ndarray) -> np.ndarray:
        """Remap a single position to its minimum image relative to QM centroid."""
        lx, ly, lz = box[0], box[1], box[2]
        return np.array([
            pos[0] - np.round((pos[0] - qm_center[0]) / lx) * lx,
            pos[1] - np.round((pos[1] - qm_center[1]) / ly) * ly,
            pos[2] - np.round((pos[2] - qm_center[2]) / lz) * lz,
        ])

    def _remap_positions_array(self, pos: np.ndarray, qm_center: np.ndarray,
                               box: np.ndarray) -> np.ndarray:
        """Vectorised remap of an (N, 3) position array to minimum image."""
        lx, ly, lz = box[0], box[1], box[2]
        out = pos.copy()
        out[:, 0] -= np.round((pos[:, 0] - qm_center[0]) / lx) * lx
        out[:, 1] -= np.round((pos[:, 1] - qm_center[1]) / ly) * ly
        out[:, 2] -= np.round((pos[:, 2] - qm_center[2]) / lz) * lz
        return out

    # ------------------------------------------------------------------
    # Supercell image tiling
    # ------------------------------------------------------------------

    def _image_shells(self, cutoff: float, box: np.ndarray,
                      expand: Tuple[bool, bool, bool]) -> Tuple[int, int, int]:
        """ceil(cutoff / L) per requested axis; 0 for suppressed axes."""
        return tuple(
            int(np.ceil(cutoff / box[i])) if do_expand else 0
            for i, do_expand in enumerate(expand)
        )

    def _tile_images(self, charges: list, qm_pos: np.ndarray,
                     cutoff: float, box: np.ndarray,
                     expand: Tuple[bool, bool, bool]
                     ) -> Tuple[list, Tuple[int, int, int], int]:
        """
        Generate periodic images along requested axes.
        Primary charges must already be remapped to minimum image positions
        relative to QM centroid so the (0,0,0) skip matches the primary set.

        Returns (image_charges, shells, n_candidates).
        """
        if not charges:
            return [], (0, 0, 0), 0

        lx, ly, lz    = box[0], box[1], box[2]
        nx, ny, nz    = self._image_shells(cutoff, box, expand)
        image_charges = []
        n_candidates  = 0

        rq       = np.array([[x, y, z] for _, x, y, z in charges])   # (N, 3)
        rcharges = np.array([q for q, *_ in charges])                  # (N,)

        for ix in range(-nx, nx + 1):
            for iy in range(-ny, ny + 1):
                for iz in range(-nz, nz + 1):
                    if ix == 0 and iy == 0 and iz == 0:
                        continue
                    shifted      = rq + np.array([ix*lx, iy*ly, iz*lz])
                    n_candidates += len(shifted)
                    dists        = distances.distance_array(
                        qm_pos, shifted, box=None
                    ).min(axis=0)
                    for idx in np.where(dists <= cutoff)[0]:
                        image_charges.append((
                            float(rcharges[idx]),
                            float(shifted[idx, 0]),
                            float(shifted[idx, 1]),
                            float(shifted[idx, 2]),
                        ))

        return image_charges, (nx, ny, nz), n_candidates

    # ------------------------------------------------------------------
    # Point-charge extraction
    # ------------------------------------------------------------------

    def extract_point_charges(self, qm_selection: str, cutoff: float,
                               frame: int, boundary_scheme: str,
                               switchdist: Optional[float] = None,
                               expand: Tuple[bool, bool, bool] = (False, False, False)):
        """
        Returns
        -------
        charges     : list of (q, x, y, z) — primary + image, switched
        mods        : list of ChargeMod
        switch_recs : list of SwitchRecord
        image_info  : dict with image statistics
        mm_ag       : MDAnalysis AtomGroup of MM atoms within cutoff
        qm_center   : np.ndarray (3,) — QM centroid used for remapping
        box         : np.ndarray (6,) — box dimensions at this frame
        """
        self.universe.trajectory[frame]
        qm_atoms   = self.universe.select_atoms(qm_selection)
        all_atoms  = self.universe.select_atoms("all")
        qm_idx_set = set(qm_atoms.indices)
        mm_atoms   = [a for a in all_atoms if a.index not in qm_idx_set]
        qm_pos     = qm_atoms.positions
        box        = self.universe.dimensions

        if not mm_atoms:
            return [], [], [], {}, self.universe.atoms[[]], qm_pos.mean(axis=0), box

        mm_positions  = np.array([a.position for a in mm_atoms])

        # PBC-aware primary selection — minimum image convention via box=dimensions
        dist_matrix   = distances.distance_array(qm_pos, mm_positions, box=box)
        min_distances = dist_matrix.min(axis=0)

        # Whole-residue inclusion
        residues_in = set()
        for i, atom in enumerate(mm_atoms):
            if min_distances[i] <= cutoff:
                residues_in.add((atom.segid, atom.resid))

        mm_cut = [a for a in mm_atoms if (a.segid, a.resid) in residues_in]
        mm_ag  = self.universe.atoms[np.array([a.index for a in mm_cut])]

        boundary_bonds = self._find_boundary_bonds(qm_atoms)

        if not boundary_bonds or boundary_scheme == 'NONE':
            primary_charges = [(a.charge, *a.position) for a in mm_cut]
            raw_mods        = []
        else:
            primary_charges, raw_mods = self._apply_boundary_scheme(
                mm_cut, boundary_bonds, boundary_scheme
            )

        # Remap primary charges to minimum image positions relative to QM centroid.
        # distance_array(box=dimensions) correctly finds atoms across PBC but keeps
        # raw trajectory coordinates. Remapping brings them physically near QM so:
        #  (a) written coordinates match the geometry seen by the QM program, and
        #  (b) the (0,0,0) skip in _tile_images matches the primary set exactly.
        qm_center = qm_pos.mean(axis=0)
        lx, ly, lz = box[0], box[1], box[2]
        primary_charges = [
            (q,
             x - np.round((x - qm_center[0]) / lx) * lx,
             y - np.round((y - qm_center[1]) / ly) * ly,
             z - np.round((z - qm_center[2]) / lz) * lz)
            for q, x, y, z in primary_charges
        ]

        image_info = {}
        if any(expand):
            image_charges, shells, n_cand = self._tile_images(
                primary_charges, qm_pos, cutoff, box, expand
            )
            nx, ny, nz = shells
            image_info = {
                'nx': nx, 'ny': ny, 'nz': nz,
                'n_images': len(image_charges),
                'n_candidates': n_cand,
                'lx': box[0], 'ly': box[1], 'lz': box[2],
            }
            all_charges = primary_charges + image_charges
        else:
            all_charges = primary_charges

        # Switching applied to all charges; box=None since positions are explicit
        switch_recs = []
        if switchdist is not None:
            all_charges, switch_recs = self._apply_switching_to_charges(
                all_charges, qm_pos, switchdist, cutoff,
                box=None, frame=frame, n_primary=len(primary_charges),
            )

        mods = self._build_charge_mods(raw_mods, frame, qm_center, box)
        return all_charges, mods, switch_recs, image_info, mm_ag, qm_center, box

    # ------------------------------------------------------------------
    # Boundary helpers
    # ------------------------------------------------------------------

    def _find_boundary_bonds(self, qm_atoms):
        qm_idx = set(qm_atoms.indices)
        boundary = []
        for atom in qm_atoms:
            if hasattr(atom, 'bonds'):
                for bond in atom.bonds:
                    other = bond.partner(atom)
                    if other.index not in qm_idx:
                        boundary.append((atom.index, other.index))
        return boundary

    def _place_link_atom(self, qm_idx: int, mm_idx: int, frame: int):
        self.universe.trajectory[frame]
        qm_pos = self.universe.atoms[qm_idx].position
        mm_pos = self.universe.atoms[mm_idx].position
        vec    = mm_pos - qm_pos
        vlen   = np.linalg.norm(vec)
        if vlen < 0.1:
            raise ValueError("Bond too short")
        return qm_pos + (vec / vlen) * 1.09

    def _get_bonded_atoms(self, atom_idx: int):
        atom = self.universe.atoms[atom_idx]
        return [bond.partner(atom).index for bond in atom.bonds] \
               if hasattr(atom, 'bonds') else []

    def _apply_switching_to_charges(self, charges, qm_pos, sw, cut, box,
                                    frame, n_primary: int = None):
        """
        Vectorised quintic switching. box=None when positions are explicit.
        n_primary marks the boundary between primary and image charges.
        Only charges with scale < 1 (i.e. in the switching zone) are recorded.
        """
        if not charges:
            return [], []

        positions = np.array([[x, y, z] for _, x, y, z in charges])
        qs        = np.array([q for q, *_ in charges])

        all_dists = distances.distance_array(
            qm_pos, positions, box=box
        ).min(axis=0)

        # Fully vectorised quintic: 1 - 10t^3 + 15t^4 - 6t^5
        sw_range = cut - sw
        t        = np.clip((all_dists - sw) / sw_range, 0.0, 1.0)
        scales   = np.where(all_dists <= sw, 1.0,
                   np.where(all_dists >= cut, 0.0,
                   1.0 - 10*t**3 + 15*t**4 - 6*t**5))

        scaled_qs = qs * scales
        scaled = [
            (float(scaled_qs[i]), float(positions[i, 0]),
             float(positions[i, 1]), float(positions[i, 2]))
            for i in range(len(charges))
        ]

        recs = []
        for i in np.where(scales < 1.0)[0]:
            is_img = (n_primary is not None) and (int(i) >= n_primary)
            recs.append(SwitchRecord(
                frame         = frame,
                psf_charge    = float(qs[i]),
                scaled_charge = float(scaled_qs[i]),
                scale         = float(scales[i]),
                dist          = float(all_dists[i]),
                position      = positions[i],
                is_image      = is_img,
            ))
        return scaled, recs

    def _apply_boundary_scheme(self, mm_atoms, boundary_bonds, scheme):
        atom_map         = {a.index: a for a in mm_atoms}
        modified_charges = {idx: a.charge for idx, a in atom_map.items()}
        virtual_charges  = []
        removed_atoms    = set()
        charge_mods      = []

        for qm_idx, mm1_idx in boundary_bonds:
            if mm1_idx not in atom_map:
                continue

            mm1_atom   = self.universe.atoms[mm1_idx]
            mm1_charge = mm1_atom.charge
            mm1_pos    = mm1_atom.position
            mm2_atoms  = self._get_bonded_atoms(mm1_idx)
            mm2_cut    = [i for i in mm2_atoms if i in atom_map]

            if scheme == 'Z1':
                removed_atoms.add(mm1_idx)
                charge_mods.append({
                    'type': 'removed', 'atom': mm1_atom,
                    'old_charge': mm1_charge, 'new_charge': 0.0,
                    'reason': 'MM1 removed (Z1)',
                })

            elif scheme == 'Z2':
                # Zero MM1 and all MM2 atoms
                removed_atoms.add(mm1_idx)
                charge_mods.append({
                    'type': 'removed', 'atom': mm1_atom,
                    'old_charge': mm1_charge, 'new_charge': 0.0,
                    'reason': 'MM1 removed (Z2)',
                })
                for mm2_idx in mm2_cut:
                    mm2_atom = self.universe.atoms[mm2_idx]
                    old_q    = mm2_atom.charge
                    removed_atoms.add(mm2_idx)
                    charge_mods.append({
                        'type': 'removed', 'atom': mm2_atom,
                        'old_charge': old_q, 'new_charge': 0.0,
                        'reason': 'MM2 removed (Z2)',
                    })

            elif scheme == 'Z3':
                # Zero MM1, all MM2, and all MM3 atoms
                removed_atoms.add(mm1_idx)
                charge_mods.append({
                    'type': 'removed', 'atom': mm1_atom,
                    'old_charge': mm1_charge, 'new_charge': 0.0,
                    'reason': 'MM1 removed (Z3)',
                })
                seen_mm3 = set()  # guard against duplicate log entries
                for mm2_idx in mm2_cut:
                    mm2_atom = self.universe.atoms[mm2_idx]
                    old_q    = mm2_atom.charge
                    removed_atoms.add(mm2_idx)
                    charge_mods.append({
                        'type': 'removed', 'atom': mm2_atom,
                        'old_charge': old_q, 'new_charge': 0.0,
                        'reason': 'MM2 removed (Z3)',
                    })
                    mm3_atoms = self._get_bonded_atoms(mm2_idx)
                    mm3_cut   = [i for i in mm3_atoms
                                 if i in atom_map
                                 and i != mm1_idx
                                 and i not in mm2_cut   # skip other MM2 atoms
                                 and i not in seen_mm3]  # skip already-logged MM3
                    for mm3_idx in mm3_cut:
                        seen_mm3.add(mm3_idx)
                        mm3_atom = self.universe.atoms[mm3_idx]
                        old_q3   = mm3_atom.charge
                        removed_atoms.add(mm3_idx)
                        charge_mods.append({
                            'type': 'removed', 'atom': mm3_atom,
                            'old_charge': old_q3, 'new_charge': 0.0,
                            'reason': 'MM3 removed (Z3)',
                        })

            elif scheme == 'RCD':
                removed_atoms.add(mm1_idx)
                charge_mods.append({
                    'type': 'removed', 'atom': mm1_atom,
                    'old_charge': mm1_charge, 'new_charge': 0.0,
                    'reason': 'MM1 removed (RCD)',
                })
                n = len(mm2_cut)
                if n:
                    for mm2_idx in mm2_cut:
                        mm2_pos  = self.universe.atoms[mm2_idx].position
                        midpoint = (mm1_pos + mm2_pos) * 0.5
                        # Factor of 2: preserves MM1-MM2 bond dipole
                        # q*d (charge at MM1) == 2q*(d/2) (charge at midpoint)
                        vq = 2.0 * mm1_charge / n
                        virtual_charges.append((vq, *midpoint))
                        charge_mods.append({
                            'type': 'virtual', 'position': midpoint,
                            'charge': vq, 'reason': 'Virtual RCD midpoint',
                        })
                        mm2_atom = self.universe.atoms[mm2_idx]
                        old_q    = mm2_atom.charge
                        modified_charges[mm2_idx] -= mm1_charge / n
                        charge_mods.append({
                            'type': 'modified', 'atom': mm2_atom,
                            'old_charge': old_q,
                            'new_charge': modified_charges[mm2_idx],
                            'reason': 'MM2 adjusted (RCD)',
                        })

            elif scheme == 'CS':
                removed_atoms.add(mm1_idx)
                charge_mods.append({
                    'type': 'removed', 'atom': mm1_atom,
                    'old_charge': mm1_charge, 'new_charge': 0.0,
                    'reason': 'MM1 removed (CS)',
                })
                n = len(mm2_cut)
                if n:
                    split = mm1_charge / n
                    for mm2_idx in mm2_cut:
                        mm2_pos = self.universe.atoms[mm2_idx].position
                        vec     = mm2_pos - mm1_pos
                        vn      = np.linalg.norm(vec)
                        if vn > 0.1:
                            u = vec / vn
                            virtual_charges.append((split,  *(mm2_pos - u * 0.3)))
                            virtual_charges.append((-split, *(mm2_pos + u * 0.3)))
                        mm2_atom = self.universe.atoms[mm2_idx]
                        old_q    = mm2_atom.charge
                        modified_charges[mm2_idx] += split
                        charge_mods.append({
                            'type': 'modified', 'atom': mm2_atom,
                            'old_charge': old_q,
                            'new_charge': modified_charges[mm2_idx],
                            'reason': 'MM2 shifted (CS)',
                        })

        charges = [
            (modified_charges[idx], *atom_map[idx].position)
            for idx in modified_charges if idx not in removed_atoms
        ]
        charges.extend(virtual_charges)
        return charges, charge_mods

    # ------------------------------------------------------------------
    # Build typed ChargeMod objects
    # ------------------------------------------------------------------

    def _build_charge_mods(self, raw_mods: list, frame: int,
                           qm_center: np.ndarray,
                           box: np.ndarray) -> List[ChargeMod]:
        out = []
        for d in raw_mods:
            if d['type'] == 'virtual':
                out.append(ChargeMod(
                    frame          = frame,
                    mod_type       = 'virtual',
                    reason         = d['reason'],
                    psf_charge     = 0.0,
                    applied_charge = d['charge'],
                    position       = self._remap_position(
                        np.array(d['position']), qm_center, box),
                ))
            else:
                atom  = d['atom']
                psf_q = self._psf_charges.get(atom.index, d['old_charge'])
                out.append(ChargeMod(
                    frame          = frame,
                    mod_type       = d['type'],
                    reason         = d['reason'],
                    psf_charge     = psf_q,
                    applied_charge = d['new_charge'],
                    position       = self._remap_position(
                        atom.position.copy(), qm_center, box),
                    atom_index     = atom.index,
                    segid          = atom.segid,
                    resid          = atom.resid,
                    resname        = atom.resname,
                    name           = atom.name,
                ))
        return out

    # ------------------------------------------------------------------
    # Structure writer (PSF + PDB)
    # ------------------------------------------------------------------

    def _write_structure(self, frame: int, qm_atoms, mm_ag, base: Path,
                         psf_source: str, qm_center: np.ndarray,
                         box: np.ndarray):
        """
        Write a full-system PDB and a paired PSF for every frame.

        MM atoms within cutoff are temporarily remapped to their minimum image
        positions relative to the QM centroid before writing — matching the
        remapping applied to the point charges in the .in file — so the PDB
        coordinates are consistent with what was written to the QM/MM input.

        Beta (tempfactor) column encodes region:
          9 = QM atom
          8 = MM point charge (within cutoff)
          0 = outside cutoff

        In VMD: load *_struct.psf first, then *_struct.pdb as coordinates.
        Color by Beta to visualise QM/MM partitioning.
        """
        qm_idx = set(qm_atoms.indices)
        mm_idx = set(mm_ag.indices)

        all_atoms = self.universe.select_atoms("all")

        # Temporarily remap mm_ag positions to minimum image relative to QM
        orig_mm_pos    = mm_ag.positions.copy()
        mm_ag.positions = self._remap_positions_array(orig_mm_pos, qm_center, box)

        # Vectorised beta assignment; QM after MM so it overwrites any overlap
        beta = np.zeros(len(all_atoms))
        if mm_idx:
            beta[np.array(list(mm_idx), dtype=int)] = 8.0
        if qm_idx:
            beta[np.array(list(qm_idx), dtype=int)] = 9.0

        orig_tf = all_atoms.tempfactors.copy()
        all_atoms.tempfactors = beta

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            all_atoms.write(f"{base}_struct.pdb")

        # Restore universe state so subsequent frame reads are unaffected
        all_atoms.tempfactors = orig_tf
        mm_ag.positions        = orig_mm_pos

        shutil.copy(psf_source, f"{base}_struct.psf")

    # ------------------------------------------------------------------
    # Log writers
    # ------------------------------------------------------------------

    def _write_boundary_log(self, fh, all_mods: List[ChargeMod]):
        fh.write("=" * 72 + "\n")
        fh.write("ezQMMM Boundary Charge Modification Detail\n")
        fh.write("=" * 72 + "\n\n")
        fh.write(f"{'Frame':<8} {'Type':<10} {'Atom / Site':<28} "
                 f"{'PSF q':>8} {'Applied q':>10} {'Delta q':>10}  Reason\n")
        fh.write("-" * 72 + "\n")

        for m in sorted(all_mods, key=lambda x: (x.frame, x.mod_type)):
            tag = (f"{m.segid}:{m.resid} {m.resname} {m.name}"
                   if m.atom_index is not None
                   else f"virtual@{m.position[0]:.1f},{m.position[1]:.1f}")
            fh.write(
                f"{m.frame:<8d} {m.mod_type:<10s} {tag:<28s} "
                f"{m.psf_charge:+8.4f} {m.applied_charge:+10.4f} "
                f"{m.delta:+10.4f}  {m.reason}\n"
            )

        n_rem = sum(1 for m in all_mods if m.mod_type == 'removed')
        n_mod = sum(1 for m in all_mods if m.mod_type == 'modified')
        n_vir = sum(1 for m in all_mods if m.mod_type == 'virtual')
        fh.write("\n" + "-" * 72 + "\n")
        fh.write(f"Total: removed={n_rem}  modified={n_mod}  "
                 f"virtual={n_vir}  total={len(all_mods)}\n")

    def _write_switching_log(self, fh, all_switch: List[SwitchRecord],
                             switchdist: float, cutoff: float,
                             expand: Tuple[bool, bool, bool]):
        supercell_on = any(expand)
        axis_labels  = ','.join(l for l, e in zip(('x', 'y', 'z'), expand) if e)
        fh.write("=" * 72 + "\n")
        fh.write("ezQMMM Switching-Function Charge Modifications\n")
        fh.write(f"  Switching zone  : {switchdist:.2f} - {cutoff:.2f} Ang\n")
        fh.write(f"  Function        : quintic (NAMD-style)\n")
        fh.write(f"  Distance metric : minimum distance to any QM atom "
                 f"(not center of mass)\n")
        if supercell_on:
            fh.write(f"  Image charges   : included (axes: {axis_labels}); "
                     f"marked in Src column as IMG\n")
        fh.write(f"  Recorded        : only charges with scale < 1 "
                 f"(dist > {switchdist:.2f} Ang)\n")
        fh.write("=" * 72 + "\n\n")

        src_col = "  Src" if supercell_on else ""
        fh.write(f"{'Frame':<8} {'Dist(Ang)':>9} {'Scale':>8} "
                 f"{'q_orig':>10} {'q_scaled':>10} {'Delta q':>10}  "
                 f"x(Ang)    y(Ang)    z(Ang){src_col}\n")
        fh.write("-" * 72 + "\n")

        prev_frame = None
        for r in sorted(all_switch, key=lambda x: (x.frame, x.dist)):
            if r.frame != prev_frame:
                fh.write(f"\n--- Frame {r.frame} ---\n")
                prev_frame = r.frame
            delta   = r.scaled_charge - r.psf_charge
            src_tag = "  IMG" if (supercell_on and r.is_image) else "     "
            fh.write(
                f"{r.frame:<8d} {r.dist:>9.3f} {r.scale:>8.5f} "
                f"{r.psf_charge:>+10.4f} {r.scaled_charge:>+10.4f} "
                f"{delta:>+10.4f}  "
                f"{r.position[0]:>8.3f}  {r.position[1]:>8.3f}  "
                f"{r.position[2]:>8.3f}{src_tag}\n"
            )

        n_prim = sum(1 for r in all_switch if not r.is_image)
        n_img  = sum(1 for r in all_switch if r.is_image)
        fh.write("\n" + "-" * 72 + "\n")
        if supercell_on:
            fh.write(f"Total switching-zone events: {len(all_switch)}  "
                     f"(primary={n_prim}, image={n_img})\n")
        else:
            fh.write(f"Total switching-zone charge events: {len(all_switch)}\n")

    # ------------------------------------------------------------------
    # Main generate loop
    # ------------------------------------------------------------------

    def generate(self, config: Dict):
        qm_sel        = config['qm_selection']
        mm_cutoff     = config.get('mm_cutoff', 40.0)
        default_sw    = mm_cutoff - 5.0 if mm_cutoff > 5.0 else mm_cutoff * 0.75
        mm_switchdist = config.get('mm_switchdist', default_sw)
        expand        = self._parse_axes(config.get('supercell_axes', []))
        supercell_on  = any(expand)
        pdb_stride    = self._parse_pdb_stride(config.get('pdb_stride', None))

        first  = config.get('first_frame', 0)
        last   = config.get('last_frame', -1)
        if last == -1 or last >= len(self.universe.trajectory):
            last = len(self.universe.trajectory) - 1
        stride = config.get('stride', 1)

        method  = config.get('method', 'B3LYP')
        basis   = config.get('basis', '6-31G*')
        charge  = config.get('charge', 0)
        mult    = config.get('multiplicity', 1)
        bscheme = config.get('boundary_scheme', 'RCD').upper()

        output_dir = Path(config.get('output_dir', '.'))
        prefix     = config.get('output_prefix', 'qmmm')

        if 'program' not in config:
            raise ValueError("'program' required (orca/qchem/psi4)")

        valid_schemes = {'RCD', 'CS', 'Z1', 'Z2', 'Z3', 'NONE'}
        if bscheme not in valid_schemes:
            raise ValueError(
                f"boundary_scheme '{bscheme}' not recognised. "
                f"Valid options: {', '.join(sorted(valid_schemes))}"
            )

        program       = config['program'].lower()
        keywords      = config.get(f'{program}_keywords', '')
        custom_blocks = config.get(f'{program}_blocks', '')

        output_dir.mkdir(parents=True, exist_ok=True)

        frames = list(range(first, last + 1, stride))
        print(f"\nSettings:")
        print(f"  Program     : {program.upper()}")
        print(f"  QM          : {method}/{basis}  charge={charge}  mult={mult}")
        print(f"  Boundary    : {bscheme}")
        print(f"  MM cutoff   : {mm_cutoff} Ang,  switch: {mm_switchdist} Ang")
        if supercell_on:
            self.universe.trajectory[first]
            box = self.universe.dimensions
            nx, ny, nz  = self._image_shells(mm_cutoff, box, expand)
            axis_labels = [l for l, e in zip(('x', 'y', 'z'), expand) if e]
            print(f"  Supercell   : axes={','.join(axis_labels)}  |  "
                  f"shells x={nx} (Lx={box[0]:.2f} Ang), "
                  f"y={ny} (Ly={box[1]:.2f} Ang), "
                  f"z={nz} (Lz={box[2]:.2f} Ang)  [first frame]")
        if pdb_stride:
            print(f"  PDB/PSF     : every {pdb_stride} frame(s)")
        print(f"  Frames      : {len(frames)}")

        all_mods:   List[ChargeMod]    = []
        all_switch: List[SwitchRecord] = []
        generated = []

        for i, frame in enumerate(frames, 1):
            coords = self.extract_coordinates(qm_sel, frame)

            charges, mods, sw_recs, img_inf, mm_ag, qm_center, box = \
                self.extract_point_charges(
                    qm_sel, mm_cutoff, frame, bscheme, mm_switchdist, expand
                )

            all_mods.extend(mods)
            all_switch.extend(sw_recs)

            img_str = (f"  img={img_inf.get('n_images', 0):5d}"
                       if supercell_on else "")
            print(f"  [{i}/{len(frames)}] frame {frame:5d}: "
                  f"QM={len(coords):4d}  MM={len(charges):6d}{img_str}  "
                  f"mods={len(mods):3d}  switched={len(sw_recs):4d}")

            base = output_dir / f"{prefix}_frame{frame}"

            if program == 'orca':
                fname = f"{base}_orca.inp"
                self._write_orca(fname, coords, charges, method, basis,
                                 charge, mult, keywords, custom_blocks)
            elif program == 'qchem':
                fname = f"{base}_qchem.in"
                self._write_qchem(fname, coords, charges, method, basis,
                                  charge, mult, keywords, custom_blocks)
            elif program == 'psi4':
                fname = f"{base}_psi4.dat"
                self._write_psi4(fname, coords, charges, method, basis,
                                 charge, mult, keywords, custom_blocks)
            else:
                raise ValueError(f"Unknown program: {program}")
            generated.append(fname)

            if pdb_stride and (i % pdb_stride == 0 or i == 1):
                qm_atoms = self.universe.select_atoms(qm_sel)
                self._write_structure(frame, qm_atoms, mm_ag, base,
                                      config['psf_file'], qm_center, box)

        if all_mods:
            bpath = output_dir / f"{prefix}_boundary.log"
            with open(bpath, 'w') as fh:
                self._write_boundary_log(fh, all_mods)
            print(f"\n  Boundary log  -> {bpath}")

        spath = output_dir / f"{prefix}_switching.log"
        with open(spath, 'w') as fh:
            self._write_switching_log(fh, all_switch, mm_switchdist,
                                      mm_cutoff, expand)
        print(f"  Switching log -> {spath}")

        print(f"\nGenerated {len(generated)} input files")
        return generated

    # ------------------------------------------------------------------
    # QM/MM input writers
    # ------------------------------------------------------------------

    def _write_orca(self, fname, coords, charges, method, basis,
                    charge, mult, keywords, custom_blocks):
        with open(fname, 'w') as f:
            f.write(f"! {method} {basis}\n")
            for line in (keywords or '').strip().split('\n'):
                if line.strip():
                    f.write(f"! {line.strip()}\n")
            f.write("\n")
            if custom_blocks and custom_blocks.strip():
                f.write(custom_blocks.strip() + "\n\n")
            if charges:
                pc = fname.replace('.inp', '_charges.pc')
                f.write(f'%pointcharges "{Path(pc).name}"\n\n')
                with open(pc, 'w') as pf:
                    pf.write(f"{len(charges)}\n")
                    for q, x, y, z in charges:
                        pf.write(f"{q:.4f}  {x:.3f}  {y:.3f}  {z:.3f}\n")
            f.write(f"* xyz {charge} {mult}\n")
            for elem, x, y, z in coords:
                f.write(f"{elem:<4s}  {x:.3f}  {y:.3f}  {z:.3f}\n")
            f.write("*\n")

    def _write_qchem(self, fname, coords, charges, method, basis,
                     charge, mult, keywords, custom_blocks):
        with open(fname, 'w') as f:
            f.write("$molecule\n")
            f.write(f"{charge} {mult}\n")
            for elem, x, y, z in coords:
                f.write(f"{elem:<4s}  {x:.3f}  {y:.3f}  {z:.3f}\n")
            f.write("$end\n\n")
            f.write("$rem\n")
            f.write("   jobtype              sp\n")
            f.write(f"   method               {method}\n")
            f.write(f"   basis                {basis}\n")
            if charges:
                f.write("   qm_mm                true\n")
            for line in (keywords or '').strip().split('\n'):
                if line.strip():
                    f.write(f"   {line.strip()}\n")
            f.write("$end\n")
            if custom_blocks and custom_blocks.strip():
                f.write("\n" + custom_blocks.strip() + "\n")
            if charges:
                f.write("\n$external_charges\n")
                for q, x, y, z in charges:
                    f.write(f"{x:.3f}  {y:.3f}  {z:.3f}  {q:.4f}\n")
                f.write("$end\n")

    def _write_psi4(self, fname, coords, charges, method, basis,
                    charge, mult, keywords, custom_blocks):
        with open(fname, 'w') as f:
            f.write("memory 4 GB\n\n")
            f.write("molecule qmmm {\n")
            f.write(f"  {charge} {mult}\n")
            for elem, x, y, z in coords:
                f.write(f"  {elem:<4s}  {x:.3f}  {y:.3f}  {z:.3f}\n")
            if charges:
                f.write("  no_com\n  no_reorient\n")
            f.write("}\n\n")
            f.write("set {\n")
            f.write(f"  basis {basis}\n")
            for line in (keywords or '').strip().split('\n'):
                if line.strip():
                    f.write(f"  {line.strip()}\n")
            f.write("}\n\n")
            if custom_blocks and custom_blocks.strip():
                f.write(custom_blocks.strip() + "\n\n")
            if charges:
                f.write("Chrgfield = QMMM()\n")
                for q, x, y, z in charges:
                    xb, yb, zb = x * 1.88973, y * 1.88973, z * 1.88973
                    f.write(f"Chrgfield.extern.addCharge"
                            f"({q:.4f}, {xb:.3f}, {yb:.3f}, {zb:.3f})\n")
                f.write("psi4.set_global_option_python('EXTERN', Chrgfield.extern)\n\n")
            f.write(f"energy('{method}')\n")


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

def create_example_config():
    config = """# ezQMMM 2.0 Configuration
psf_file: system.psf
dcd_file: trajectory.dcd
qm_selection: "resid 100 and not backbone"
mm_cutoff: 40.0
mm_switchdist: 35.0
first_frame: 0
last_frame: 100
stride: 10
method: B3LYP
basis: 6-31G*
charge: 0
multiplicity: 1
boundary_scheme: RCD
output_dir: ./qmmm_calculations
output_prefix: qmmm
program: qchem

# Axes along which to tile periodic images (any combo of x/y/z or a/b/c).
# Leave empty or omit to disable. Examples:
#   supercell_axes: [x, y]       # membrane, normal along z
#   supercell_axes: [x, y, z]    # full 3D periodic
#   supercell_axes: z            # expand along z only
supercell_axes: []

# Write PSF/PDB structure files (QM + MM atoms within cutoff).
# Options: all, half, tenth, or any integer stride over generated frames.
# Leave absent or null to disable.
pdb_stride: all

qchem_keywords: |
  scf_convergence    8
  mem_total          8000

qchem_blocks: |
  $basis
  C 0
  cc-pVDZ
  ****
  H 0
  cc-pVDZ
  ****
  $end
"""
    with open('config_example.yaml', 'w') as f:
        f.write(config)
    print("Created: config_example.yaml")


def main():
    if len(sys.argv) < 2:
        print("ezQMMM 2.0 - Easy QM/MM Input Generator")
        print("\nUsage:")
        print("  python ezQMMM2.py config.yaml")
        print("  python ezQMMM2.py --example")
        sys.exit(1)

    if sys.argv[1] == '--example':
        create_example_config()
        return

    try:
        with open(sys.argv[1]) as f:
            config = yaml.safe_load(f)
        gen = QMMMGenerator(config['psf_file'], config['dcd_file'])
        gen.generate(config)
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
