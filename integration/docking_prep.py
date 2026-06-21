"""
Docking preparation for OpenTargetAI.

Assembles a ready-to-dock package for AutoDock Vina from a predicted target and
the query ligand:

    <output>/
        protein/        receptor coordinate file (PDB, if available)
        ligand/         3D-embedded, energy-minimised query ligand (SDF + PDB)
        config/         vina_config.txt with a suggested grid box
        README.txt      step-by-step instructions to finish preparation & dock

The suggested search box is centred on a co-crystallised ligand if one is found
in the receptor PDB (the most reliable binding-site proxy), otherwise on the
receptor centroid with a generous default size. Conversion of PDB → PDBQT is
left to the user's AutoDockTools / Meeko install (documented in the README),
since those binaries are not bundled.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)

# Solvent / ion / cryoprotectant / glycan / common-cofactor residues to ignore
# when searching for the genuine binding-site ligand.
_IGNORE_HET = {
    # water & ions
    "HOH", "WAT", "DOD", "NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE",
    "CD", "NI", "CU", "CO", "HG", "BR", "IOD", "FE2", "SR", "CS", "BA",
    # buffers / cryoprotectants / additives
    "SO4", "PO4", "GOL", "EDO", "ACT", "DMS", "PEG", "PG4", "PGE", "1PE",
    "MPD", "TRS", "FMT", "MES", "EPE", "BME", "IMD", "CIT", "TLA", "ACY",
    "MRD", "BU3", "P6G", "12P", "15P", "NH4", "FLC", "UNX", "UNL",
    # glycans (glycosylation, not the inhibitor)
    "NAG", "MAN", "BMA", "FUC", "GAL", "GLC", "BGC", "NDG", "SIA", "XYP",
    # ubiquitous cofactors (rarely the docking target of interest)
    "HEM", "NAD", "NAP", "FAD", "FMN", "SAH", "SAM",
}


@dataclass
class GridBox:
    """A suggested AutoDock Vina search box (Ångströms)."""

    center: tuple[float, float, float]
    size: tuple[float, float, float]
    source: str  # how the box was derived (e.g. "bound ligand HETATM ...")


@dataclass
class DockingPackage:
    """Summary of a generated docking package."""

    output_dir: str
    protein_file: Optional[str]
    ligand_sdf: Optional[str]
    ligand_pdb: Optional[str]
    config_file: str
    grid_box: GridBox
    zip_file: Optional[str] = None


class DockingPreparation:
    """Builds AutoDock Vina docking packages."""

    DEFAULT_BOX_SIZE = 24.0

    def prepare(
        self,
        output_dir: str,
        ligand_smiles: str,
        target_name: str = "",
        uniprot: str = "",
        receptor_pdb_path: Optional[str] = None,
        exhaustiveness: int = 8,
        num_modes: int = 9,
    ) -> DockingPackage:
        """Generate a complete docking package and return its description."""
        out = Path(output_dir)
        protein_dir = out / "protein"
        ligand_dir = out / "ligand"
        config_dir = out / "config"
        for d in (protein_dir, ligand_dir, config_dir):
            d.mkdir(parents=True, exist_ok=True)

        # 1. Receptor.
        protein_file: Optional[str] = None
        if receptor_pdb_path and Path(receptor_pdb_path).exists():
            protein_file = str(protein_dir / Path(receptor_pdb_path).name)
            shutil.copy(receptor_pdb_path, protein_file)

        # 2. Ligand: 3D embed + minimise, write SDF and PDB.
        ligand_sdf, ligand_pdb = self._prepare_ligand(ligand_smiles, ligand_dir)

        # 3. Grid box suggestion.
        grid_box = self._suggest_grid_box(protein_file)

        # 4. Vina config + README.
        config_file = self._write_vina_config(
            config_dir, grid_box, exhaustiveness, num_modes,
            has_receptor=protein_file is not None,
        )
        self._write_readme(
            out, target_name, uniprot, grid_box, protein_file is not None
        )

        logger.info("Docking package prepared at %s", out)
        return DockingPackage(
            output_dir=str(out),
            protein_file=protein_file,
            ligand_sdf=ligand_sdf,
            ligand_pdb=ligand_pdb,
            config_file=config_file,
            grid_box=grid_box,
        )

    def export_zip(self, package: DockingPackage) -> str:
        """Zip the prepared package directory; return the archive path."""
        out = Path(package.output_dir)
        archive = shutil.make_archive(str(out), "zip", root_dir=str(out))
        package.zip_file = archive
        logger.info("Docking package zipped to %s", archive)
        return archive

    # ── Ligand preparation ────────────────────────────────────────────────

    def _prepare_ligand(
        self, smiles: str, ligand_dir: Path
    ) -> tuple[Optional[str], Optional[str]]:
        """3D-embed and energy-minimise the ligand; write SDF and PDB."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("Invalid ligand SMILES for docking prep: %s", smiles)
            return None, None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) != 0:
            # Fall back to a random-coordinate embedding.
            AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:  # noqa: BLE001
            try:
                AllChem.UFFOptimizeMolecule(mol)
            except Exception:  # noqa: BLE001
                pass

        sdf_path = ligand_dir / "ligand.sdf"
        pdb_path = ligand_dir / "ligand.pdb"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()
        Chem.MolToPDBFile(mol, str(pdb_path))
        return str(sdf_path), str(pdb_path)

    # ── Grid box ──────────────────────────────────────────────────────────

    def _suggest_grid_box(self, protein_file: Optional[str]) -> GridBox:
        """Derive a search box from the receptor structure if possible."""
        if not protein_file:
            return GridBox(
                center=(0.0, 0.0, 0.0),
                size=(self.DEFAULT_BOX_SIZE,) * 3,
                source="No receptor supplied — using a placeholder box. "
                "Centre it on the binding site after downloading a structure.",
            )

        try:
            het_coords, ca_coords, het_name = self._parse_pdb_coords(protein_file)
        except Exception as exc:  # noqa: BLE001
            logger.debug("PDB parse error: %s", exc)
            het_coords, ca_coords, het_name = [], [], ""

        if het_coords:
            center = _centroid(het_coords)
            extent = _extent(het_coords)
            size = tuple(min(30.0, max(16.0, e + 8.0)) for e in extent)
            return GridBox(
                center=center,
                size=size,  # type: ignore[arg-type]
                source=f"Centred on co-crystallised ligand '{het_name}' "
                f"({len(het_coords)} atoms).",
            )
        if ca_coords:
            return GridBox(
                center=_centroid(ca_coords),
                size=(self.DEFAULT_BOX_SIZE,) * 3,
                source="Centred on the receptor Cα centroid (no bound ligand "
                "detected — verify the box covers the binding site).",
            )
        return GridBox(
            center=(0.0, 0.0, 0.0),
            size=(self.DEFAULT_BOX_SIZE,) * 3,
            source="Could not parse receptor coordinates — set the box manually.",
        )

    @staticmethod
    def _parse_pdb_coords(
        pdb_path: str,
    ) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]], str]:
        """
        Extract candidate bound-ligand and Cα coordinates from a PDB file.

        HETATM records are grouped by *individual residue instance*
        ``(resname, chain, resSeq)`` so that, e.g., dozens of separate glycan
        sugars are not merged into one giant pseudo-ligand. The largest single
        drug-like instance (preferring those in the 12–90 atom range) is
        returned. Falls back to all Cα atoms if no ligand is found.

        Returns ``(ligand_coords, ca_coords, ligand_residue_name)``.
        """
        het: dict[tuple[str, str, str], list[tuple[float, float, float]]] = {}
        ca: list[tuple[float, float, float]] = []
        with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                record = line[:6].strip()
                if record == "HETATM":
                    res = line[17:20].strip()
                    if res in _IGNORE_HET:
                        continue
                    key = (res, line[21:22], line[22:26].strip())
                    coord = _read_xyz(line)
                    if coord:
                        het.setdefault(key, []).append(coord)
                elif record == "ATOM" and line[12:16].strip() == "CA":
                    coord = _read_xyz(line)
                    if coord:
                        ca.append(coord)
        if het:
            # Prefer a drug-like sized instance; otherwise the largest one.
            def score(item: tuple) -> tuple:
                n = len(item[1])
                drug_like = 12 <= n <= 90
                return (drug_like, n)

            best_key, best_coords = max(het.items(), key=lambda kv: score(kv))
            return best_coords, ca, best_key[0]
        return [], ca, ""

    # ── Config & docs ─────────────────────────────────────────────────────

    def _write_vina_config(
        self,
        config_dir: Path,
        box: GridBox,
        exhaustiveness: int,
        num_modes: int,
        has_receptor: bool,
    ) -> str:
        """Write the AutoDock Vina configuration file."""
        cx, cy, cz = box.center
        sx, sy, sz = box.size
        receptor = "../protein/receptor.pdbqt" if has_receptor else "receptor.pdbqt"
        content = f"""# AutoDock Vina configuration generated by OpenTargetAI
# Grid box derivation: {box.source}

receptor = {receptor}
ligand   = ../ligand/ligand.pdbqt

center_x = {cx:.3f}
center_y = {cy:.3f}
center_z = {cz:.3f}

size_x = {sx:.1f}
size_y = {sy:.1f}
size_z = {sz:.1f}

exhaustiveness = {exhaustiveness}
num_modes = {num_modes}
energy_range = 4
"""
        path = config_dir / "vina_config.txt"
        path.write_text(content, encoding="utf-8")
        return str(path)

    def _write_readme(
        self,
        out: Path,
        target_name: str,
        uniprot: str,
        box: GridBox,
        has_receptor: bool,
    ) -> None:
        """Write step-by-step docking instructions."""
        cx, cy, cz = box.center
        receptor_step = (
            "A receptor PDB has been copied to protein/. "
            if has_receptor
            else "No receptor structure was supplied. Download one from the "
            "Structures tab (or RCSB PDB) and place it in protein/. "
        )
        text = f"""OpenTargetAI — AutoDock Vina docking package
=============================================

Target:   {target_name or 'N/A'}
UniProt:  {uniprot or 'N/A'}

Contents
--------
  protein/   receptor structure (PDB)
  ligand/    query ligand, 3D-embedded & minimised (ligand.sdf, ligand.pdb)
  config/    vina_config.txt (suggested grid box)

Suggested search box
--------------------
  center = ({cx:.3f}, {cy:.3f}, {cz:.3f})
  size   = ({box.size[0]:.1f}, {box.size[1]:.1f}, {box.size[2]:.1f}) Angstrom
  basis  = {box.source}

How to finish preparation and dock
----------------------------------
1. {receptor_step}
2. Prepare the receptor as PDBQT (remove waters, add H, assign charges), e.g.
   with AutoDockTools:  prepare_receptor -r protein/receptor.pdb -o protein/receptor.pdbqt
   or with Meeko/Open Babel:  obabel protein/receptor.pdb -O protein/receptor.pdbqt -xr
3. Prepare the ligand as PDBQT, e.g. with Meeko:
   mk_prepare_ligand.py -i ligand/ligand.sdf -o ligand/ligand.pdbqt
   or:  obabel ligand/ligand.sdf -O ligand/ligand.pdbqt
4. Verify / adjust the grid box in config/vina_config.txt so it encloses the
   binding site.
5. Run docking:
   cd config && vina --config vina_config.txt --out ../ligand/docked.pdbqt

NOTE: OpenTargetAI does not bundle AutoDock Vina or AutoDockTools. Install them
separately (https://vina.scripps.edu / https://ccsb.scripps.edu/mgltools).
"""
        (out / "README.txt").write_text(text, encoding="utf-8")


# ── Geometry helpers ───────────────────────────────────────────────────────


def _read_xyz(line: str) -> Optional[tuple[float, float, float]]:
    """Parse the x/y/z columns of a PDB ATOM/HETATM record."""
    try:
        return (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )
    except (ValueError, IndexError):
        return None


def _centroid(
    coords: list[tuple[float, float, float]]
) -> tuple[float, float, float]:
    """Geometric centroid of a list of points."""
    n = len(coords)
    return (
        sum(c[0] for c in coords) / n,
        sum(c[1] for c in coords) / n,
        sum(c[2] for c in coords) / n,
    )


def _extent(
    coords: list[tuple[float, float, float]]
) -> tuple[float, float, float]:
    """Axis-aligned bounding-box extent of a list of points."""
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    return (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
