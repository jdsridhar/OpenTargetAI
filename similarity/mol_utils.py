"""
Molecule utilities for OpenTargetAI.

A single place for parsing input molecules (SMILES / MOL / SDF), computing
physicochemical descriptors and drug-likeness rules, deriving scaffolds and
maximum-common-substructures, and rendering dark-themed 2D depictions as SVG.

All functions are defensive: invalid input returns ``None`` / empty results
rather than raising, so the UI layer never crashes on bad user input.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Draw,
    QED,
    rdMolDescriptors,
)
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold

logger = logging.getLogger(__name__)


# ── Parsing ───────────────────────────────────────────────────────────────


def parse_smiles(smiles: str) -> Optional[Chem.Mol]:
    """Parse and sanitise a SMILES string. Returns ``None`` if invalid."""
    if not smiles or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    return mol


def mol_from_molblock(block: str) -> Optional[Chem.Mol]:
    """Parse a MOL/SDF block of text into a molecule."""
    if not block or not block.strip():
        return None
    try:
        return Chem.MolFromMolBlock(block)
    except Exception as exc:  # noqa: BLE001
        logger.debug("MolFromMolBlock error: %s", exc)
        return None


def mol_from_file(path: str) -> Optional[Chem.Mol]:
    """
    Load a single molecule from a ``.mol``, ``.sdf`` or ``.smi`` file.
    For multi-record SDF/SMILES files the first valid record is returned.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    suffix = p.suffix.lower()
    try:
        if suffix in (".sdf", ".sd"):
            supplier = Chem.SDMolSupplier(str(p))
            for mol in supplier:
                if mol is not None:
                    return mol
            return None
        if suffix == ".mol":
            return Chem.MolFromMolFile(str(p))
        if suffix in (".smi", ".smiles", ".txt"):
            text = p.read_text(encoding="utf-8", errors="replace").strip()
            first = text.splitlines()[0].split()[0] if text else ""
            return parse_smiles(first)
    except Exception as exc:  # noqa: BLE001
        logger.warning("mol_from_file error for %s: %s", path, exc)
    return None


def canonical_smiles(mol: Chem.Mol) -> str:
    """Return the canonical SMILES for a molecule."""
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol)


def inchikey(mol: Chem.Mol) -> Optional[str]:
    """Return the InChIKey for a molecule (or ``None`` on failure)."""
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:  # noqa: BLE001
        return None


# ── Descriptors & drug-likeness ────────────────────────────────────────────


@dataclass
class MoleculeProfile:
    """Computed descriptor and drug-likeness profile of a molecule."""

    smiles: str
    formula: str
    mol_weight: float
    exact_mass: float
    logp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    rings: int
    aromatic_rings: int
    heavy_atoms: int
    fraction_csp3: float
    qed: float
    # Drug-likeness rule outcomes: name -> (passed, detail)
    rules: dict[str, tuple[bool, str]] = field(default_factory=dict)
    alerts: list[str] = field(default_factory=list)

    def as_rows(self) -> list[tuple[str, str]]:
        """Return (label, value) rows for display in a property table."""
        return [
            ("Molecular Formula", self.formula),
            ("Molecular Weight", f"{self.mol_weight:.2f} g/mol"),
            ("Exact Mass", f"{self.exact_mass:.4f}"),
            ("LogP (cLogP)", f"{self.logp:.2f}"),
            ("TPSA", f"{self.tpsa:.2f} Å²"),
            ("H-Bond Donors", str(self.hbd)),
            ("H-Bond Acceptors", str(self.hba)),
            ("Rotatable Bonds", str(self.rotatable_bonds)),
            ("Rings", str(self.rings)),
            ("Aromatic Rings", str(self.aromatic_rings)),
            ("Heavy Atoms", str(self.heavy_atoms)),
            ("Fraction Csp3", f"{self.fraction_csp3:.2f}"),
            ("QED (drug-likeness)", f"{self.qed:.3f}"),
        ]


def compute_profile(mol: Chem.Mol) -> Optional[MoleculeProfile]:
    """Compute the full descriptor + drug-likeness profile for a molecule."""
    if mol is None:
        return None
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        try:
            qed_val = QED.qed(mol)
        except Exception:  # noqa: BLE001
            qed_val = 0.0

        profile = MoleculeProfile(
            smiles=canonical_smiles(mol),
            formula=rdMolDescriptors.CalcMolFormula(mol),
            mol_weight=mw,
            exact_mass=Descriptors.ExactMolWt(mol),
            logp=logp,
            tpsa=tpsa,
            hbd=hbd,
            hba=hba,
            rotatable_bonds=rotb,
            rings=rdMolDescriptors.CalcNumRings(mol),
            aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(mol),
            heavy_atoms=mol.GetNumHeavyAtoms(),
            fraction_csp3=rdMolDescriptors.CalcFractionCSP3(mol),
            qed=qed_val,
        )
        profile.rules = _druglikeness_rules(mw, logp, tpsa, hbd, hba, rotb)
        profile.alerts = structural_alerts(mol)
        return profile
    except Exception as exc:  # noqa: BLE001
        logger.warning("compute_profile error: %s", exc)
        return None


def _druglikeness_rules(
    mw: float, logp: float, tpsa: float, hbd: int, hba: int, rotb: int
) -> dict[str, tuple[bool, str]]:
    """Evaluate common drug-likeness rule sets."""
    # Lipinski's Rule of Five (≤1 violation is generally accepted).
    lip_viol = sum(
        [mw > 500, logp > 5, hbd > 5, hba > 10]
    )
    lipinski = (
        lip_viol <= 1,
        f"{lip_viol} violation(s): MW≤500, LogP≤5, HBD≤5, HBA≤10",
    )

    # Veber oral bioavailability rules.
    veber_ok = rotb <= 10 and tpsa <= 140
    veber = (veber_ok, "RotB ≤ 10 and TPSA ≤ 140 Å²")

    # Ghose filter.
    ghose_ok = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6)
    ghose = (ghose_ok, "160 ≤ MW ≤ 480 and -0.4 ≤ LogP ≤ 5.6")

    # Lead-likeness (Teague).
    lead_ok = mw <= 350 and logp <= 3.5
    lead = (lead_ok, "MW ≤ 350 and LogP ≤ 3.5")

    return {
        "Lipinski (Ro5)": lipinski,
        "Veber": veber,
        "Ghose": ghose,
        "Lead-like": lead,
    }


_PAINS_CATALOG = None


def _get_pains_catalog():
    """Lazily build (and cache) the PAINS filter catalogue."""
    global _PAINS_CATALOG
    if _PAINS_CATALOG is None:
        try:
            from rdkit.Chem import FilterCatalog
            from rdkit.Chem.FilterCatalog import FilterCatalogParams

            params = FilterCatalogParams()
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
            _PAINS_CATALOG = FilterCatalog.FilterCatalog(params)
        except Exception as exc:  # noqa: BLE001
            logger.debug("PAINS catalogue unavailable: %s", exc)
            _PAINS_CATALOG = False
    return _PAINS_CATALOG


def structural_alerts(mol: Chem.Mol) -> list[str]:
    """Return the names of any PAINS structural alerts matched by the molecule."""
    catalog = _get_pains_catalog()
    if not catalog:
        return []
    try:
        matches = catalog.GetMatches(mol)
        return [m.GetDescription() for m in matches]
    except Exception:  # noqa: BLE001
        return []


# ── Scaffolds & substructure ───────────────────────────────────────────────


def murcko_scaffold(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Return the Bemis-Murcko scaffold of a molecule."""
    if mol is None:
        return None
    try:
        return MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:  # noqa: BLE001
        return None


def murcko_scaffold_smiles(mol: Chem.Mol) -> str:
    """Return the SMILES of the Bemis-Murcko scaffold."""
    scaffold = murcko_scaffold(mol)
    return canonical_smiles(scaffold) if scaffold else ""


def maximum_common_substructure(mols: list[Chem.Mol], timeout: int = 5):
    """
    Compute the maximum common substructure (MCS) across a list of molecules.
    Returns an RDKit ``MCSResult`` or ``None``.
    """
    from rdkit.Chem import rdFMCS

    valid = [m for m in mols if m is not None]
    if len(valid) < 2:
        return None
    try:
        return rdFMCS.FindMCS(
            valid,
            timeout=timeout,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
        )
    except Exception as exc:  # noqa: BLE001
        logger.debug("MCS error: %s", exc)
        return None


def match_substructure_atoms(mol: Chem.Mol, smarts: str) -> list[int]:
    """Return atom indices in ``mol`` matching a SMARTS pattern."""
    if not smarts:
        return []
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        return []
    match = mol.GetSubstructMatch(patt)
    return list(match)


# ── 2D rendering (dark theme SVG) ──────────────────────────────────────────

# Palette shared with the Qt stylesheet so depictions blend into the UI.
_DARK_BG = (0.086, 0.106, 0.133)         # #161b22
_HIGHLIGHT = (0.122, 0.435, 0.922)       # #1f6feb (blue)
_HIGHLIGHT_2 = (0.247, 0.725, 0.314)     # #3fb950 (green)
_HIGHLIGHT_SOFT = (0.36, 0.78, 0.45)     # softer green for substructure tint


def render_svg(
    mol: Chem.Mol,
    width: int = 420,
    height: int = 320,
    highlight_atoms: Optional[list[int]] = None,
    dark: bool = True,
    legend: str = "",
) -> str:
    """
    Render a molecule to an SVG string with an optional set of highlighted
    atoms. Designed to drop straight into a ``QSvgWidget``.
    """
    if mol is None:
        return _empty_svg(width, height, dark)
    try:
        draw_mol = Chem.Mol(mol)
        try:
            AllChem.Compute2DCoords(draw_mol)
        except Exception:  # noqa: BLE001
            pass

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        if dark:
            opts.setBackgroundColour(_DARK_BG)
            _apply_dark_atom_palette(opts)
        opts.bondLineWidth = 2
        # Slightly tighter highlight so it reads as an outlined region, not blobs.
        try:
            opts.highlightRadius = 0.32
            opts.highlightBondWidthMultiplier = 12
        except Exception:  # noqa: BLE001 - older RDKit
            pass

        highlight_atoms = highlight_atoms or []
        highlight_bonds: list[int] = []
        if highlight_atoms:
            hset = set(highlight_atoms)
            for bond in draw_mol.GetBonds():
                if (
                    bond.GetBeginAtomIdx() in hset
                    and bond.GetEndAtomIdx() in hset
                ):
                    highlight_bonds.append(bond.GetIdx())
        atom_colors = {a: _HIGHLIGHT_SOFT for a in highlight_atoms}
        bond_colors = {b: _HIGHLIGHT_SOFT for b in highlight_bonds}

        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            draw_mol,
            highlightAtoms=highlight_atoms,
            highlightBonds=highlight_bonds,
            highlightAtomColors=atom_colors,
            highlightBondColors=bond_colors,
            legend=legend,
        )
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as exc:  # noqa: BLE001
        logger.warning("render_svg error: %s", exc)
        return _empty_svg(width, height, dark)


def _apply_dark_atom_palette(opts) -> None:
    """Brighten hetero-atom colours so they are legible on a dark background."""
    palette = {
        6: (0.83, 0.86, 0.85),    # C – light grey
        7: (0.35, 0.55, 1.0),     # N – blue
        8: (1.0, 0.42, 0.40),     # O – red
        9: (0.45, 0.85, 0.55),    # F – green
        15: (1.0, 0.62, 0.30),    # P – orange
        16: (0.95, 0.82, 0.30),   # S – yellow
        17: (0.45, 0.85, 0.55),   # Cl – green
        35: (0.80, 0.50, 0.35),   # Br
        53: (0.65, 0.45, 0.80),   # I
    }
    try:
        opts.updateAtomPalette(palette)
    except Exception:  # noqa: BLE001
        pass


def render_png(
    mol: Chem.Mol, width: int = 320, height: int = 240
) -> Optional[bytes]:
    """Render a molecule to PNG bytes (used in exports / tables when needed)."""
    if mol is None:
        return None
    try:
        img = Draw.MolToImage(mol, size=(width, height))
        import io

        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    except Exception as exc:  # noqa: BLE001
        logger.debug("render_png error: %s", exc)
        return None


def _empty_svg(width: int, height: int, dark: bool) -> str:
    """Return a placeholder SVG when no molecule is available."""
    bg = "#161b22" if dark else "#ffffff"
    fg = "#484f58"
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">'
        f'<rect width="100%" height="100%" fill="{bg}"/>'
        f'<text x="50%" y="50%" fill="{fg}" font-family="Segoe UI, Arial" '
        f'font-size="14" text-anchor="middle" dominant-baseline="middle">'
        f'No molecule loaded</text></svg>'
    )
