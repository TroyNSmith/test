import re
from pathlib import Path

import automol
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdchem import Mol
from sqlalchemy import func
from sqlmodel import Field, Relationship, Session, SQLModel, create_engine, select

from inputs import (
    CALC_INP,
    FREQ_INP,
    GOAT_INP,
    OPT_INP,
    SCAN_INP,
    STATIONARY_SH,
    TRANS_FREQ_INP,
    TRANSITION_SH,
)

# === Constants ===
A_ = "Cv4,Ov2"
C_ = "Cv4"
O_ = "Ov2"
As = "CX4,OX2"
Cs = "CX4"
Os = "OX2"
Ar = "Cv3,Ov1"
Cr = "Cv3"
Or = "Ov1"
Au = "CX3,OX1"
Cu = "CX3"
Ou = "OX1"

RADICAL_TOKENS = [Ar, Cr, Or]


# === SQL ===
class ReactionLink(SQLModel, table=True):
    stationary_id: int | None = Field(
        default=None, foreign_key="stationary.id", primary_key=True
    )
    transition_id: int | None = Field(
        default=None, foreign_key="transition.id", primary_key=True
    )
    role: str  # NOTE: "reactant", "product"


class Stationary(SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    amchi: str = Field(index=True, unique=True)
    smiles: str
    directory_name: str

    transitions: list["Transition"] = Relationship(
        back_populates="stationary_points", link_model=ReactionLink
    )


class Transition(SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    amchi: str = Field(index=True, unique=True)
    scan_string: str
    active_atoms: str
    directory_name: str

    stationary_points: list["Stationary"] = Relationship(
        back_populates="transitions", link_model=ReactionLink
    )


class SQLReactionRegistry:
    def __init__(self, db_url="sqlite://reactions.db", root_dir="../results"):
        self.engine = create_engine(db_url)
        SQLModel.metadata.create_all(self.engine)
        self.root = Path(root_dir)
        self.root.mkdir(exist_ok=True)

    def get_or_add_stationary(self, amchi: str, smiles: str, xyz: str):
        with Session(self.engine) as session:
            statement = select(Stationary).where(Stationary.amchi == amchi)
            existing = session.exec(statement).first()

            if existing:
                return existing.id

            count_statement = select(func.count()).select_from(Stationary)
            count = session.exec(count_statement).one() + 1
            dir_name = f"S{count}"

            new_s = Stationary(amchi=amchi, smiles=smiles, directory_name=dir_name)
            session.add(new_s)
            session.commit()

            self._write_files(dir_name, xyz, "stationary")
            return new_s.id

    def get_or_add_transition(
        self,
        amchi: str,
        xyz: str,
        scan: str,
        atoms: str,
        r_ids: list[int],
        p_ids: list[int],
    ):
        with Session(self.engine) as session:
            statement = select(Transition).where(Transition.amchi == amchi)
            existing = session.exec(statement).first()

            if existing:
                return existing.id

            count_statement = select(func.count()).select_from(Transition)
            count = session.exec(count_statement).one() + 1
            dir_name = f"T{count}"

            new_t = Transition(
                amchi=amchi,
                scan_string=scan,
                active_atoms=atoms,
                directory_name=dir_name,
            )
            session.add(new_t)
            session.flush()

            for s_id in r_ids:
                stat = session.get(Stationary, s_id)
                if stat:
                    link = ReactionLink(
                        stationary_id=s_id, transition_id=new_t.id, role="reactant"
                    )
                    session.add(link)

            for s_id in p_ids:
                stat = session.get(Stationary, s_id)
                if stat:
                    link = ReactionLink(
                        stationary_id=s_id, transition_id=new_t.id, role="product"
                    )
                    session.add(link)

            session.add(new_t)
            session.commit()

            self._write_files(dir_name, xyz, "transition", scan, atoms)
            return new_t.id

    def _write_files(self, dir_name, xyz, role, scan=None, active_atoms=None):
        path = self.root / dir_name
        path.mkdir(exist_ok=True, parents=True)
        (path / "guess.xyz").write_text(xyz)
        if role == "stationary":
            (path / "goat.inp").write_text(GOAT_INP)
            (path / "opt.inp").write_text(OPT_INP)
            (path / "freq.inp").write_text(FREQ_INP)
            (path / "calc.inp").write_text(CALC_INP)
            (path / "submit.sh").write_text(STATIONARY_SH)

        if role == "transition":
            TF_INP = TRANS_FREQ_INP.replace("[ATOMS]", active_atoms)
            SC_INP = SCAN_INP.replace("[SCAN]", scan).replace("[ATOMS]", active_atoms)
            (path / "scan.inp").write_text(SC_INP)
            (path / "freq.inp").write_text(TF_INP)
            (path / "calc.inp").write_text(CALC_INP)
            (path / "submit.sh").write_text(TRANSITION_SH)


# === RDKit ===
def mol_from_smiles(smiles: str, with_coords: bool = False) -> Mol:
    """
    Generate an RDKit Mol from a SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    if with_coords:
        raise NotImplementedError

    return mol


def mol_to_smiles(molecule: Mol, ignore_map_numbers: bool = True) -> str:
    """Converts RDKit molecule(s) to SMILES string(s)."""

    if ignore_map_numbers:
        molecule = Chem.RemoveHs(molecule)
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)

        Chem.rdmolops.AssignRadicals(molecule)

    return Chem.MolToSmiles(molecule, ignoreAtomMapNumbers=True)


def smarts_radicals(smarts: str) -> set[int]:
    """Deduce radical electron positions from SMARTS string."""
    radicals = set()

    atom_re = re.compile(r"\[([^\]:]+):(\d+)\]")
    for valence_identity, map_number in atom_re.findall(smarts):
        if valence_identity in ["Cv3,Ov1", "Cv3", "Ov1"]:
            radicals.add(int(map_number))

    return radicals


def isomorphic(molecules_1: Mol, molecules_2: Mol) -> bool:
    """Determines whether two sets of molecules are isomorphic to each other."""
    matcher = graph_matcher(molecules_1, molecules_2)
    return matcher.is_isomorphic()


def graph_matcher(molecules_1: list[Mol], molecules_2: list[Mol]):
    """Determines isomorphism of one molecule onto another."""
    graph_1 = molecular_graph(molecules_1)
    graph_2 = molecular_graph(molecules_2)

    def _node_match(node_1, node_2) -> bool:
        return node_1["symbol"] == node_2["symbol"]

    return GraphMatcher(graph_1, graph_2, node_match=_node_match)


def molecular_graph(mols: list[Mol]) -> nx.Graph:
    """Generates a NetworkX graph from RDKit molecule."""
    graph = nx.Graph()

    for idx, mol in enumerate(mols):
        for atom in mol.GetAtoms():
            key = (idx, atom.GetIdx())
            graph.add_node(key, symbol=atom.GetSymbol())

        for bond in mol.GetBonds():
            key_1 = (idx, bond.GetBeginAtomIdx())
            key_2 = (idx, bond.GetEndAtomIdx())
            graph.add_edge(key_1, key_2)

    return graph


def unique_molecules(molecules: list[Mol]) -> list[Mol]:
    "Identifies unique molecules from a sequence."
    unique_set = []

    for molecule in molecules:
        if not any(isomorphic(molecule, m) for m in unique_set):
            unique_set.append(molecule)

    return unique_set


def reaction(reactants: tuple[Mol], smarts: str, isomorphs: bool = False) -> list[Mol]:
    """
    Perform a reaction on reactant(s) following SMARTS template.
    """
    rxn = rdChemReactions.ReactionFromSmarts(smarts)

    lhs_smarts, rhs_smarts = smarts.split(">>")
    rhs_radicals = smarts_radicals(rhs_smarts)

    if isinstance(reactants, Mol):
        reactants = (reactants,)

    for reactant in reactants:
        for atom in reactant.GetAtoms():
            atom.SetIntProp("molAtomMapNumber", atom.GetIdx())

    product_sets = list(rxn.RunReactants(reactants))
    product_sets = product_sets if isomorphs else unique_molecules(product_sets)

    for products in product_sets:
        for product in products:
            for atom in product.GetAtoms():
                if atom.HasProp("old_mapno"):
                    if int(atom.GetProp("old_mapno")) in rhs_radicals:
                        atom.SetNumRadicalElectrons(1)

    return product_sets


# === AutoMol ===
def reaction_graphs(
    reactants: str | list[str],
    products: str | list[str],
):
    """Generates a reaction Automol graph from AMChI strings."""
    reactants = (reactants,) if isinstance(reactants, str) else reactants
    products = (products,) if isinstance(products, str) else products

    return automol.reac.from_amchis(reactants, products)


def stationary_graph(smiles: str = None, amchi: str = None, canonical: str = True):
    """Generates a stationary AutoMol graph from SMILES or AMChI string."""
    assert smiles is not None or amchi is not None, (
        "A SMILES or AMChI string must be provided to create a graph."
    )

    graph = (
        automol.smiles.graph(smiles)
        if smiles is not None
        else automol.amchi.graph(amchi)
    )
    if canonical:
        return automol.graph.canonical(graph)

    return graph


def transition_graph(reaction, canonical: bool = True):
    """Generates a transition AutoMol graph from AutoMol reaction."""
    graph = automol.reac.ts_graph(reaction)

    if canonical:
        return automol.graph.canonical(graph)

    return graph


def canonical_enantiomer(graph):
    """Returns canonical amchi and graph of stationary species."""
    amchi = automol.graph.amchi(graph)

    if not automol.amchi.is_canonical_enantiomer(amchi):
        amchi = automol.amchi.canonical_enantiomer(amchi)
        graph = stationary_graph(amchi=amchi)

    return amchi, graph


def process_rdkit_reaction(reactant: Mol, product_sets: list[Mol]):
    """Processes RDKit reactions through the AutoMol package."""
    registry = SQLReactionRegistry(
        db_url="sqlite:///species.db", root_dir="results"
    )
    reactant_smi = mol_to_smiles(reactant)
    reactant_amchi, reactant_graph = canonical_enantiomer(
        stationary_graph(reactant_smi)
    )
    reactant_geo = automol.graph.geometry(reactant_graph)
    reactant_xyz = automol.geom.xyz_string(reactant_geo)

    reactant_id = registry.get_or_add_stationary(
        amchi=reactant_amchi, smiles=reactant_smi, xyz=reactant_xyz
    )

    for products in product_sets:
        product_amchis, product_ids = [], []
        for product in products:
            product_smi = mol_to_smiles(product)
            product_amchi, product_graph = canonical_enantiomer(
                stationary_graph(product_smi)
            )
            product_amchis.append(product_amchi)

            product_geo = automol.graph.geometry(product_graph)
            product_xyz = automol.geom.xyz_string(product_geo)

            product_id = registry.get_or_add_stationary(
                amchi=product_amchi, smiles=product_smi, xyz=product_xyz
            )
            product_ids.append(product_id)

            for reaction in reaction_graphs(reactant_amchi, tuple(product_amchis)):
                transition_gra = transition_graph(reaction)
                try:
                    transition_amchi = automol.graph.amchi(transition_gra)

                    with Session(registry.engine) as session:
                        statement = select(Transition).where(
                            Transition.amchi == transition_amchi
                        )
                        if session.exec(statement).first():
                            continue

                    transition_geo = automol.graph.geometry(transition_gra)
                    transition_xyz = automol.geom.xyz_string(transition_geo)

                    formed = automol.graph.ts.forming_bond_keys(transition_gra)
                    broken = automol.graph.ts.breaking_bond_keys(transition_gra)

                    dmat_angstrom = (
                        automol.geom.distance_matrix(transition_geo) * 0.529177
                    )

                    for broken_bond in broken:
                        a, b = broken_bond
                        if len(formed) > 0:
                            for formed_bond in formed:
                                shared = broken_bond & formed_bond
                                shared_idx = next(iter(shared))
                                if len(shared) == 1:
                                    c = next(iter(formed_bond - shared))
                                    dist = dmat_angstrom[shared_idx, c]
                                    active_atoms = f"{shared_idx} {c}"
                                    scan = (
                                        f"scan B {active_atoms} = {dist:.3f}, 0.7, 100"
                                    )

                        else:
                            active_atoms = f"{a} {b}"
                            scan = f"scan B {active_atoms} = {dmat_angstrom[a, b]:.3f}, 2.5, 100"

                        registry.get_or_add_transition(
                            transition_amchi,
                            transition_xyz,
                            scan,
                            active_atoms,
                            [
                                reactant_id,
                            ],
                            product_ids,
                        )

                except Exception:
                    continue