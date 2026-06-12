import numpy as np
import math
from collections import defaultdict
import biotite.structure.io as strucio
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
from rich.console import Console
from rich.panel import Panel
from rich import box
from rich.table import Table

def clean_structure (path, keep_ligands = False, cif_path = ""):

    #Step 1: Filter to model 1 only
    structure = strucio.load_structure(path, altloc="all", extra_fields=["occupancy", "b_factor"])
    if isinstance(structure, struc.AtomArrayStack):
        structure = structure[0]





    altloc_count = int(np.sum(structure.altloc_id != "."))
    insertion_codes = sorted(set(
        f"{structure.chain_id[i]}{structure.res_id[i]}{structure.ins_code[i]}"
        for i in range(structure.array_length())
        if structure.ins_code[i] != ""
    ))


    #Step 2: Report any Microheterogenity
    micro = {}
    residue_names = defaultdict(set)
    for i in range (structure.array_length()):
        key = (structure.chain_id[i], structure.res_id[i], structure.ins_code[i])
        residue_names[key].add(structure.res_name[i])

        if len(residue_names[key]) > 1:
            micro[key] = residue_names[key]
    
    #Step 3: Filter to altloc's with the highest possible occupancy
    structure = structure[struc.filter_highest_occupancy_altloc(structure, structure.altloc_id, structure.occupancy)].copy()

    
    #Step 4: Set it so that if keep_ligands == False, then all non-protein molecules are filtered out
    if keep_ligands == False:
        structure = structure[~structure.hetero]

    #Step 5: Filter all MSE into MET for standardized naming convention
    mse_mask = structure.res_name == "MSE"
    structure.res_name[mse_mask] = "MET"

    #Step 6: Find missing residues
    if cif_path != "":
        import biotite.structure.io.pdbx as pdbx
        cif_file = pdbx.CIFFile.read(cif_path)
        seqres = pdbx.get_sequence(cif_file)
        missing = {}
        coverage = {}
        all_chain_ids = np.unique(structure.chain_id)

        for chain_id in all_chain_ids:
            if chain_id in seqres:
                ca_atoms = structure[(structure.atom_name == "CA") & (structure.chain_id == chain_id)]    
                n_observed = len(ca_atoms)
                n_expected = len(seqres[chain_id])
                missing[chain_id] = n_expected - n_observed
                coverage[chain_id] = round(n_observed / n_expected * 100, 1) if n_expected > 0 else 0.0
                       
    else:
        missing = {}
        coverage = {}

    


    mse_count = int(np.sum(mse_mask))


    qc = {
    "missing_residues": missing,
    "microheterogeneities": micro,
    "n_atoms": structure.array_length(),
    "altloc_count": altloc_count,
    "insertion_codes": insertion_codes,
    "mse_count": mse_count,
    "coverage" : coverage
    }

    return qc, structure





console = Console()

def format_report(qc):
    console.print(Panel("[bold]Structure Quality Report[/bold]", box=box.DOUBLE))
    
    # Basic stats
    console.print(f"[cyan]Atoms:[/cyan] {qc['n_atoms']}")
    
    # Altlocs
    if qc['altloc_count'] > 0:
        console.print(f"[yellow]Altlocs detected:[/yellow] {qc['altloc_count']} atoms")
    else:
        console.print(f"[green]Altlocs detected:[/green] none")
    
    # Insertion codes
    if qc['insertion_codes']:
        console.print(f"[yellow]Insertion codes:[/yellow] {', '.join(qc['insertion_codes'])}")
    else:
        console.print(f"[green]Insertion codes:[/green] none")
    
    # Microheterogeneity
    if qc['microheterogeneities']:
        for key, names in qc['microheterogeneities'].items():
            chain, res_id, ins = key
            console.print(f"[yellow]Microheterogeneity:[/yellow] {chain}:{res_id} ({'/'.join(sorted(names))})")
    else:
        console.print(f"[green]Microheterogeneity:[/green] none")
    
    # Missing residues + coverage
    if qc['missing_residues']:
        table = Table(box=box.SIMPLE, header_style="bold cyan")
        table.add_column("Chain")
        table.add_column("Missing Residues")
        table.add_column("Coverage")
        for chain, count in qc['missing_residues'].items():
            cov = qc['coverage'].get(chain, "N/A")
            color = "green" if count == 0 else "yellow" if count < 10 else "red"
            table.add_row(chain, f"[{color}]{count}[/{color}]", f"{cov}%")
        console.print(table)
    else:
        console.print("[dim]Missing residues: no CIF provided[/dim]")
    
    # MSE
    if qc['mse_count'] > 0:
        console.print(f"[yellow]MSE residues converted:[/yellow] {qc['mse_count']}")
    else:
        console.print(f"[green]MSE residues converted:[/green] none")