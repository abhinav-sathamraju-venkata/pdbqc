import argparse
from .structure import clean_structure, format_report

def main():
    parser = argparse.ArgumentParser(description="Run structure QC on a PDB or mmCIF file")
    parser.add_argument("path", help="Path to PDB or mmCIF file")
    parser.add_argument("--cif", default="", help="Path to mmCIF file for SEQRES data")
    parser.add_argument("--keep-ligands", action="store_true", help="Keep ligands in cleaned structure")
    args = parser.parse_args()

    qc, structure = clean_structure(args.path, keep_ligands=args.keep_ligands, cif_path=args.cif)
    format_report(qc)

if __name__ == "__main__":
    main()