import pyscreener as ps
from pathlib import Path
import numpy as np
from typing import List
def load_smiles(file_path: str) -> List[str]:
    """Load SMILES from a text file"""
    with open(file_path) as f:
        return [line.strip() for line in f]

def test_dock3_virtual_screen():
    """Test DOCK3VirtualScreen with a set of SMILES"""
    print("Starting test_dock3_virtual_screen")
    # Initialize virtual screen
    virtual_screen = ps.DOCK3VirtualScreen(
        screen_type="dock3.8",
        dockfiles="/nfs/home/rmaksymiuk/.conda/envs/molpal/lib/python3.8/site-packages/pyscreener/docking/dock3/running_scripts/dockfiles",
        library_file="/nfs/home/rmaksymiuk/molpal/data/libraries/nsp3_mac1/H24P060_Mac1_Everted_ZINC22_docking_results.csv",
        output_dir="/nfs/home/rmaksymiuk/molpal/tests/output_virtual_screen_test"
    )

    # Load test SMILES
    smiles_file = "smiles_76k.txt"
    smis = load_smiles(smiles_file)
    print(f"Loaded {len(smis)} SMILES strings")

    # Run docking
    try:
        scores = virtual_screen(smis)
        print(f"Docking completed. Got {len(scores)} scores")
        
        # Print results
        results = virtual_screen.results()
        print(f"\nDetailed results ({len(results)} compounds):")
        for result in results[:5]:
            print(f"SMILES: {result.smiles[:30]}... ZINC ID: {result.zinc_id} Score: {result.score:.2f}")
        
        # Copy output files
        virtual_screen.collect_files()
        print(f"\nOutput files saved to: {virtual_screen.path}")
        
    except Exception as e:
        print(f"Error during docking: {e}")

def main():
    test_dock3_virtual_screen()

if __name__ == "__main__":
    main()