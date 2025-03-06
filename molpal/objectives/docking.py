import atexit
import dataclasses
import csv
from typing import Dict, Iterable, Optional

import numpy as np
import pyscreener as ps

from molpal.objectives.base import Objective


class DockingObjective(Objective):
    """A DockingObjective calculates the objective function by calculating the
    docking score of a molecule

    Attributes
    ----------
    c : int
        the min/maximization constant, depending on the objective
    virtual_screen : pyscreener.docking.DockingVirtualScreen
        the VirtualScreen object that calculated docking scores of molecules against a given
        receptor with specfied docking parameters

    Parameters
    ----------
    objective_config : str
        the path to a pyscreener config file containing the options for docking calculations
    path : str, default="."
        the path under which docking inputs/outputs should be collected
    verbose : int, default=0
        the verbosity of pyscreener
    minimize : bool, default=True
        whether this objective should be minimized
    **kwargs
        additional and unused keyword arguments
    """

    def __init__(
        self,
        objective_config: str,
        path: str = ".",
        verbose: int = 0,
        minimize: bool = True,
        **kwargs,
    ):
        args = ps.args.gen_args(f"--config {objective_config}")

        #Support for DOCK3.8 on the UCSF cluster
        if args.screen_type in ["dock3", "dock3.8"]:
            self.virtual_screen = ps.DOCK3VirtualScreen(
                output_dir = args.docking_output_dir,
                screen_type = args.screen_type,
                dockfiles = args.dockfiles,
                pipeline_scripts = args.pipeline_scripts,
                ncpu = args.ncpu,
                library_file = args.library,
            )
        else:
            metadata_template = ps.build_metadata(args.screen_type, args.metadata_template)
            self.virtual_screen = ps.virtual_screen(
                args.screen_type,
                args.receptors,
                args.center,
                args.size,
                metadata_template,
                args.pdbids,
                args.docked_ligand_file,
                args.buffer,
                args.ncpu,
                args.base_name,
                path,
                args.reduction,
                args.receptor_reduction,
                args.k,
                args.dockfiles,
                args.pipeline_scripts,
            )
            # Print all arguments as a dictionary
        print("\n=== All Arguments ===")
        args_dict = vars(args)
        for key, value in args_dict.items():
            print(f"{key}: {value}")
        print("==================\n")
        atexit.register(self.cleanup)

        super().__init__(minimize=minimize)

    def forward(self, smis: Iterable[str], **kwargs) -> Dict[str, Optional[float]]:
        """Calculate the docking scores for a list of SMILES strings

        Parameters
        ----------
        smis : List[str]
            the SMILES strings of the molecules to dock
        **kwargs
            additional and unused positional and keyword arguments

        Returns
        -------
        scores : Dict[str, Optional[float]]
            a map from SMILES string to docking score. Ligands that failed
            to dock will be scored as None
        """
        Y = self.c * self.virtual_screen(smis)
        Y = np.where(np.isnan(Y), None, Y)

        # Only do detailed matching for DOCK3.8
        if self.virtual_screen.screen_type in ["dock3", "dock3.8"]:
            # Get the stored results which contain exact SMILES-ZINC_ID-score mappings
            docking_results = self.virtual_screen._results
            
            # Create mapping of SMILES to scores using stored results
            smiles_to_score = {}
            for result in docking_results:
                smiles_to_score[result.smiles] = self.c * result.score
            
            # Map each input SMILES to its score, using None for unmatched SMILES
            return {smi: smiles_to_score.get(smi, None) for smi in smis}
        else:
            # For other docking methods, use simple zip
            return dict(zip(smis, Y))

    def cleanup(self):
        results = self.virtual_screen.results()
        self.virtual_screen.collect_files()

        with open(self.virtual_screen.path / "extended.csv", "w") as fid:
            writer = csv.writer(fid)
            writer.writerow(field.name for field in dataclasses.fields(results[0]))
            writer.writerows(dataclasses.astuple(r) for r in results)
