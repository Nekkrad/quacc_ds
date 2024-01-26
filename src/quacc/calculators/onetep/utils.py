from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.io.espresso import Namelist

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms



def onetep_copy_files(
    prev_dir: str | Path, cp_files: list[str]
) -> dict[str, list[str | Path]]:
    """
    Function that take care of copying the correct files from a previous ONETEP
    to a current ONETEP calculation
    Parameters
    ----------
    prev_dir
        Outdir of the previously ran onetep calculation. This is used to partially
        copy the tree structure of that directory to the working directory
        of this calculation.
    cp_files
        List of files to copy
    Returns
    -------
    dict
        Dictionary of files to copy. The key is the directory to copy from and
        the value is a list of files to copy from that directory.
    """

    file_to_copy = {prev_dir: []}


    for file in cp_files:
        file_to_copy[prev_dir].append(file)


    return file_to_copy