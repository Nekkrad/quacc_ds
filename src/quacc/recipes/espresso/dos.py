"""
This module, 'dos.py', contains recipes for performing phonon calculations using the
dos.x binary from Quantum ESPRESSO via the quacc library.

The recipes provided in this module are jobs and flows that can be used to perform
dos calculations.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from quacc import flow, job
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.recipes.espresso._base import base_fn
from quacc.recipes.espresso.core import non_scf_job, static_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema


@job
def dos_job(
    prev_dir: str | Path,
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic dos.x calculation (density of states).
    It is mainly used to extract the charge density and wavefunction from a previous pw.x calculation.
    It generates the total density of states. Fore more details please see
    https://www.quantum-espresso.org/Doc/INPUT_DOS.html

    Parameters
    ----------
    prev_dir
        Outdir of the previously ran pw.x calculation. This is used to copy
        the entire tree structure of that directory to the working directory
        of this calculation.
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
    **calc_kwargs
        calc_kwargs dictionary possibly containing the following keys:

        - input_data: dict
        - additional_fields: list[str] | str

        See the docstring of ase.io.espresso.write_fortran_namelist for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "input_data": {
            "DOS": {"fileout": "total_dos.dos"},
            },
        }

    return base_fn(
        template=EspressoTemplate("dos", test_run=test_run),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "dos.x Density-of-States"},
        copy_files=prev_dir,
    )


@flow
def dos_flow(
    atoms: Atoms,
    job_decorators: dict[str, Callable | None] | None = None,
    job_params: dict[str, Any] | None = None,
) -> RunSchema:
    """
    This function performs a total density of states calculations.

    Consists of following jobs that can be modified:

    1. pw.x static
        - name: "static_job"
        - job: [quacc.recipes.espresso.core.static_job][]

    2. pw.x non self-consistent
        - name: "non_scf_job"
        - job: [quacc.recipes.espresso.core.non_scf_job][]

    3. dos.x total density of states
        - name: "ph_job"
        - job: [quacc.recipes.espresso.dos.dos_job][]

    Parameters
    ----------
    atoms
        Atoms object
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """
    calc_defaults = {
        "static_job": {
            "kspacing":0.05,
            "input_data":{
                "system":{
                    "occupations":"tetrahedra"
                }
            }
            },
        "non_scf_job": {
            "kspacing":0.05,
            "input_data":{
                "system":{
                    "occupations":"tetrahedra"
                }
            }
        },
        "dos_job": {
            "input_data": {
                "DOS": {
                    "fileout": "total_dos.dos",
                    "bz_sum":"tetrahedra"
                    }
            }
        },
        }
    job_params = recursive_dict_merge(calc_defaults, job_params)

    pw_job, nscf_job, tdos_job= customize_funcs(
        ["static_job", "non_scf_job", "dos_job"],
        [static_job, non_scf_job, dos_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    return [pw_job(atoms),nscf_job(atoms),tdos_job(atoms)]
