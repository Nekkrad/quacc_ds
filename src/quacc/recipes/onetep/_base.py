"""Base jobs for Onetep."""
from __future__ import annotations

from typing import TYPE_CHECKING
from ase import Atoms
from quacc.calculators.onetep.onetep import Onetep, OnetepProfile
from ase.calculators.onetep import OnetepTemplate


from quacc import SETTINGS
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas._aliases.ase import RunSchema


def base_fn(
    atoms: Atoms = Atoms(),
    directory: str = ".",
    template: OnetepTemplate | None = None,
    profile: OnetepProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    parallel_info: dict[str] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:

    """
    Base function to carry out Onetep recipes.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the ONETEP calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.onetep.Onetep` calculator.
    opt_defaults
        The default optimization parameters.
    opt_params
        Dictionary of parameters to pass to the optimizer. pass "optimizer"
        to change the optimizer being used. "fmax" and "max_steps" are commonly
        used keywords. See the ASE documentation for more information.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)
    atoms.calc = Onetep(
        input_atoms=atoms,
        directory=directory,
        template=template,
        profile=profile,
        calc_defaults=calc_defaults,
        parallel_info=parallel_info,
        pseudo_path=str(SETTINGS.ONETEP_PP_PATH) if SETTINGS.ONETEP_PP_PATH else ".",
        **calc_flags,
    )


    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(final_atoms, atoms, additional_fields=additional_fields)


def base_opt_fn(
    atoms: Atoms = Atoms(),
    directory : str =".",
    template: OnetepTemplate | None = None,
    profile: OnetepProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: dict[str, Any] | None = None,
    parallel_info: dict[str] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:

    """
    Base function to carry out Onetep recipes.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the ONETEP calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.onetep.Onetep` calculator.
    opt_defaults
        The default optimization parameters.
    opt_params
        Dictionary of parameters to pass to the optimizer. pass "optimizer"
        to change the optimizer being used. "fmax" and "max_steps" are commonly
        used keywords. See the ASE documentation for more information.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    atoms.calc = Onetep(
        input_atoms=atoms,
        directory=directory,
        template=template,
        profile=profile,
        calc_defaults=calc_defaults,
        parallel_info=parallel_info,
        pseudo_path=str(SETTINGS.ONETEP_PP_PATH) if SETTINGS.ONETEP_PP_PATH else ".",
        **calc_flags,
    )

    final_atoms = run_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(final_atoms, additional_fields=additional_fields)
