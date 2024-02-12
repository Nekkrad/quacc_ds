from __future__ import annotations

import time
from typing import TYPE_CHECKING

from ase.calculators.onetep import Onetep as Onetep_
from ase.calculators.onetep import OnetepProfile
from ase.calculators.onetep import OnetepTemplate as OnetepTemplate_
from ase.io import read

from quacc import SETTINGS

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms


class OnetepTemplate(OnetepTemplate_):
    def __init__(self, append, timeout, **kwargs):
        super().__init__(append=append, timeout=timeout)

        self.created_time = time.time()
        self.max_walltime = timeout

        self.job_error = None
        self.read_error = None

    def read_results(self, directory):
        try:
            atoms = read(directory / self.outputname, format="onetep-out")

            return self._error_handler(dict(atoms.calc.properties()))
        except Exception as e:
            self.read_error = e

            return self._error_handler({"energy": None})

    def _error_handler(self, results):
        errors = {"job_error": self.job_error, "read_error": self.read_error}
        return {**results, **errors}

    def execute(self, directory, profile):
        corrected_time = self.max_walltime - 300

        try:
            profile.run(
                directory,
                self.inputname,
                self.outputname,
                self.errorname,
                append=self.append,
                timeout=corrected_time,
            )
        except Exception as e:
            self.job_error = e


class Onetep(Onetep_):
    def __init__(
        self,
        input_atoms: Atoms = None,
        parallel_info: dict[str | Any] | None = None,
        opt_restart: bool = False,
        **kwargs,
    ):
        self.input_atoms = input_atoms

        self._binary = str(SETTINGS.ONETEP_CMD)
        profile = OnetepProfile(self._binary, parallel_info=parallel_info)

        super().__init__(profile=profile, **kwargs)

        self.template = OnetepTemplate(append=True, timeout=SETTINGS.WALLTIME)

        if opt_restart:
            try:
                self.atoms = self.input_atoms.copy()
                self.results = self.input_atoms.calc.results.copy()
            except AttributeError:
                self.atoms = None
                self.results = {}
