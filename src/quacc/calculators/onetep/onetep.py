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
    def __init__(self, *args, max_walltime, **kwargs):
        super().__init__(*args, **kwargs)

        self.atoms_list = []

        self.max_walltime = max_walltime
        self.created_time = time.time()

    def read_results(self, directory):
        try:
            output_path = directory / self.output
            atoms = read(output_path, format="onetep-out")
            self.atoms_list.append(atoms)
            return dict(atoms.calc.properties())
        except Exception as e:
            raise e(directory, self.label, self.atoms_list) from e

    def execute(self, directory, profile):
        try:
            profile.run(directory, self.input, self.output, self.error, self.append)
        except Exception as e:
            raise e(directory, self.label, self.atoms_list) from e


class Onetep(Onetep_):
    def __init__(
        self,
        input_atoms: Atoms = None,
        profile: OnetepProfile = None,
        parallel_info: dict[str | Any] | None = None,
        **kwargs,
    ):
        self.input_atoms = input_atoms

        self._binary = str(SETTINGS.ONETEP_CMD)

        profile = OnetepProfile(self._binary, parallel_info=parallel_info)

        super().__init__(profile=profile, **kwargs)
