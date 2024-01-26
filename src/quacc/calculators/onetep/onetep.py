from typing import TYPE_CHECKING

from ase import Atoms
from ase.calculators.onetep import Onetep as Onetep_
from ase.calculators.onetep import OnetepProfile,OnetepTemplate

from quacc import SETTINGS
from quacc.utils.dicts import recursive_dict_merge
from typing import Any
if TYPE_CHECKING:
    from typing import Any


class Onetep(Onetep_):
    def __init__(
        self,
        input_atoms: Atoms = None,
        directory: str = ".",
        profile: OnetepProfile = None,
        template: OnetepTemplate =None,
        calc_defaults: dict[str | Any ] | None = None,
        parallel_info: dict[str | Any ]| None = None,
        **kwargs,
    ):
        profile = OnetepProfile(str(SETTINGS.ONETEP_CMD), parallel_info=parallel_info)
        self.input_atoms = input_atoms
        self.calc_defaults = calc_defaults
        self.directory =directory
        self.template = template
        kwargs = recursive_dict_merge(self.calc_defaults, kwargs)

        super().__init__(
            profile=profile, template=template,directory=directory, parallel_info=parallel_info, **kwargs
        )