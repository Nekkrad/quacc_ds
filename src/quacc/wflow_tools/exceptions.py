"""Exceptions for Quacc"""


class QuaccException(Exception):
    def __init__(self, job_error=None, read_error=None, current_atoms=None):
        self.job_error = job_error
        self.read_error = read_error
        self.current_atoms = current_atoms
