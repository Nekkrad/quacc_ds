"""Exceptions for Quacc"""


class QuaccException(Exception):
    def __init__(self, job_error=None, read_error=None, current_state=None):
        self.job_error = job_error
        self.read_error = read_error
        self.current_state = current_state

    def __str__(self):
        if self.job_error:
            return f"failed when running {self.job_error}"
        if self.read_error:
            return f"failed when reading {self.read_error}"
