import time
import numpy as np
from ase.parallel import rank

class Timer(object):

    def __init__(self, record_time_usage, sections):
        self.start_time = 0
        self.end_time = 0
        self.record_time_usage = False
        self.sections = sections
        self.current_section = None

    def start(self, section_name):
        if rank == 0:
            self.current_section = section_name
            self.start_time = time.time()

    def end(self):
        if rank == 0:
            self.end_time = time.time()
            elapsed = self.end_time - self.start_time
            self.sections[self.current_section] += elapsed
            self.current_section = None
    
    def get_total_time(self):
        return np.sum(self.sections.values())

