#! /usr/bin/env python
"""A module for tracking time usage."""

import time


#===============================================================================
class Timer(object):
    """Keeps track of time usage for different sections.

    Usage: Define the section names in the constructor dictionary, control time
    tracking with start() and stop().
    """

    def __init__(self, section_names):
        self.start_time = 0
        self.end_time = 0
        self.sections = {}
        for name in section_names:
            self.sections[name] = 0
        self.current_section = None

    def start(self, section_name):
        self.current_section = section_name
        self.start_time = time.time()

    def stop(self):
        self.end_time = time.time()
        elapsed = self.end_time - self.start_time
        self.sections[self.current_section] += elapsed
        self.current_section = None

    def get_total_time(self):
        return sum(self.sections.values())
