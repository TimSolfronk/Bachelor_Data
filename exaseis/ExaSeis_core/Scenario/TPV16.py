# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from .Scenario import Scenario

import os

class TPV16(Scenario):
    """
    Part of a series of benchmarks by the Statewide California Earthquake Center (SCEC)

    Spontaneous rupture on a vertical strike-slip fault in a homogeneous halfspace.
    There are slightly heterogeneous initial stress conditions. 

    The description of the scenario can be found at: https://strike.scec.org/cvws/tpv16_17docs.html
    """

    domain_offset   = [-5.6,  0., -5.6] #centers the fault which is at 20.
    domain_size     = [51.2, 51.2, 51.2] #leads to a cell size of 0.632 for 81 cells

    end_time = 15.0

    tracer_coordinates = [
          
    ]

    def initial_conditions(self):
        return """
      Q[Shortcuts::rho] = 2.67;
      Q[Shortcuts::cp ] = 6.0;
      Q[Shortcuts::cs ] = 3.464;
"""
    
    def generate_required_files(self, order):

        dictionary = { "MODE": order+1 }
        self.generate_file_from_template(
            os.path.dirname(os.path.realpath(__file__))+"/specs/tpv16/tpv16.yaml.template",
            "tpv16.yaml", dictionary
        )

        self.copy_file_to_current_folder(
            os.path.dirname(os.path.realpath(__file__))+"/specs/tpv16/tpv16_fault.yaml", 
            "tpv16_fault.yaml"
        )

        self.copy_file_to_current_folder(
            os.path.dirname(os.path.realpath(__file__))+"/specs/tpv16/tpv16_input_file.txt", 
            "tpv16_input_file.txt"
        )

        return