# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from .Scenario import Scenario

import os

class TPV8(Scenario):
    """
    Part of a series of benchmarks by the Statewide California Earthquake Center (SCEC)

    Spontaneous rupture on a vertical strike-slip fault in a homogeneous halfspace. 
    Initial stress conditions are linearly dependent on depth. Subshear rupture conditions. 

    The description of the scenario can be found at: https://strike.scec.org/cvws/tpv89docs.html
    """

    domain_offset   = [3.8, 0.0, 3.8] #centers the fault which is at 20.
    domain_size     = [32.4, 32.4, 32.4] #leads to a cell size of 0.333 for 81 cells

    end_time = 15.0

    tracer_sets = {
        "l3_o5": [
            # off-fault at surface
            [17.0, 0.0, 20.0],
            [18.0, 0.0, 20.0],
            [19.0, 0.0, 20.0],
            [21.0, 0.0, 20.0],
            [22.0, 0.0, 20.0],
            [23.0, 0.0, 20.0],
            [17.0, 0.0, 32.0],
            [23.0, 0.0, 32.0],
            # off-fault at 0.3km depth
            [18.9, 0.3, 20.0],
            [19.4, 0.3, 20.0],
            [20.6, 0.3, 20.0],
            [21.1, 0.3, 20.0],

            # off-fault at hypocenter depth
            [17.0,12.0, 32.0],
            [23.0,12.0, 32.0],

            # near-side on-fault
            # on-fault at surface
            [19.952, 0.0, 20.0],
            [19.952, 0.0, 24.5],
            [19.952, 0.0, 32.0],
            #on-fault below surface
            [19.952, 4.5, 20.0], 
            [19.952, 7.5, 20.0],
            [19.952, 7.5, 24.5],
            [19.952, 7.5, 32.0],
            [19.952,12.0, 20.0], #center of fault hypocenter

            # far-side on-fault
            # on-fault at surface
            [19.972, 0.0, 20.0],
            [19.972, 0.0, 24.5],
            [19.972, 0.0, 32.0],
            #on-fault below surface
            [19.972, 4.5, 20.0], 
            [19.972, 7.5, 20.0],
            [19.972, 7.5, 24.5],
            [19.972, 7.5, 32.0],
            [19.972,12.0, 20.0], #center of fault hypocenter
        ]
    }

    tracer_coordinates = tracer_sets["l3_o5"]

    def initial_conditions(self):
        return """
      Q[Shortcuts::rho] = 2.7;
      Q[Shortcuts::cp ] = 5.716;
      Q[Shortcuts::cs ] = 3.3;
"""
    
    def generate_required_files(self, order):

        dictionary = { "MODE": order+1 }
        self.generate_file_from_template(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV8/tpv8.yaml.template",
            "tpv8.yaml", dictionary
        )

        self.copy_file_to_current_folder(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV8/tpv8_fault.yaml", 
            "tpv8_fault.yaml"
        )

        return
