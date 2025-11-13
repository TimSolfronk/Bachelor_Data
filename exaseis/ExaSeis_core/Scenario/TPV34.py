# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from .Scenario import Scenario

import os

class TPV5(Scenario):
    """
    Part of a series of benchmarks by the Statewide California Earthquake Center (SCEC)

    Spontaneous rupture on a vertical strike-slip fault in a homogeneous halfspace.
    There are slightly heterogeneous initial stress conditions. 

    The description of the scenario can be found at: https://strike.scec.org/cvws/tpv5docs.html
    """

    domain_offset   = [6.5,  0., 6.5] #centers the fault which is at 20.
    domain_size     = [27., 27., 27.] #leads to a cell size of 0.333 for 81 cells

    end_time = 10.0

    tracer_coordinates = [
        # on-fault tracers (35)
        [20.0, 0.0,14.0],
        [20.0, 1.0,14.0],
        [20.0, 2.4,14.0],
        [20.0, 5.0,14.0],
        [20.0, 7.5,14.0],
        [20.0,10.0,14.0],
        [20.0,12.0,14.0],

        [20.0, 0.0, 8.0],
        [20.0, 1.0, 8.0],
        [20.0, 2.4, 8.0],
        [20.0, 5.0, 8.0],
        [20.0, 7.5, 8.0],
        [20.0,10.0, 8.0],
        [20.0,12.0, 8.0],

        [20.0, 0.0,20.0],
        [20.0, 1.0,20.0],
        [20.0, 2.4,20.0],
        [20.0, 5.0,20.0],
        [20.0, 7.5,20.0],
        [20.0,10.0,20.0],
        [20.0,12.0,20.0],

        [20.0, 0.0,26.0],
        [20.0, 1.0,26.0],
        [20.0, 2.4,26.0],
        [20.0, 5.0,26.0],
        [20.0, 7.5,26.0],
        [20.0,10.0,26.0],
        [20.0,12.0,26.0],

        [20.0, 0.0,32.0],
        [20.0, 1.0,32.0],
        [20.0, 2.4,32.0],
        [20.0, 5.0,32.0],
        [20.0, 7.5,32.0],
        [20.0,10.0,32.0],
        [20.0,12.0,32.0],

        # off-fault tracers (56)
        [17.0, 0.0,10.0],
        [17.0, 2.4,10.0],
        [17.0, 0.0, 0.0],
        [17.0, 2.4, 0.0],
        [17.0, 0.0,20.0],
        [17.0, 2.4,20.0],
        [17.0, 0.0,30.0],
        [17.0, 2.4,30.0],
        [17.0, 0.0,40.0],
        [17.0, 2.4,40.0],

        [11.0, 0.0,10.0],
        [11.0, 2.4,10.0],
        [11.0, 0.0, 0.0],
        [11.0, 2.4, 0.0],
        [11.0, 0.0,20.0],
        [11.0, 2.4,20.0],
        [11.0, 0.0,30.0],
        [11.0, 2.4,30.0],
        [11.0, 0.0,40.0],
        [11.0, 2.4,40.0],

        [ 5.0, 0.0, 5.0],
        [ 5.0, 2.4, 5.0],
        [ 5.0, 0.0,20.0],
        [ 5.0, 2.4,20.0],
        [ 5.0, 0.0,35.0],
        [ 5.0, 2.4,35.0],

        [20.0, 0.0, 0.0],
        [20.0, 2.4, 0.0],
        [20.0, 0.0,40.0],
        [20.0, 2.4,40.0],

        [23.0, 0.0,10.0],
        [23.0, 2.4,10.0],
        [23.0, 0.0, 0.0],
        [23.0, 2.4, 0.0],
        [23.0, 0.0,20.0],
        [23.0, 2.4,20.0],
        [23.0, 0.0,30.0],
        [23.0, 2.4,30.0],
        [23.0, 0.0,40.0],
        [23.0, 2.4,40.0],

        [29.0, 0.0,10.0],
        [29.0, 2.4,10.0],
        [29.0, 0.0, 0.0],
        [29.0, 2.4, 0.0],
        [29.0, 0.0,20.0],
        [29.0, 2.4,20.0],
        [29.0, 0.0,30.0],
        [29.0, 2.4,30.0],
        [29.0, 0.0,40.0],
        [29.0, 2.4,40.0],

        [35.0, 0.0, 5.0],
        [35.0, 2.4, 5.0],
        [35.0, 0.0,20.0],
        [35.0, 2.4,20.0],
        [35.0, 0.0,35.0],
        [35.0, 2.4,35.0]
    ]

    def initial_conditions(self):
        return """
      //Q[Shortcuts::rho] = 2.67; # has to be gotten by algorithm
      //Q[Shortcuts::cp ] = 6.0; # has to be gotten by algorithm
      //Q[Shortcuts::cs ] = 3.464; # has to be gotten by algorithm
"""
    
    def generate_required_files(self, order):

        dictionary = { "MODE": order+1 }
        self.generate_file_from_template(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV34/tpv34.yaml.template",
            "tpv34.yaml", dictionary
        )

        self.copy_file_to_current_folder(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV5/tpv34_fault.yaml", 
            "tpv34_fault.yaml"
        )

        return