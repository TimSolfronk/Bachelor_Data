# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from .Scenario import Scenario

import os

class TPV6(Scenario):
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
          # Near Side of Fault
          # Surface
          [19.8, 0.0,  8.0],
          [19.8, 0.0,  0.0],
          [19.8, 0.0, 32.0],
          # Underground
          [19.8, 7.5,  8.0],
          [19.8, 7.5, 32.0],
          
          # Far Side of Fault
          # Surface
          [20.2, 0.0,  8.0],
          [20.2, 0.0,  0.0],
          [20.2, 0.0, 32.0],
          # Underground
          [20.2, 7.5,  8.0],
          [20.2, 7.5, 32.0],
    ]

    def initial_conditions(self):
        return """
        if(Q[Shortcuts::curve_grid + 0] >= 20.0)   //Look on Curve_grid if x coord is >=20
        {                                           //Far Side
            Q[Shortcuts::rho] = 2.225; 
            Q[Shortcuts::cp ] = 3.75; 
            Q[Shortcuts::cs ] = 2.165; 
        } else {                                    //Near Side
            Q[Shortcuts::rho] = 2.67; 
            Q[Shortcuts::cp ] = 6.0; 
            Q[Shortcuts::cs ] = 3.464; 
        }
"""
    
    def generate_required_files(self, order):

        dictionary = { "MODE": order+1 }
        self.generate_file_from_template(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV6/tpv6.yaml.template",
            "tpv6.yaml", dictionary
        )

        self.copy_file_to_current_folder(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV6/tpv6_fault.yaml", 
            "tpv6_fault.yaml"
        )

        return