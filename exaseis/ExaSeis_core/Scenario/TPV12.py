# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from .Scenario import Scenario

import os
import math

class TPV12(Scenario):
    """
    Part of a series of benchmarks by the Statewide California Earthquake Center (SCEC)

    Spontaneous rupture on a vertical strike-slip fault in a homogeneous halfspace.
    There are slightly heterogeneous initial stress conditions. 

    The description of the scenario can be found at: https://strike.scec.org/cvws/tpv5docs.html
    """

    domain_offset   = [6.5,  0., 6.5] #centers the fault which is at 20.
    domain_size     = [27., 27., 27.] #leads to a cell size of 0.333 for 81 cells

    end_time = 8.0


    angle_to_y_axis = 30 * (3.14159265359/180) # 30° in radians
    body_mult = math.tan(angle_to_y_axis)
    off_fault_tracer_coordinates = [
        [17.0, 0.0,20.0],
        [18.0, 0.0,20.0],
        [19.0, 0.0,20.0],
        [21.0, 0.0,20.0],
        [22.0, 0.0,20.0],
        [23.0, 0.0,20.0],
        [19.0, 0.3,20.0],
        [19.5, 0.3,20.0],
        [20.5, 0.3,20.0],
        [21.0, 0.3,20.0],
        [17.0, 0.0,32.0],
        [23.0, 0.0,32.0]
    ]
    for tracer in off_fault_tracer_coordinates:
        tracer[0] = tracer[0] + tracer[1] * body_mult

    body_mult = math.sin(angle_to_y_axis)
    dip_mult = math.cos(angle_to_y_axis)
    on_fault_tracer_coordinates = [
          [20.0, 0.0,20.0],
          [20.0, 0.0,24.5],
          [20.0, 0.0,32.0],
          [20.0, 1.5,20.0],
          [20.0, 3.0,20.0],
          [20.0, 4.5,20.0],
          [20.0, 7.5,20.0],
          [20.0, 7.5,24.5],
          [20.0, 7.5,32.0],
          [20.0,12.0,20.0]
    ]
    for tracer in on_fault_tracer_coordinates:
        tracer[0] = tracer[0] + tracer[1] * body_mult
        tracer[1] = tracer[1] * dip_mult


    tracer_coordinates = off_fault_tracer_coordinates + on_fault_tracer_coordinates

    
    def initial_conditions(self):
        return """
      Q[Shortcuts::rho] = 2.7;
      Q[Shortcuts::cp ] = 5.716;
      Q[Shortcuts::cs ] = 3.3;
"""
    
    def generate_required_files(self, order):

        dictionary = { "MODE": order+1 }
        self.generate_file_from_template(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV12/tpv12.yaml.template",
            "tpv12.yaml", dictionary
        )

        self.copy_file_to_current_folder(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV12/tpv12_fault.yaml", 
            "tpv12_fault.yaml"
        )

        return