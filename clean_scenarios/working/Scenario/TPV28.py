# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from .Scenario import Scenario

import os

class TPV28(Scenario):
    """
    Part of a series of benchmarks by the Statewide California Earthquake Center (SCEC)

    Single non-planar fault with two hills.
    On the plane, the initial stresses are defined by the surface.
    There are homogeneous initial stress conditions outside of the plane.

    The description of the scenario can be found at: https://strike.scec.org/cvws/tpv28docs.html
    """

    domain_offset   = [0.56,  0., 0.56]
    domain_size     = [38.88, 38.88, 38.88] #leads to a cell size of 0.3333 for 81 cells

    end_time = 13.0
    tracer_sets = {
        "l3_o7": [
            # off-fault at surface
            [16.93, 0.0,  5.0],
            [23.08, 0.0,  5.0],
            [16.93, 0.0, 20.0],
            [23.08, 0.0, 20.0],
            [16.93, 0.0, 35.0],
            [23.08, 0.0, 35.0],

            # on-fault near side
            # surface
            [19.895, 0.0,  8.0],
            [19.895, 0.0, 20.0], #surface above fault hypocenter
            [19.895, 0.0, 32.0],
            #on-fault, general (21)
            [19.895,  3.0, 20.0], #above fault hypocenter
            [19.598,  6.0,  9.5], #above left patch center    #should be at 19.7
            [19.598,  6.0, 30.5], #above right patch center   #should be at 19.7

            [19.895,  7.5,  5.0],
            [19.598,  7.5,  8.0],                             #should be at 19.7
            [19.302,  7.5,  9.5], #center of left patch       #should be at 19.4
            [19.598,  7.5, 11.0],                             #should be at 19.7
            [19.895,  7.5, 14.0],
            [19.895,  7.5, 20.0], #on fault hypocenter
            [19.895,  7.5, 26.0],
            [19.598,  7.5, 29.0],                             #should be at 19.7
            [19.302,  7.5, 30.5], #center of right patch      #should be at 19.4
            [19.598,  7.5, 32.0],                             #should be at 19.7
            [19.895,  7.5, 35.0],
            [19.598,  9.0,  9.5], #below left patch center    #should be at 19.7
            [19.598,  9.0, 30.5],  #below right patch center  #should be at 19.7
            [19.895, 12.0, 20.0], #below fault hypocenter

            # on-fault far side
            # surface
            [19.915, 0.0,  8.0],
            [19.915, 0.0, 20.0], #surface above fault hypocenter
            [19.915, 0.0, 32.0],
            #on-fault, general (21)
            [19.915,  3.0, 20.0], #above fault hypocenter
            [19.617,  6.0,  9.5], #above left patch center    #should be at 19.7
            [19.617,  6.0, 30.5], #above right patch center   #should be at 19.7

            [19.915,  7.5,  5.0],
            [19.617,  7.5,  8.0],                             #should be at 19.7
            [19.321,  7.5,  9.5], #center of left patch       #should be at 19.4
            [19.617,  7.5, 11.0],                             #should be at 19.7
            [19.915,  7.5, 14.0],
            [19.915,  7.5, 20.0], #on fault hypocenter
            [19.915,  7.5, 26.0],
            [19.617,  7.5, 29.0],                             #should be at 19.7
            [19.321,  7.5, 30.5], #center of right patch      #should be at 19.4
            [19.617,  7.5, 32.0],                             #should be at 19.7
            [19.915,  7.5, 35.0],
            [19.617,  9.0,  9.5], #below left patch center    #should be at 19.7
            [19.617,  9.0, 30.5],  #below right patch center  #should be at 19.7
            [19.915, 12.0, 20.0] #below fault hypocenter
        ],
    }

    tracer_coordinates = tracer_sets["l3_o7"]

    def initial_conditions(self):
        return """
      Q[Shortcuts::rho] = 2.67;
      Q[Shortcuts::cp ] = 6.0;
      Q[Shortcuts::cs ] = 3.464;
"""
    
    def generate_required_files(self, order):

        dictionary = { "MODE": order+1 }
        self.generate_file_from_template(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV28/tpv28.yaml.template",
            "tpv28.yaml", dictionary
        )

        self.copy_file_to_current_folder(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV28/tpv28_fault.yaml",
            "tpv28_fault.yaml"
        )

        return
