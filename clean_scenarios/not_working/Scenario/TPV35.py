# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
from .Scenario import Scenario

import os

class TPV35(Scenario):
    """
    Part of a series of benchmarks by the Statewide California Earthquake Center (SCEC)

    Parkfield 2004 M6 Earthquake. Spontaneous rupture on a vertical strike-slip fault. The first SCEC validation benchmark. 
    Also has real-world data to compare to. 
    TPV35 is based on a spontaneous rupture model from Ma, Custódio, Archuleta, and Liu (2008), 
    Dynamic modeling of the 2004 Mw 6.0 Parkfield, California, earthquake, J. Geophys. Res., 113, B02301, doi:10.1029/2007JB005216.

    The description of the scenario can be found at: https://strike.scec.org/cvws/tpv35docs.html
    """

    domain_offset   = [-1.6,  0., -1.6] #centers the fault which is at 20.
    domain_size     = [43.2, 43.2, 43.2] #leads to a cell size of 0.333 for 81 cells

    end_time = 10.0

    tracer_sets = {
        "l3_o5": [
            #on-fault at far side
            [19.87,8.1, 5.0],
            [19.87,8.1,10.0], 
            [19.87,8.1,15.0], 
            [19.87,8.1,20.0], 
            [19.87,8.1,25.0], 
            [19.87,8.1,30.0], 
            [19.87,8.1,35.0],

            #on-fault at near side
            [19.84,8.1, 5.0],
            [19.84,8.1,10.0],
            [19.84,8.1,15.0],
            [19.84,8.1,20.0],
            [19.84,8.1,25.0],
            [19.84,8.1,30.0],
            [19.84,8.1,35.0],
        ]
    }

    tracer_coordinates = tracer_sets["l3_o5"]

    def initial_conditions(self):
        return """
     if(x[0] >= 20.0)   // if x coord is >=20
        {                                           //Far Side
            if(x[1] <= 1.0)
            {
                Q[Shortcuts::rho] = 2.0;
                Q[Shortcuts::cp ] = 2.0;
                Q[Shortcuts::cs ] = 1.1;
            } else if(x[1] <= 1.8)
            {
                Q[Shortcuts::rho] = 2.3;
                Q[Shortcuts::cp ] = 3.5;
                Q[Shortcuts::cs ] = 2.2;
            } else if(x[1] <= 2.1)
            {
                Q[Shortcuts::rho] = 2.3;
                Q[Shortcuts::cp ] = 4.2;
                Q[Shortcuts::cs ] = 2.8;
            } else if(x[1] <= 3.4)
            {
                Q[Shortcuts::rho] = 2.3;
                Q[Shortcuts::cp ] = 4.8;
                Q[Shortcuts::cs ] = 2.7;
            } else if(x[1] <= 3.9)
            {
                Q[Shortcuts::rho] = 2.3;
                Q[Shortcuts::cp ] = 5.2;
                Q[Shortcuts::cs ] = 2.8;
            } else if(x[1] <= 8.3)
            {
                Q[Shortcuts::rho] = 2.7;
                Q[Shortcuts::cp ] = 5.3;
                Q[Shortcuts::cs ] = 3.2;
            } else if(x[1] <= 12.7)
            {
                Q[Shortcuts::rho] = 2.8;
                Q[Shortcuts::cp ] = 5.7;
                Q[Shortcuts::cs ] = 3.7;
            } else if(x[1] <= 17.5)
            {
                Q[Shortcuts::rho] = 2.8;
                Q[Shortcuts::cp ] = 6.5;
                Q[Shortcuts::cs ] = 3.8;
            } else if(x[1] <= 20.3)
            {
                Q[Shortcuts::rho] = 2.8;
                Q[Shortcuts::cp ] = 6.7;
                Q[Shortcuts::cs ] = 4.3;
            } else
            {
                Q[Shortcuts::rho] = 2.8;
                Q[Shortcuts::cp ] = 7.3;
                Q[Shortcuts::cs ] = 4.3;
            }
        } else {                                    //Near Side
            if(x[1] <= 1.0)
            {
                Q[Shortcuts::rho] = 2.0;
                Q[Shortcuts::cp ] = 2.0;
                Q[Shortcuts::cs ] = 1.1;
            } else if(x[1] <= 2.0)
            {
                Q[Shortcuts::rho] = 2.3;
                Q[Shortcuts::cp ] = 3.5;
                Q[Shortcuts::cs ] = 2.0;
            } else if(x[1] <= 3.0)
            {
                Q[Shortcuts::rho] = 2.3;
                Q[Shortcuts::cp ] = 4.5;
                Q[Shortcuts::cs ] = 2.5;
            } else if(x[1] <= 3.5)
            {
                Q[Shortcuts::rho] = 2.5;
                Q[Shortcuts::cp ] = 5.2;
                Q[Shortcuts::cs ] = 3.0;
            } else if(x[1] <= 5.8)
            {
                Q[Shortcuts::rho] = 2.7;
                Q[Shortcuts::cp ] = 5.7;
                Q[Shortcuts::cs ] = 3.2;
            } else if(x[1] <= 14.1)
            {
                Q[Shortcuts::rho] = 2.7;
                Q[Shortcuts::cp ] = 6.2;
                Q[Shortcuts::cs ] = 3.6;
            } else if(x[1] <= 17.1)
            {
                Q[Shortcuts::rho] = 2.8;
                Q[Shortcuts::cp ] = 6.8;
                Q[Shortcuts::cs ] = 3.6;
            } else if(x[1] <= 20.4)
            {
                Q[Shortcuts::rho] = 2.8;
                Q[Shortcuts::cp ] = 6.8;
                Q[Shortcuts::cs ] = 4.3;
            } else
            {
                Q[Shortcuts::rho] = 2.8;
                Q[Shortcuts::cp ] = 7.3;
                Q[Shortcuts::cs ] = 4.3;
            }
        }

"""
    
    def generate_required_files(self, order):

        dictionary = { "MODE": order+1 }
        self.generate_file_from_template(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV35/tpv35.yaml.template",
            "tpv35.yaml", dictionary
        )

        self.copy_file_to_current_folder(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV35/tpv35_fault.yaml", 
            "tpv35_fault.yaml"
        )

        self.copy_file_to_current_folder(
            os.path.dirname(os.path.realpath(__file__))+"/specs/TPV35/tpv35_muS_sXZ.nc", 
            "tpv35_muS_sXZ.nc"
        )

        return
