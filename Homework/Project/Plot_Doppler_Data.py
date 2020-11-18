""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import Statements
# ============================================================================
import matplotlib.pyplot as plt

from Project_Utilities import Nuclide_Data
            
# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
if __name__ == "__main__":
    H1 = Nuclide_Data('H1', 1.008, [1, 1, 0])
    # O16 = Nuclide_Data('O16', 15.995, [1, 1, 0])
    # U235 = Nuclide_Data('U235', 235.044, [1, 1, 1])
    # U238 = Nuclide_Data('U238', 238.051, [1, 1, 1])
    
    H1.Load_Doppler_Data([600, 900, 1200])