####################################################################################################
#                                             main.py                                              #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/02/25                                                                                #
#                                                                                                  #
# Purpose: Main script for the BasisREMY tool.                                                     #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
from core.basisremy import BasisREMY


if __name__ == "__main__":
    # BasisREMY().run('./example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR', './output/')
    BasisREMY().run_gui()
    
