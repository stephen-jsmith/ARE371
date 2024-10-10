# =============================================================================
# Copyright (c) 2024 Stephen Smith
# All rights reserved.
# =============================================================================

# Import necessary libraries
import json
import pandas as pd
from sympy import Eq, symbols as sym, solve

# Import scripts for data loading and manipulation
from scripts.load_xlsx import *
from scripts.xlsx_to_json import *
from scripts.ovr_solve import *

# Import solver scripts for each system
from scripts.solvers.in_win import *
from scripts.solvers.out_win import *
from scripts.solvers.qhvac import *
from scripts.solvers.swall_in import *
from scripts.solvers.swall_out import *


def main():
    # Your code here
    pass

if __name__ == "__main__":
    main()
