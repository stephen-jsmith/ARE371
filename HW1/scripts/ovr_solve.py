import numpy as np
from load_xlsx import load_excel_file
import pandas as pd

FILENAME = None


def eq_solver(df:pd.DataFrame):
    # Placeholder for your function logic
    for i in 10:
        print(i)

    # Define


def unknown_counter(Q_list: list, vars_dict: dict, defaults: dict) -> bool:
    # Placeholder for your function logic
    for Q in Q_list:
        counter = 0
        for key in vars_dict.keys():
            if vars_dict[key] != float or vars_dict[key] != int:
                counter += 1
        if Q == None:
            if counter == 0:
                pass # All variables are known, solve for Q
            else:
                break # At least one unknown, can't solve at this time. Come back later
        elif type(Q) == float or type(Q) == int:
            if counter == 0:
                break # All variables are known, no need to solve
            elif counter == 1:
                pass # Solve for the unknown variable
            else:
                break # More than one unknown, can't solve at this time. Come back later






if __name__ == "__main__":
    eq_solver(load_excel_file(FILENAME))
