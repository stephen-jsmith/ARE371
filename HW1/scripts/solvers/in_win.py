import json
import pandas as pd
from sympy import Eq, symbols as sym, solve
import os
import os.path
import shutil

# Set path for getting data in higher/different directories
os.chdir('..')
os.chdir('..')
d = os.getcwd()
for item in os.listdir(d):
    print(item)

from scripts.load_xlsx import *


def solve_q_sol_in_win(alpha_glass_in=None, A_w_in=None, q_dir_in=None, q_diff_in=None, q_refl_in=None, tau_in=None, A_s_in=None, A_i_in=None, q_sol_in=None):
    """
    Solves the equation for solar heat gain through a window.
    
    The function solves the following equation for solar heat gain through a window:
    
    q_{\text{sol\_w}} \cdot \tau \cdot \left( \frac{A_s}{A_s + A_w + A_i} \right) - q_{\text{sol}} = 0

    Parameters:
    alpha_glass_in (float or None): The solar absorptance of the glass (alpha_glass). If None, it will be treated as a symbolic variable.
    A_w_in (float or None): The area of the window (A_w). If None, it will be treated as a symbolic variable.
    q_dir_in (float or None): The direct solar radiation (q_dir). If None, it will be treated as a symbolic variable.
    q_diff_in (float or None): The diffuse solar radiation (q_diff). If None, it will be treated as a symbolic variable.
    q_refl_in (float or None): The reflected solar radiation (q_refl). If None, it will be treated as a symbolic variable.
    tau_in (float or None): The transmittance of the window (tau). If None, it will be treated as a symbolic variable.
    A_s_in (float or None): The area of the sunlit surface (A_s). If None, it will be treated as a symbolic variable.
    A_i_in (float or None): The area of the internal surface (A_i). If None, it will be treated as a symbolic variable.
    q_sol_in (float or None): The total solar heat gain (q_sol). If None, it will be treated as a symbolic variable.

    Returns:
    list: A list of solutions for the equation.
    """
    alpha_glass, A_w, q_dir, q_diff, q_refl, tau, A_s, A_i, q_sol = sym('alpha_glass A_w q_dir q_diff q_refl tau A_s A_i q_sol')
    f = Eq(0.5 * alpha_glass * A_w * (q_dir + q_diff + q_refl) + tau * A_w * (q_dir + q_diff + q_refl) * (A_s/(A_s + A_w + A_i)) - q_sol, 0)
    solved_var = None
    variables = {
        'alpha_glass': alpha_glass_in, 'A_w': A_w_in, 'q_dir': q_dir_in, 
        'q_diff': q_diff_in, 'q_refl': q_refl_in, 'tau': tau_in, 
        'A_s': A_s_in, 'A_i': A_i_in, 'q_sol': q_sol_in
    }
    for var, value in variables.items():
        if value is not None:
            f = f.subs(var, value)
        else:
            solved_var = var  # Assumes only one variable is unknown
    return [solved_var, solve(f)]

def solve_q_cond_wo_wi(A_w_in, k_in, L_in, T_wo_in, T_wi_in, q_cond_wo_wi_in):
    """
    Solves the equation for the conductive heat transfer through a wall without insulation.

    The function solves the following equation for the conductive heat transfer through a wall without insulation:

    .. math::
        q_{\text{cond\_wo\_wi}} = \frac{k \cdot A_w \cdot (T_{\text{wo}} - T_{\text{wi}})}{L}

    Parameters:
    A_w_in (float or None): The area of the wall (A_w). If None, it will be treated as a symbolic variable.
    k_in (float or None): The thermal conductivity of the wall (k). If None, it will be treated as a symbolic variable.
    L_in (float or None): The thickness of the wall (L). If None, it will be treated as a symbolic variable.
    T_wo_in (float or None): The temperature on the outside of the wall (T_wo). If None, it will be treated as a symbolic variable.
    T_wi_in (float or None): The temperature on the inside of the wall (T_wi). If None, it will be treated as a symbolic variable.
    q_cond_wo_wi_in (float or None): The conductive heat transfer through the wall without insulation (q_cond_wo_wi). If None, it will be treated as a symbolic variable.

    Returns:
    list: A list of solutions for the equation.
    """
    A_w, k, L, T_wo, T_wi, q_cond_wo_wi = sym('A_w k L T_wo T_wi q_cond_wo_wi')
    f = Eq(k * A_w * (T_wo - T_wi) / L - q_cond_wo_wi, 0)
    solved_var = None
    variables = {'A_w': A_w_in, 'k': k_in, 'L': L_in, 'T_wo': T_wo_in, 'T_wi': T_wi_in, 'q_cond_wo_wi': q_cond_wo_wi_in}
    for var, value in variables.items():
        if value is not None:
            f = f.subs(var, value)
        else:
            solved_var = var  # Assumes only one variable is unknown
    return [solved_var, solve(f)]

def solve_q_conv_ia_wi(A_w_in, h_in, T_ia_in, T_wi_in, q_conv_ia_wi_in):
    """
    Solves the equation for convective heat transfer between indoor air and a wall with insulation.

    .. math::
        q_{\text{conv\_ia\_wi}} = A_w \cdot \left(\frac{1.8}{h}\right) \cdot \left(\left|T_{\text{ia}} - T_{\text{wi}}\right|\right)^{0.29} \cdot \left(T_{\text{ia}} - T_{\text{wi}}\right)

    Parameters:
    A_w_in (float or None): Area of the wall with insulation (m^2). If None, it will be solved for.
    h_in (float or None): Convective heat transfer coefficient (W/m^2K). If None, it will be solved for.
    T_ia_in (float or None): Indoor air temperature (K). If None, it will be solved for.
    T_wi_in (float or None): Wall temperature with insulation (K). If None, it will be solved for.
    q_conv_ia_wi_in (float or None): Convective heat transfer rate (W). If None, it will be solved for.
    Returns:
    list: A list containing the name of the solved variable as a string and the solution as a sympy expression.
    Example:
    >>> solve_q_conv_ia_wi(10, 5, 300, 290, None)
    ['q_conv_ia_wi', 10 * (1.8/5) * ((abs(300 - 290))**0.29) * (300 - 290)]
    """

    A_w, h, T_ia, T_wi, q_conv_ia_wi = sym('A_w h T_ia T_wi q_conv_ia_wi')
    f = Eq(A_w * (1.8/h) * ((abs(T_ia - T_wi))**0.29) * (T_ia - T_wi) - q_conv_ia_wi, 0)
    solved_var = None
    variables = {'A_w': A_w_in, 'h': h_in, 'T_ia': T_ia_in, 'T_wi': T_wi_in, 'q_conv_ia_wi': q_conv_ia_wi_in}
    for var, value in variables.items():
        if value is not None:
            f = f.subs(var, value)
        else:
            solved_var = var # Assumes only one variable is unknown
    return [solved_var, solve(f)]

def solve_unknowns_in_win(df:pd.Series, defaults:dict):
    """
    Placeholder function for finding unknown variables in the system.
    """
    # Load Defaults
    vars_involved = {
        'q_sol' : None,
        'q_cond_wo_wi' : None,
        'q_conv_ia_wi' : None,
        'q_lwr_i_wi' : None,
        'alpha_glass' : defaults["Window"]["alpha_swg"],
        'A_w' : defaults["Window"]["A_w"], # Area of the window
        'q_dir' : df.get('q_dir'),
        'q_diff' : df.get('q_diff'),
        'q_refl' : df.get('q_refl'),
        'tau' : defaults["Window"]["tau"], # Transmittance of the window. Solved from sum of parts must equal 1. 
        'A_s' : defaults["SW"]["A_s"], # Area of south wall
        'A_i' : defaults["Internal"]["A_i"], # Area of internal surfaces
        'k' : defaults["Window"]["k"], # Thermal conductivity of the window
        'L' : defaults["Window"]["L"], # Length of the window
        'T_wo' : None, # Temperature of the wall, outside. Expected unknown
        'T_wi' : None, # Temperature of the wall, inside. Expected unknown
        'h' : defaults["Window"]["H_w"], # Height of the element being considered. This case, window
        'T_ia' : defaults["T_ia"], # Temperature of internal air
        'F_wi_i' : None, # View factor between inside wall and internal surfaces
        'E_wi' : defaults["window"]["E_wi"], # Emissivity of window to internal surfaces
        'E_i' : defaults["Window"]["E_i"], # Emissivity of internal surfaces
        'sigma' : None,
        'T_i' : defaults["T_ia"],
    }

    # Sum the unknowns
    init_unknowns = 0
    for i in vars_involved.keys():
        if vars_involved[i] == None:
            init_unknowns += 1
    print(f"Number of initial unknowns: {init_unknowns}")
    if init_unknowns == 0:
        print("All variables are known. No need to solve.")
        return [df, init_unknowns, 0]

    # Take a pass at solving the system
    q_sol_in_win_vars = [vars_involved['alpha_glass'], vars_involved['A_w'], vars_involved['q_dir'], vars_involved['q_diff'], vars_involved['q_refl'], vars_involved['tau'], vars_involved['A_s'], vars_involved['A_i'], vars_involved['q_sol']]
    counter = 0
    for i in q_sol_in_win_vars:
        if i == None:
            counter += 1
    if counter == 1:
        print("Solving for q_sol_in...")
        res = solve_q_sol_in_win(vars_involved['alpha_glass'], vars_involved['A_w'], vars_involved['q_dir'], vars_involved['q_diff'], vars_involved['q_refl'], vars_involved['tau'], vars_involved['A_s'], vars_involved['A_i'], vars_involved['q_sol'])
        print(f"Solved q_sol_in EQ. Results: \n Var Solved for: {res[0]}, Solution: {res[1]}. \n Updating dataframe...")
        vars_involved[res[0]] = res[1]
    else:
        print("Cannot solve for q_sol_in_win at this time. Progressing to next equation.")

    q_cond_wo_wi_vars = [vars_involved['A_w'], vars_involved['k'], vars_involved['L'], vars_involved['T_wo'], vars_involved['T_wi'], vars_involved['q_cond_wo_wi']]
    counter = 0
    for i in q_cond_wo_wi_vars:
        if i == None:
            counter += 1
    if counter == 1:
        print("Solving for q_cond_wo_wi...")
        res = solve_q_cond_wo_wi(vars_involved['A_w'], vars_involved['k'], vars_involved['L'], vars_involved['T_wo'], vars_involved['T_wi'], vars_involved['q_cond_wo_wi'])
        print(f"Solved q_cond_wo_wi EQ. Results: \n Var Solved for: {res[0]}, Solution: {res[1]}. \n Updating dataframe...")
        vars_involved[res[0]] = res[1]
    else:
        print("Cannot solve for q_cond_wo_wi at this time. Progressing to next equation.")

    q_conv_ia_wi_vars = [vars_involved['A_w'], vars_involved['h'], vars_involved['T_ia'], vars_involved['T_wi'], vars_involved['q_conv_ia_wi']]
    counter = 0
    for i in q_conv_ia_wi_vars:
        if i == None:
            counter += 1
    if counter == 1:
        print("Solving for q_conv_ia_wi...")
        res = solve_q_conv_ia_wi(vars_involved['A_w'], vars_involved['h'], vars_involved['T_ia'], vars_involved['T_wi'], vars_involved['q_conv_ia_wi'])
        print(f"Solved q_conv_ia_wi EQ. Results: \n Var Solved for: {res[0]}, Solution: {res[1]}. \n Updating dataframe...")
        vars_involved[res[0]] = res[1]
    else:
        print("Cannot solve for q_conv_ia_wi at this time. Progressing to next equation.")

    fin_unknowns = 0
    for i in vars_involved.keys():
        if vars_involved[i] == None:
            fin_unknowns += 1
    print(f"Number of final unknowns: {init_unknowns}")
    # Return the updated dataframe, which should now have some solved values. Also return the number of unknowns started with, and number of unknowns now.

    # UPDATE DATAFRAME HERE
    # df['alpha_glass'] = vars_involved['alpha_glass']  <--- Something like this, but figure it out first
    return [df, init_unknowns, fin_unknowns]  # <-- Allows for the logic of if the fin_unknowns < init_unknowns, run the function again, as some values might have been solved for. If not, move on to the next system.