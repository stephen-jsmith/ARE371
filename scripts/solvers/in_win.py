import json
import pandas as pd
from sympy import Eq, symbols as sym, solve

def solve_q_sol_in_win(alpha_glass_in, A_w_in, q_dir_in, q_diff_in, q_refl_in, tau_in, A_s_in, A_i_in, q_sol_in):
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
    f = Eq(0.5 * alpha_glass * A_w * (q_dir + q_diff + q_refl) + tau * A_w * (q_dir + q_diff + q_ref) * (A_s/(A_s + A_w + A_i)) - q_sol, 0)
    if alpha_glass_in != None:
        f.subs(alpha_glass, alpha_glass_in)
    if A_w_in != None:
        f.subs(A_w, A_w_in)
    if q_dir_in != None:
        f.subs(q_dir, q_dir_in)
    if q_diff_in != None:
        f.subs(q_diff, q_diff_in)
    if q_refl_in != None:
        f.subs(q_refl, q_refl_in)
    if tau_in != None:
        f.subs(tau, tau_in)
    if A_s_in != None:
        f.subs(A_s, A_s_in)
    if A_i_in != None:
        f.subs(A_i, A_i_in)
    if q_sol_in != None:
        f.subs(q_sol, q_sol_in)
    return solve(f)

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
    if A_w_in != None:
        f.subs(A_w, A_w_in)
    if k_in != None:
        f.subs(k, k_in)
    if L_in != None:
        f.subs(L, L_in)
    if T_wo_in != None:
        f.subs(T_wo, T_wo_in)
    if T_wi_in != None:
        f.subs(T_wi, T_wi_in)
    if q_cond_wo_wi_in != None:
        f.subs(q_cond_wo_wi, q_cond_wo_wi_in)
    return solve(f)

def solve_q_conv_ia_wi(A_w_in, h_in, T_ia_in, T_wi_in, q_conv_ia_wi_in):
    """
    Solves the equation for convective heat transfer between indoor air and a wall with insulation.
    """

    A_w, h, T_ia, T_wi, q_conv_ia_wi = sym('A_w h T_ia T_wi q_conv_ia_wi')
    f = Eq(A_w * (1.8/h) * ((abs(T_ia - T_wi))**0.29) * (T_ia - T_wi) - q_conv_ia_wi, 0)
    if A_w_in != None:
        f.subs(A_w, A_w_in)
    if h_in != None:
        f.subs(h, h_in)
    if T_ia_in != None:
        f.subs(T_ia, T_ia_in)
    if T_wi_in != None:
        f.subs(T_wi, T_wi_in)
    if q_conv_ia_wi_in != None:
        f.subs(q_conv_ia_wi, q_conv_ia_wi_in)
    return solve(f)
