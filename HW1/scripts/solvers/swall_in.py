import json
import pandas as pd
from sympy import Eq, symbols as sym, solve

def solve_q_sol_swall_in(q_sol_w_in, tau_in, A_s_in, A_w_in, A_i_in, q_sol_in):
    """
    Solves the equation for solar heat gain through a wall.

    The function solves the following equation for solar heat gain through a wall:
    
    q_{\text{sol\_w}} \cdot \tau \cdot \left( \frac{A_s}{A_s + A_w + A_i} \right) - q_{\text{sol}} = 0

    Parameters:
    q_sol_w_in (float or None): The solar heat gain through the wall (q_sol_w). If None, it will be treated as a symbolic variable.
    tau_in (float or None): The transmittance of the wall (tau). If None, it will be treated as a symbolic variable.
    A_s_in (float or None): The area of the sunlit surface (A_s). If None, it will be treated as a symbolic variable.
    A_w_in (float or None): The area of the wall (A_w). If None, it will be treated as a symbolic variable.
    A_i_in (float or None): The area of the internal surface (A_i). If None, it will be treated as a symbolic variable.
    q_sol_in (float or None): The total solar heat gain (q_sol). If None, it will be treated as a symbolic variable.

    Returns:
    list: A list of solutions for the equation.
    """
    q_sol_w, tau, A_s, A_w, A_i, q_sol = sym('q_sol_w tau A_s A_w A_i q_sol')
    f = Eq(q_sol_w * tau * (A_s/(A_s + A_w + A_i)) - q_sol, 0)
    solved_var = None
    if q_sol_w_in != None:
        f.subs(q_sol_w, q_sol_w_in)
    else:
        solved_var = q_sol_w
    if tau_in != None:
        f.subs(tau, tau_in)
    if A_s_in != None:
        f.subs(A_s, A_s_in)
    if A_w_in != None:
        f.subs(A_w, A_w_in)
    if A_i_in != None:
        f.subs(A_i, A_i_in)
    if q_sol_in != None:
        f.subs(q_sol, q_sol_in)
    return solve(f)

def solve_q_sol_w(A_w_in, q_dir_in, q_diff_in, q_refl_in, q_sol_w):
    """
    Solves the equation for the solar heat gain through a wall.

    The function solves the following equation for the given inputs:

    .. math::
        Q_{sol} + Q_{conv,oa,so} + Q_{lwr,gr,so} + Q_{lwr,sky,so} + Q_{cond,si,so} = 0

    Parameters:
    A_w_in (float or None): Area of the wall. If None, it will be treated as a symbolic variable.
    q_dir_in (float or None): Direct solar radiation. If None, it will be treated as a symbolic variable.
    q_diff_in (float or None): Diffuse solar radiation. If None, it will be treated as a symbolic variable.
    q_refl_in (float or None): Reflected solar radiation. If None, it will be treated as a symbolic variable.
    q_sol_w (float or None): Solar heat gain through the wall. If None, it will be treated as a symbolic variable.
    Returns:
    list: Solutions to the equation.
    """

    A_w, q_dir, q_diff, q_refl, q_sol_w = sym('A_w q_dir q_diff q_refl q_sol_w')
    f = Eq((q_dir + q_diff + q_refl) * A_w - q_sol_w, 0)
    if A_w_in != None:
        f = f.subs(A_w, A_w_in)
    if q_dir_in is not None:
        f = f.subs(q_dir, q_dir_in)
    if q_diff_in != None:
        f = f.subs(q_diff, q_diff_in)
    if q_refl_in != None:
        f = f.subs(q_refl, q_refl_in)
    if q_sol_w != None:
        f = f.subs(q_sol_w, q_sol_w)
    return solve(f)

def solve_q_cond_so_si(A_s_in, k_in, L_in, T_so_in, T_si_in, Q_cond_so_si_in):
    """
    Solves the heat conduction equation for a wall given the inputs. 
    
    The function solves the following equation for the given inputs:

    .. math::
    \frac{k \cdot A_s \cdot (T_{so} - T_{si})}{L} - Q_{cond\_so\_si} = 0

    Parameters:
    A_s_in (float or None): Cross-sectional area of the wall (A_s). If None, it will be treated as a variable.
    k_in (float or None): Thermal conductivity of the wall material (k). If None, it will be treated as a variable.
    L_in (float or None): Thickness of the wall (L). If None, it will be treated as a variable.
    T_so_in (float or None): Temperature on the outside surface of the wall (T_so). If None, it will be treated as a variable.
    T_si_in (float or None): Temperature on the inside surface of the wall (T_si). If None, it will be treated as a variable.
    Q_cond_so_si_in (float or None): Heat transfer rate through the wall (Q_cond_so_si). If None, it will be treated as a variable.

    Returns:
    list: Solutions for the unknown variable(s) in the equation.
    """

    A_s, k, L, T_so, T_si, Q_cond_so_si = sym('A_s k L T_so T_si Q_cond_so_si')
    f = Eq((k * A_s * (T_so - T_si) / L) - Q_cond_so_si, 0)
    if A_s_in is not None:
        f = f.subs(A_s, A_s_in)
    if k_in is not None:
        f = f.subs(k, k_in)
    if L_in is not None:
        f = f.subs(L, L_in)
    if T_so_in is not None:
        f = f.subs(T_so, T_so_in)
    if T_si_in is not None:
        f = f.subs(T_si, T_si_in)
    if Q_cond_so_si_in is not None:
        f = f.subs(Q_cond_so_si, Q_cond_so_si_in)
    return solve(f)

def solve_q_conv_ia_si(A_s_in, h_in, T_ia_in, T_si_in, Q_conv_ia_si_in):
    """
    Solves for the unknown variable in the convective heat transfer equation between indoor air and a surface.

    The function solves the following equation for the given inputs:

    .. math::
    (A_s * (1.8 / h) * (|T_{ia} - T_{si}|^{0.29}) * (T_{ia} - T_{si})) - Q_{conv\_ia\_si} = 0

    Parameters:
    A_s_in (float or None): Surface area (A_s). If None, it will be treated as a variable to solve for.
    h_in (float or None): Convective heat transfer coefficient (h). If None, it will be treated as a variable to solve for.
    T_ia_in (float or None): Indoor air temperature (T_ia). If None, it will be treated as a variable to solve for.
    T_si_in (float or None): Surface temperature (T_si). If None, it will be treated as a variable to solve for.
    Q_conv_ia_si_in (float or None): Convective heat transfer rate (Q_conv_ia_si). If None, it will be treated as a variable to solve for.

    Returns:
    list: Solutions for the unknown variable(s).
    """

    A_s, h, T_ia, T_si, Q_conv_ia_si = sym('A_s h T_ia T_si Q_conv_ia_si')
    f = Eq((A_s*(1.8/h)*(abs(T_ia - T_si)**0.29)*(T_ia-T_si)) - Q_conv_ia_si, 0)
    if A_s_in is not None:
        f = f.subs(A_s, A_s_in)
    if h_in is not None:
        f = f.subs(h, h_in)
    if T_ia_in is not None:
        f = f.subs(T_ia, T_ia_in)
    if T_si_in is not None:
        f = f.subs(T_si, T_si_in)
    if Q_conv_ia_si_in is not None:
        f = f.subs(Q_conv_ia_si, Q_conv_ia_si_in)
    return solve(f)

def solve_q_lwr_i_si(A_s_in, F_si_i_in, E_si_in, E_i_in, sigma_in, T_i_in, T_si_in,q_lwr_i_si_in):
    """
    Solves for the lower radiative heat flux between surfaces i and si.

    The function solves the following equation for the given inputs:

    .. math::
    A_s * F_si_i * E_si * E_i * sigma * (T_i^4 + T_si^4) - q_lwr_i_si = 0

    Parameters:
    A_s_in (float or None): Surface area of surface s.
    F_si_i_in (float or None): View factor from surface si to surface i.
    E_si_in (float or None): Emissivity of surface si.
    E_i_in (float or None): Emissivity of surface i.
    sigma_in (float or None): Stefan-Boltzmann constant.
    T_i_in (float or None): Temperature of surface i in Kelvin.
    T_si_in (float or None): Temperature of surface si in Kelvin.
    q_lwr_i_si_in (float or None): Lower radiative heat flux between surfaces i and si.

    Returns:
    list: Solutions for the lower radiative heat flux (q_lwr_i_si).
    """

    A_s, F_si_i, E_si, E_i, sigma, T_i, T_si, q_lwr_i_si = sym('A_s F_si_i E_si E_i sigma T_i T_si q_lwr_i_si')
    f = Eq(A_s*F_si_i*E_si*E_i*sigma*(T_i**4 + T_si**4)-q_lwr_i_si, 0)
    if A_s_in is not None:
        f = f.subs(A_s, A_s_in)
    if F_si_i_in is not None:
        f = f.subs(F_si_i, F_si_i_in)
    if E_si_in is not None:
        f = f.subs(E_si, E_si_in)
    if E_i_in is not None:
        f = f.subs(E_i, E_i_in)
    if sigma_in is not None:
        f = f.subs(sigma, sigma_in)
    if T_i_in is not None:
        f = f.subs(T_i, T_i_in)
    if T_si_in is not None:
        f = f.subs(T_si, T_si_in)
    return solve(f)

def solve_ovr_swall_out(q_sol_in, q_cond_so_si_in, q_conv_ia_si_in, q_lwr_i_si_in):
    """
    Solves the equation q_sol + q_cond_so_si + q_conv_ia_si + q_lwr_i_si = 0 for the unknown variables.

    The function solves the following equation for the given inputs:
    
    .. math::
    q_{sol} + q_{cond\_so\_si} + q_{conv\_ia\_si} + q_{lwr\_i\_si} = 0

    Parameters:
    q_sol_in (float or None): The value of q_sol. If None, it will be treated as an unknown.
    q_cond_so_si_in (float or None): The value of q_cond_so_si. If None, it will be treated as an unknown.
    q_conv_ia_si_in (float or None): The value of q_conv_ia_si. If None, it will be treated as an unknown.
    q_lwr_i_si_in (float or None): The value of q_lwr_i_si. If None, it will be treated as an unknown.
    Returns:
    list: A list of solutions for the unknown variables.
    """

    q_sol, q_cond_so_si, q_conv_ia_si, q_lwr_i_si = sym('q_sol q_cond_so_si q_conv_ia_si q_lwr_i_si')
    f = Eq(q_sol + q_cond_so_si + q_conv_ia_si + q_lwr_i_si, 0)
    if q_sol_in is not None:
        f = f.subs(q_sol, q_sol_in)
    if q_cond_so_si_in is not None:
        f = f.subs(q_cond_so_si, q_cond_so_si_in)
    if q_conv_ia_si_in is not None:
        f = f.subs(q_conv_ia_si, q_conv_ia_si_in)
    if q_lwr_i_si_in is not None:
        f = f.subs(q_lwr_i_si, q_lwr_i_si_in)
    return solve(f)
