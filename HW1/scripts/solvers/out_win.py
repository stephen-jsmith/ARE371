import json
import pandas as pd
from sympy import Eq, symbols as sym, solve

def solve_q_sol_out_win(alpha_glass_in, A_w_in, q_dir_in, q_diff_in, q_refl_in, q_sol):
    """
    Solves the equation for the solar heat gain through a window.

    The function solves the following equation for the given inputs:

    .. math::
        Q_{sol} + Q_{conv,oa,wi} + Q_{lwr,gr,wi} + Q_{lwr,sky,wi} + Q_{cond,si,wi} = 0

    Parameters:
    alpha_glass_in (float or None): The solar absorptance of the glass. If None, it will be treated as a symbolic variable.
    q_dir_in (float or None): Direct solar radiation. If None, it will be treated as a symbolic variable.
    q_diff_in (float or None): Diffuse solar radiation. If None, it will be treated as a symbolic variable.
    q_refl_in (float or None): Reflected solar radiation. If None, it will be treated as a symbolic variable.
    q_sol (float or None): Solar heat gain through the window. If None, it will be treated as a symbolic variable.
    Returns:
    list: Solutions to the equation.
    """

    alpha_glass, A_w, q_dir, q_diff, q_refl, q_sol = sym('alpha_glass A_w q_dir q_diff q_refl q_sol')
    f = Eq(0.5*alpha_glass*A_w(q_dir + q_diff + q_refl) - q_sol, 0)
    var_solved = None
    for var, val in zip([alpha_glass, A_w, q_dir, q_diff, q_refl], [alpha_glass_in, A_w_in, q_dir_in, q_diff_in, q_refl_in]):
        if val is not None:
            f = f.subs(var, val)
        else:
            var_solved = var
    return (var_solved, solve(f))

def solve_q_conv_oa_wo(A_w_in=None, h_wo_in=None, T_oa_in=None, T_wo_in=None, q_conv_oa_wo_in=None):
    """
    Solves the equation for the convective heat transfer between the outdoor air and the window.

    The function solves the following equation for the given inputs:

    .. math::
        h_{wo}A_{w}(T_{oa} - T_{wo}) = Q_{conv,oa,wo}

    Parameters:
    A_w_in (float or None): Area of the window. If None, it will be treated as a symbolic variable.
    h_wo_in (float or None): Convective heat transfer coefficient between the outdoor air and the window. If None, it will be treated as a symbolic variable.
    T_oa_in (float or None): Outdoor air temperature. If None, it will be treated as a symbolic variable.
    T_wo_in (float or None): Window temperature. If None, it will be treated as a symbolic variable.
    q_conv_oa_wo_in (float or None): Convective heat transfer between the outdoor air and the window. If None, it will be treated as a symbolic variable.
    Returns:
    list: Solutions to the equation.
    """

    A_w, h_wo, T_oa, T_wo, q_conv_oa_wo = sym('A_w h_wo T_oa T_wo q_conv_oa_wo')
    f = Eq(h_wo*A_w*(T_oa - T_wo) - q_conv_oa_wo, 0)
    var_solved = None
    for var, val in zip([A_w, h_wo, T_oa, T_wo, q_conv_oa_wo], [A_w_in, h_wo_in, T_oa_in, T_wo_in, q_conv_oa_wo_in]):
        if val is not None:
            f = f.subs(var, val)
        else:
            var_solved = var
    return [solve(f), var_solved]

def solve_q_lwr_gr_wo(A_w_in=None, F_wo_gr_in=None, E_wo_in=None, E_gr_in=None, sigma_in=None, T_gr_in=None, T_wo_in=None, q_lwr_gr_wo_in=None):
    """
    Solves the equation for the longwave radiation between the ground and the window.

    The function solves the following equation for the given inputs:

    .. math::
        F_{wo,gr}E_{wo}E_{gr}\\sigma(T_{gr}^4 - T_{wo}^4) = Q_{lwr,gr,wo}

    Parameters:
    A_w_in (float or None): Area of the window. If None, it will be treated as a symbolic variable.
    F_wo_gr_in (float or None): View factor between the window and the ground. If None, it will be treated as a symbolic variable.
    E_wo_in (float or None): Emissivity of the window. If None, it will be treated as a symbolic variable.
    E_gr_in (float or None): Emissivity of the ground. If None, it will be treated as a symbolic variable.
    sigma_in (float or None): Stefan-Boltzmann constant. If None, it will be treated as a symbolic variable.
    T_gr_in (float or None): Ground temperature. If None, it will be treated as a symbolic variable.
    T_wo_in (float or None): Window temperature. If None, it will be treated as a symbolic variable.
    q_lwr_gr_wo_in (float or None): Longwave radiation between the ground and the window. If None, it will be treated as a symbolic variable.
    Returns:
    list: Solutions to the equation.
    """

    A_w, F_wo_gr, E_wo, E_gr, sigma, T_gr, T_wo, q_lwr_gr_wo = sym('A_w F_wo_gr E_wo E_gr sigma T_gr T_wo q_lwr_gr_wo')
    f = Eq(A_w*F_wo_gr*E_wo*E_gr*sigma*(T_gr**4 - T_wo**4) - q_lwr_gr_wo, 0)
    var_solved = None
    for var, val in zip([A_w, F_wo_gr, E_wo, E_gr, sigma, T_gr, T_wo, q_lwr_gr_wo], 
                        [A_w_in, F_wo_gr_in, E_wo_in, E_gr_in, sigma_in, T_gr_in, T_wo_in, q_lwr_gr_wo_in]):
        if val is not None:
            f = f.subs(var, val)
        else:
            var_solved = var
    return (var_solved,solve(f))

def solve_q_lwr_sky_wo(A_w_in, F_wo_sky_in, E_wo_in, E_sky_in, sigma_in, T_sky_in, T_wo_in, q_lwr_sky_wo_in):
    """
    Solves the equation for the longwave radiation between the sky and the window.

    The function solves the following equation for the given inputs:

    .. math::
        F_{wo,sky}E_{wo}E_{sky}\\sigma(T_{sky}^4 - T_{wo}^4) = Q_{lwr,sky,wo}

    Parameters:
    A_w_in (float or None): Area of the window. If None, it will be treated as a symbolic variable.
    F_wo_sky_in (float or None): View factor between the window and the sky. If None, it will be treated as a symbolic variable.
    E_wo_in (float or None): Emissivity of the window. If None, it will be treated as a symbolic variable.
    E_sky_in (float or None): Emissivity of the sky. If None, it will be treated as a symbolic variable.
    sigma_in (float or None): Stefan-Boltzmann constant. If None, it will be treated as a symbolic variable.
    T_sky_in (float or None): Sky temperature. If None, it will be treated as a symbolic variable.
    T_wo_in (float or None): Window temperature. If None, it will be treated as a symbolic variable.
    q_lwr_sky_wo_in (float or None): Longwave radiation between the sky and the window. If None, it will be treated as a symbolic variable.
    Returns:
    list: Solutions to the equation.
    """

    A_w, F_wo_sky, E_wo, E_sky, sigma, T_sky, T_wo, q_lwr_sky_wo = sym('A_w F_wo_sky E_wo E_sky sigma T_sky T_wo q_lwr_sky_wo')
    f = Eq(A_w*F_wo_sky*E_wo*E_sky*sigma*(T_sky**4 - T_wo**4) - q_lwr_sky_wo, 0)
    solved_var = None
    for var, val in zip([A_w, F_wo_sky, E_wo, E_sky, sigma, T_sky, T_wo, q_lwr_sky_wo], 
                        [A_w_in, F_wo_sky_in, E_wo_in, E_sky_in, sigma_in, T_sky_in, T_wo_in, q_lwr_sky_wo_in]):
        if val is not None:
            f = f.subs(var, val)
        else:
            solved_var = var
    return [solved_var, solve(f)]

def solve_q_cond_wi_wo(A_w_in, k_in, L_in, T_wi_in, T_wo_in, q_cond_wi_wo):
    """
    Solves the equation for the conductive heat transfer between the window and the outdoor air.

    The function solves the following equation for the given inputs:

    .. math::
        A_w\\frac{k}{L}(T_{wi} - T_{wo}) = Q_{cond,wi,wo}

    Parameters:
    A_w_in (float or None): Area of the window. If None, it will be treated as a symbolic variable.
    k_in (float or None): Thermal conductivity of the window. If None, it will be treated as a symbolic variable.
    L_in (float or None): Thickness of the window. If None, it will be treated as a symbolic variable.
    T_wi_in (float or None): Window temperature. If None, it will be treated as a symbolic variable.
    T_wo_in (float or None): Outdoor air temperature. If None, it will be treated as a symbolic variable.
    q_cond_wi_wo (float or None): Conductive heat transfer between the window and the outdoor air. If None, it will be treated as a symbolic variable.
    Returns:
    list: Solutions to the equation.
    """

    A_w, k, L, T_wi, T_wo, q_cond_wi_wo = sym('A_w k L T_wi T_wo q_cond_wi_wo')
    f = Eq(A_w*k/L*(T_wi - T_wo) - q_cond_wi_wo, 0)
    solved_var = None
    for var, val in zip([A_w, k, L, T_wi, T_wo, q_cond_wi_wo], [A_w_in, k_in, L_in, T_wi_in, T_wo_in, q_cond_wi_wo]):
        if val is not None:
            f = f.subs(var, val)
        else:
            solved_var = var
    return [solved_var, solve(f)]

def solve_ovr_out_win(q_sol_in, q_conv_oa_wo_in, q_lwr_sky_wo_in, q_cond_wi_wo):
    """
    Solves the overall equation for the heat transfer through a window.

    The function solves the following equation for the given inputs:

    .. math::
        Q_{sol} + Q_{conv,oa,wo} + Q_{lwr,sky,wo} + Q_{cond,wi,wo} = 0

    Parameters:
    q_sol_in (float or None): Solar heat gain through the window. If None, it will be treated as a symbolic variable.
    q_conv_oa_wo_in (float or None): Convective heat transfer between the outdoor air and the window. If None, it will be treated as a symbolic variable.
    q_lwr_sky_wo_in (float or None): Longwave radiation between the sky and the window. If None, it will be treated as a symbolic variable.
    q_cond_wi_wo (float or None): Conductive heat transfer between the window and the outdoor air. If None, it will be treated as a symbolic variable.
    Returns:
    list: Solutions to the equation.
    """

    q_sol, q_conv_oa_wo, q_lwr_sky_wo, q_cond_wi_wo = sym('q_sol q_conv_oa_wo q_lwr_sky_wo q_cond_wi_wo')
    f = Eq(q_sol + q_conv_oa_wo + q_lwr_sky_wo + q_cond_wi_wo, 0)
    solved_var = None
    if q_sol_in != None:
        for var, val in zip([q_sol, q_conv_oa_wo, q_lwr_sky_wo, q_cond_wi_wo], [q_sol_in, q_conv_oa_wo_in, q_lwr_sky_wo_in, q_cond_wi_wo]):
            if val is not None:
                f = f.subs(var, val)
            else:
                solved_var = var
        return [solved_var, solve(f)]
    return solve(f)