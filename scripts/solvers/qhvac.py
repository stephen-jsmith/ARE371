import json
import pandas as pd
from sympy import Eq, symbols as sym, solve

def solve_q_inf_oa_ia(ACH_in, V_in, rho_in, C_pair_in, T_oa_in, T_ia_in, q_inf_oa_ia_in):
    """
    Solves the equation for the infiltration heat gain from outdoor air to indoor air.
    
    The function solves the following equation for the infiltration heat gain from outdoor air to indoor air:
    
    .. math::
        q_{\text{inf\_oa\_ia}} = \frac{ACH \cdot V \cdot \rho \cdot C_{\text{pair}} \cdot (T_{\text{oa}} - T_{\text{ia}})}{3600}
    
    Parameters:
    ACH_in (float or None): The air changes per hour (ACH). If None, it will be treated as a symbolic variable.
    V_in (float or None): The volume of the room (V). If None, it will be treated as a symbolic variable.
    rho_in (float or None): The density of air (rho). If None, it will be treated as a symbolic variable.
    C_pair_in (float or None): The specific heat capacity of air (C_pair). If None, it will be treated as a symbolic variable.
    T_oa_in (float or None): The temperature of outdoor air (T_oa). If None, it will be treated as a symbolic variable.
    T_ia_in (float or None): The temperature of indoor air (T_ia). If None, it will be treated as a symbolic variable.
    q_inf_oa_ia_in (float or None): The infiltration heat gain from outdoor air to indoor air (q_inf_oa_ia). If None, it will be treated as a symbolic variable.
    
    Returns:
    list: A list of solutions for the equation.
    """
    ACH, V, rho, C_pair, T_oa, T_ia, q_inf_oa_ia = sym('ACH V rho C_pair T_oa T_ia q_inf_oa_ia')
    f = Eq(ACH * V * rho * C_pair * (T_oa - T_ia) / 3600 - q_inf_oa_ia, 0)
    if ACH_in != None:
        f.subs(ACH, ACH_in)
    if V_in != None:
        f.subs(V, V_in)
    if rho_in != None:
        f.subs(rho, rho_in)
    if C_pair_in != None:
        f.subs(C_pair, C_pair_in)
    if T_oa_in != None:
        f.subs(T_oa, T_oa_in)
    if T_ia_in != None:
        f.subs(T_ia, T_ia_in)
    if q_inf_oa_ia_in != None:
        f.subs(q_inf_oa_ia, q_inf_oa_ia_in)
    return solve(f)

def solve_q_conv_wi_ia(A_w_in, h_w_in, T_wi_in, T_ia_in, q_conv_wi_ia_in):
    """
    Solves the equation for convective heat transfer between a wall with insulation and indoor air.
    
    The function solves the following equation for convective heat transfer between a wall with insulation and indoor air:
    
    .. math::
        q_{\text{conv\_wi\_ia}} = h_w \cdot A_w \cdot (T_{\text{wi}} - T_{\text{ia}})
    
    Parameters:
    A_w_in (float or None): The area of the wall (A_w). If None, it will be treated as a symbolic variable.
    h_w_in (float or None): The convective heat transfer coefficient of the wall (h_w). If None, it will be treated as a symbolic variable.
    T_wi_in (float or None): The temperature of the wall with insulation (T_wi). If None, it will be treated as a symbolic variable.
    T_ia_in (float or None): The temperature of indoor air (T_ia). If None, it will be treated as a symbolic variable.
    q_conv_wi_ia_in (float or None): The convective heat transfer between the wall with insulation and indoor air (q_conv_wi_ia). If None, it will be treated as a symbolic variable.
    
    Returns:
    list: A list of solutions for the equation.
    """
    A_w, h_w, T_wi, T_ia, q_conv_wi_ia = sym('A_w h_w T_wi T_ia q_conv_wi_ia')
    f = Eq(h_w * A_w * (T_wi - T_ia) - q_conv_wi_ia, 0)
    if A_w_in != None:
        f.subs(A_w, A_w_in)
    if h_w_in != None:
        f.subs(h_w, h_w_in)
    if T_wi_in != None:
        f.subs(T_wi, T_wi_in)
    if T_ia_in != None:
        f.subs(T_ia, T_ia_in)
    if q_conv_wi_ia_in != None:
        f.subs(q_conv_wi_ia, q_conv_wi_ia_in)
    return solve(f)

def solve_q_conv_si_ia(A_s_in, h_s_in, T_si_in, T_ia_in, q_conv_si_ia_in):
    """
    Solves the equation for convective heat transfer between a surface and indoor air.
    
    The function solves the following equation for convective heat transfer between a surface and indoor air:
    
    .. math::
        q_{\text{conv\_si\_ia}} = h_s \cdot A_s \cdot (T_{\text{si}} - T_{\text{ia}})
    
    Parameters:
    A_s_in (float or None): The area of the surface (A_s). If None, it will be treated as a symbolic variable.
    h_s_in (float or None): The convective heat transfer coefficient of the surface (h_s). If None, it will be treated as a symbolic variable.
    T_si_in (float or None): The temperature of the surface (T_si). If None, it will be treated as a symbolic variable.
    T_ia_in (float or None): The temperature of indoor air (T_ia). If None, it will be treated as a symbolic variable.
    q_conv_si_ia_in (float or None): The convective heat transfer between the surface and indoor air (q_conv_si_ia). If None, it will be treated as a symbolic variable.
    
    Returns:
    list: A list of solutions for the equation.
    """
    A_s, h_s, T_si, T_ia, q_conv_si_ia = sym('A_s h_s T_si T_ia q_conv_si_ia')
    f = Eq(h_s * A_s * (T_si - T_ia) - q_conv_si_ia, 0)
    if A_s_in != None:
        f.subs(A_s, A_s_in)
    if h_s_in != None:
        f.subs(h_s, h_s_in)
    if T_si_in != None:
        f.subs(T_si, T_si_in)
    if T_ia_in != None:
        f.subs(T_ia, T_ia_in)
    if q_conv_si_ia_in != None:
        f.subs(q_conv_si_ia, q_conv_si_ia_in)
    return solve(f)

def solve_q_conv_i_ia(A_i_in, h_in, T_i_in, T_ia_in, q_conv_i_ia_in):
    """
    Solves the equation for convective heat transfer between an indoor surface and indoor air.
    
    The function solves the following equation for convective heat transfer between an indoor surface and indoor air:
    
    .. math::
        q_{\text{conv\_i\_ia}} = h \cdot A_i \cdot (T_i - T_{\text{ia}})
    
    Parameters:
    A_i_in (float or None): The area of the indoor surface (A_i). If None, it will be treated as a symbolic variable.
    h_in (float or None): The convective heat transfer coefficient (h). If None, it will be treated as a symbolic variable.
    T_i_in (float or None): The temperature of the indoor surface (T_i). If None, it will be treated as a symbolic variable.
    T_ia_in (float or None): The temperature of indoor air (T_ia). If None, it will be treated as a symbolic variable.
    q_conv_i_ia_in (float or None): The convective heat transfer between the indoor surface and indoor air (q_conv_i_ia). If None, it will be treated as a symbolic variable.
    
    Returns:
    list: A list of solutions for the equation.
    """
    A_i, h, T_i, T_ia, q_conv_i_ia = sym('A_i h T_i T_ia q_conv_i_ia')
    f = Eq(h * A_i * (T_i - T_ia) - q_conv_i_ia, 0)
    if A_i_in != None:
        f.subs(A_i, A_i_in)
    if h_in != None:
        f.subs(h, h_in)
    if T_i_in != None:
        f.subs(T_i, T_i_in)
    if T_ia_in != None:
        f.subs(T_ia, T_ia_in)
    if q_conv_i_ia_in != None:
        f.subs(q_conv_i_ia, q_conv_i_ia_in)
    return solve(f)

def solve_ovr_qhvac(q_inf_oa_ia_in, q_conv_wi_ia_in, q_conv_si_ia_in, q_conv_i_ia_in, q_hvac_in):
    """
    Solves the equation for the overall heat gain/loss of a room.
    
    The function solves the following equation for the overall heat gain/loss of a room:
    
    .. math::
        q_{\text{hvac}} = q_{\text{inf\_oa\_ia}} + q_{\text{conv\_wi\_ia}} + q_{\text{conv\_si\_ia}} + q_{\text{conv\_i\_ia}}
    
    Parameters:
    q_inf_oa_ia_in (float or None): The infiltration heat gain from outdoor air to indoor air (q_inf_oa_ia). If None, it will be treated as a symbolic variable.
    q_conv_wi_ia_in (float or None): The convective heat transfer between the wall with insulation and indoor air (q_conv_wi_ia). If None, it will be treated as a symbolic variable.
    q_conv_si_ia_in (float or None): The convective heat transfer between the surface and indoor air (q_conv_si_ia). If None, it will be treated as a symbolic variable.
    q_conv_i_ia_in (float or None): The convective heat transfer between the indoor surface and indoor air (q_conv_i_ia). If None, it will be treated as a symbolic variable.
    q_hvac_in (float or None): The overall heat gain/loss of the room (q_hvac). If None, it will be treated as a symbolic variable.
    
    Returns:
    list: A list of solutions for the equation.
    """
    q_inf_oa_ia, q_conv_wi_ia, q_conv_si_ia, q_conv_i_ia, q_hvac = sym('q_inf_oa_ia q_conv_wi_ia q_conv_si_ia q_conv_i_ia q_hvac')
    f = Eq(q_inf_oa_ia + q_conv_wi_ia + q_conv_si_ia + q_conv_i_ia - q_hvac, 0)
    if q_inf_oa_ia_in != None:
        f.subs(q_inf_oa_ia, q_inf_oa_ia_in)
    if q_conv_wi_ia_in != None:
        f.subs(q_conv_wi_ia, q_conv_wi_ia_in)
    if q_conv_si_ia_in != None:
        f.subs(q_conv_si_ia, q_conv_si_ia_in)
    if q_conv_i_ia_in != None:
        f.subs(q_conv_i_ia, q_conv_i_ia_in)
    if q_hvac_in != None:
        f.subs(q_hvac, q_hvac_in)
    return solve(f)

def solve_hw(h_in, T_ia_in, T_wi_in, h_w_in):
    """
    Solves the equation for the heat transfer between a wall with insulation and indoor air.
    
    The function solves the following equation for the heat transfer between a wall with insulation and indoor air:
    
    .. math::
        
    
    Parameters:
    h_in (float or None): The convective heat transfer coefficient (h). If None, it will be treated as a symbolic variable.
    T_ia_in (float or None): The temperature of indoor air (T_ia). If None, it will be treated as a symbolic variable.
    T_wi_in (float or None): The temperature of the wall with insulation (T_wi). If None, it will be treated as a symbolic variable.
    
    Returns:
    list: A list of solutions for the equation.
    """
    h, T_ia, T_wi h_w = sym('h T_ia T_wi')
    f = Eq((1.8/h)*abs(T_ia-T_wi)**0.29 - h_w, 0)
    if h_in != None:
        f.subs(h, h_in)
    if T_ia_in != None:
        f.subs(T_ia, T_ia_in)
    if T_wi_in != None:
        f.subs(T_wi, T_wi_in)
    if h_w_in != None:
        f.subs(h_w, h_w_in)
    return solve(f)

def solve_h_s(h_in, T_ia_in, T_si_in, h_s_in)
    """
    Solves the equation for the heat transfer between a surface and indoor air.
    
    The function solves the following equation for the heat transfer between a surface and indoor air:
    
    .. math::
        
    
    Parameters:
    h_in (float or None): The convective heat transfer coefficient (h). If None, it will be treated as a symbolic variable.
    T_ia_in (float or None): The temperature of indoor air (T_ia). If None, it will be treated as a symbolic variable.
    T_si_in (float or None): The temperature of the surface (T_si). If None, it will be treated as a symbolic variable.
    
    Returns:
    list: A list of solutions for the equation.
    """
    h, T_ia, T_si, h_s = sym('h T_ia T_si h_s')
    f = Eq((1.8/h)*abs(T_ia-T_si)**0.29 - h_s, 0)
    if h_in != None:
        f.subs(h, h_in)
    if T_ia_in != None:
        f.subs(T_ia, T_ia_in)
    if T_si_in != None:
        f.subs(T_si, T_si_in)
    if h_s_in != None:
        f.subs(h_s, h_s_in)
    return solve(f)