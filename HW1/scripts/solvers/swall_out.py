import json
import pandas as pd
from sympy import Eq, symbols as sym, solve

def solve_qsol_swall_out(As_in, alpha_in, qdir_in, qdiff_in, qrefl_in, Qsol_in):
    """
    Solves the solar radiation balance equation for a surface.

    The function solves the following equation for the given inputs:
    
    .. math::
        A_s \cdot \alpha \cdot (q_{dir} + q_{diff} + q_{refl}) - Q_{sol} = 0

    where:
    - As: Surface area
    - alpha: Absorptivity of the surface
    - qdir: Direct solar radiation
    - qdiff: Diffuse solar radiation
    - qrefl: Reflected solar radiation
    - Qsol: Total solar radiation absorbed by the surface

    Parameters:
    - As_in (float or None): Input value for surface area (As). If None, it will be treated as a symbolic variable.
    - alpha_in (float or None): Input value for absorptivity (alpha). If None, it will be treated as a symbolic variable.
    - qdir_in (float or None): Input value for direct solar radiation (qdir). If None, it will be treated as a symbolic variable.
    - qdiff_in (float or None): Input value for diffuse solar radiation (qdiff). If None, it will be treated as a symbolic variable.
    - qrefl_in (float or None): Input value for reflected solar radiation (qrefl). If None, it will be treated as a symbolic variable.
    - Qsol_in (float or None): Input value for total solar radiation absorbed (Qsol). If None, it will be treated as a symbolic variable.

    Returns:
    - None: The function prints the solution to the equation.
    """
    
    As, alpha, qdir, qdiff, qrefl, Qsol = sym('As alpha qdir qdiff qrefl Qsol')
    f = Eq(As * alpha * (qdir + qdiff + qrefl) - Qsol, 0)
    if As_in != None:
        f = f.subs(As, As_in)
    if alpha_in != None:
        f = f.subs(alpha, alpha_in)
    if qdir_in != None:
        f = f.subs(qdir, qdir_in)
    if qdiff_in != None:
        f = f.subs(qdiff, qdiff_in)
    if qrefl_in != None:
        f = f.subs(qrefl, qrefl_in)
    if Qsol_in != None:
        f = f.subs(Qsol, Qsol_in)
    print(solve(f))
    return solve(f)

def solve_qconv_oa_so(Asw_in, h_so_in, T_so_in, T_oa_in, Qconv_oa_so_in):
    """
    Solves the convective heat transfer equation between the outside air and a surface.

    The function solves the following equation for the given inputs:
    
    .. math::
        A_{sw} \cdot h_{so} \cdot (T_{so} - T_{oa}) - Q_{conv,oa,so} = 0

    where:
    - Asw: Surface area of the wall
    - h_so: Convective heat transfer coefficient
    - T_so: Surface temperature
    - T_oa: Outside air temperature
    - Qconv_oa_so: Convective heat transfer rate

    Parameters:
    - Asw_in (float or None): Input value for Asw. If None, it will be treated as a symbolic variable.
    - h_so_in (float or None): Input value for h_so. If None, it will be treated as a symbolic variable.
    - T_so_in (float or None): Input value for T_so. If None, it will be treated as a symbolic variable.
    - T_oa_in (float or None): Input value for T_oa. If None, it will be treated as a symbolic variable.
    - Qconv_oa_so_in (float or None): Input value for Qconv_oa_so. If None, it will be treated as a symbolic variable.

    Returns:
    - None: Prints the solution to the equation.
    """
    Asw, h_so, T_so, T_oa, Qconv_oa_so = sym("Asw h_so T_so T_oa Qconv_oa_so")
    f = Eq(Asw * h_so * (T_so - T_oa) - Qconv_oa_so, 0)
    if Asw_in != None:
        f = f.subs(Asw, Asw_in)
    if h_so_in != None:
        f = f.subs(h_so, h_so_in)
    if T_so_in != None:
        f = f.subs(T_so, T_so_in)
    if T_oa_in != None:
        f = f.subs(T_oa, T_oa_in)
    if Qconv_oa_so_in != None:
        f = f.subs(Qconv_oa_so, Qconv_oa_so_in)
    print(solve(f))
    return solve(f)

def solve_qlwr_gr_so(As_in, Fso_gr_in, E_so_in, E_gr_in, sigma_in, T_gr_in, T_so_in, Qlwr_gr_so):
    As, Fso_gr, E_so, E_gr, sigma, T_gr, T_so, Qlwr_gr_so = sym("As Fso_gr E_so E_gr sigma T_gr T_so Qlwr_gr_so")
    f = Eq(As * sigma * (Fso_gr * E_so - E_gr) * (T_so**4 - T_gr**4) - Qlwr_gr_so, 0)
    if As_in != None:
        f = f.subs(As, As_in)
    if Fso_gr_in != None:
        f = f.subs(Fso_gr, Fso_gr_in)
    if E_so_in != None:
        f = f.subs(E_so, E_so_in)
    if E_gr_in != None:
        f = f.subs(E_gr, E_gr_in)
    if sigma_in != None:
        f = f.subs(sigma, sigma_in)
    if T_gr_in != None:
        f = f.subs(T_gr, T_gr_in)
    if T_so_in != None:
        f = f.subs(T_so, T_so_in)
    if Qlwr_gr_so != None:
        f = f.subs(Qlwr_gr_so, Qlwr_gr_so)
    print(solve(f))
    return solve(f)

def solve_lwr_sky_so(As_in, F_so_sky_in, E_so_in, E_gr_in, sigma_in, T_sky_in, T_so_in, Qlwr_sky_so_in):
    """
    Solves the longwave radiation exchange between the soil and the sky.

    The function solves the following equation for the given inputs:
    
    .. math::
        A_s \cdot \sigma \cdot F_{so,sky} \cdot E_{so} \cdot (T_{so}^4 - T_{sky}^4) - Q_{lwr,sky,so} = 0

    where:
    - As: Surface area of the soil
    - F_so_sky: View factor from soil to sky
    - E_so: Emissivity of the soil
    - E_gr: Emissivity of the ground (not used in the equation)
    - sigma: Stefan-Boltzmann constant
    - T_so: Temperature of the soil
    - T_sky: Temperature of the sky
    - Qlwr_sky_so: Longwave radiation from sky to soil

    Parameters:
    - As_in (float or None): Input value for As
    - F_so_sky_in (float or None): Input value for F_so_sky
    - E_so_in (float or None): Input value for E_so
    - E_gr_in (float or None): Input value for E_gr
    - sigma_in (float or None): Input value for sigma
    - T_sky_in (float or None): Input value for T_sky
    - T_so_in (float or None): Input value for T_so
    - Qlwr_sky_so_in (float or None): Input value for Qlwr_sky_so

    Returns:
    - None: Prints the solution to the equation.
    """
    As, F_so_sky, E_so, E_gr, sigma, T_so, T_sky, Qlwr_sky_so = sym("As F_so_sky E_so E_gr sigma T_so T_sky Qlwr_sky_so")
    f = Eq(As * sigma * F_so_sky * E_so * (T_so**4 - T_sky**4) - Qlwr_sky_so, 0)
    if As_in != None:
        f = f.subs(As, As_in)
    if F_so_sky_in != None:
        f = f.subs(F_so_sky, F_so_sky_in)
    if E_so_in != None:
        f = f.subs(E_so, E_so_in)
    if E_gr_in != None:
        f = f.subs(E_gr, E_gr_in)
    if sigma_in != None:
        f = f.subs(sigma, sigma_in)
    if T_sky_in != None:
        f = f.subs(T_sky, T_sky_in)
    if T_so_in != None:
        f = f.subs(T_so, T_so_in)
    if Qlwr_sky_so_in != None:
        f = f.subs(Qlwr_sky_so, Qlwr_sky_so_in)
    print(solve(f))
    return solve(f)

def solve_qcond_si_so(As_in, k_in, L_in, T_si_in, T_so_in, Q_cond_si_so_in):
    """
    Solves for the unknown variable in the heat conduction equation for a wall.

    The function solves the following equation for the given inputs:
    
    .. math::
        Q_{cond,si,so} = \frac{A_s \cdot k \cdot (T_{si}^4 - T_{so}^4)}{L}

    where:
    - As: Surface area (As) in square meters
    - k: Thermal conductivity (k) in W/(m*K)
    - L: Thickness of the wall (L) in meters
    - T_si: Temperature on the inside surface (T_si) in Kelvin
    - T_so: Temperature on the outside surface (T_so) in Kelvin
    - Q_cond_si_so: Heat conduction (Q_cond_si_so) in Watts

    Parameters:
    - As_in (float or None): Input value for As. If None, it will be treated as a symbolic variable.
    - k_in (float or None): Input value for k. If None, it will be treated as a symbolic variable.
    - L_in (float or None): Input value for L. If None, it will be treated as a symbolic variable.
    - T_si_in (float or None): Input value for T_si. If None, it will be treated as a symbolic variable.
    - T_so_in (float or None): Input value for T_so. If None, it will be treated as a symbolic variable.
    - Q_cond_si_so_in (float or None): Input value for Q_cond_si_so. If None, it will be treated as a symbolic variable.

    Returns:
    - None: Prints the solution to the equation.
    """
    
    As, k, L, T_si, T_so, Q_cond_si_so = sym("As k L T_si T_so Q_cond_si_so")
    f = Eq(((As * k * (T_si**4 - T_so**4))/L) - Q_cond_si_so, 0)
    if As_in != None:
        f = f.subs(As, As_in)
    if k_in != None:
        f = f.subs(k, k_in)
    if L_in != None:
        f = f.subs(L, L_in)
    if T_si_in != None:
        f = f.subs(T_si, T_si_in)
    if T_so_in != None:
        f = f.subs(T_so, T_so_in)
    if Q_cond_si_so_in != None:
        f = f.subs(Q_cond_si_so, Q_cond_si_so_in)
    print(solve(f))
    return solve(f)

def solve_ovr_swall_out(Q_sol_in, Q_conv_oa_so, Q_lwr_gr_so, Q_lwr_sky_so, Q_cond_si_so):
    """
    Solves the overall heat balance equation for a surface.

    The function solves the following equation for the given inputs:
    
    .. math::
        Q_{sol} + Q_{conv,oa,so} + Q_{lwr,gr,so} + Q_{lwr,sky,so} + Q_{cond,si,so} = 0

    Parameters:
    - Q_sol_in (float or None): Solar radiation heat gain. If None, it will be treated as a symbolic variable.
    - Q_conv_oa_so (float or None): Convective heat transfer from outdoor air to surface. If None, it will be treated as a symbolic variable.
    - Q_lwr_gr_so (float or None): Longwave radiation heat transfer from ground to surface. If None, it will be treated as a symbolic variable.
    - Q_lwr_sky_so (float or None): Longwave radiation heat transfer from sky to surface. If None, it will be treated as a symbolic variable.
    - Q_cond_si_so (float or None): Conductive heat transfer from surface to interior. If None, it will be treated as a symbolic variable.

    Returns:
    - None: The function prints the solution to the heat balance equation.
    """
    
    Q_sol, Q_conv_oa_so, Q_lwr_gr_so, Q_lwr_sky_so, Q_cond_si_so = sym("Q_sol Q_conv_oa_so Q_lwr_gr_so Q_lwr_sky_so Q_cond_si_so")
    f = Eq(Q_sol + Q_conv_oa_so + Q_lwr_gr_so + Q_lwr_sky_so + Q_cond_si_so, 0)
    if Q_sol_in != None:
        f = f.subs(Q_sol, Q_sol_in)
    if Q_conv_oa_so != None:
        f = f.subs(Q_conv_oa_so, Q_conv_oa_so)
    if Q_lwr_gr_so != None:
        f = f.subs(Q_lwr_gr_so, Q_lwr_gr_so)
    if Q_lwr_sky_so != None:
        f = f.subs(Q_lwr_sky_so, Q_lwr_sky_so)
    if Q_cond_si_so != None:
        f = f.subs(Q_cond_si_so, Q_cond_si_so)
    return solve(f)