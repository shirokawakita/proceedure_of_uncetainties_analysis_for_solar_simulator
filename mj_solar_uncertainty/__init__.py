"""
多接合太陽電池ソーラシミュレータ不確かさ評価パッケージ.
Baur & Bett (Fraunhofer ISE, 2005) 手法に準拠.
"""

from .core import (
    SpectralCurve,
    common_grid,
    integral_J,
    solve_eq2,
    u_J_eq4,
    u_J_independent,
    u_Ai_analytical,
    monte_carlo_subcell_currents,
    MCResult,
    IVVariation,
    combine_subcells_RSS,
    make_table1,
)
from .io import (
    load_spectral,
    load_reference_spectrum,
    load_astm_e490_am0_xls,
    export_reference_spectrum_csv,
    reference_am0_from_e490_xls_to_grid,
    load_iv_variation,
    save_table1_csv,
    save_table1_md,
)

__version__ = "0.1.0"
__all__ = [
    "SpectralCurve", "common_grid", "integral_J", "solve_eq2",
    "u_J_eq4", "u_J_independent", "u_Ai_analytical",
    "monte_carlo_subcell_currents", "MCResult",
    "IVVariation", "combine_subcells_RSS", "make_table1",
    "load_spectral", "load_reference_spectrum",
    "load_astm_e490_am0_xls", "export_reference_spectrum_csv",
    "reference_am0_from_e490_xls_to_grid",
    "load_iv_variation",
    "save_table1_csv", "save_table1_md",
]
