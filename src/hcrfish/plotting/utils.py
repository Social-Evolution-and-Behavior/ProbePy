"""
Plotting utilities for HCR-FISH probe visualization.

This module provides functions to set matplotlib themes for publication-ready
figures with both light and dark color schemes.
"""

import matplotlib.pyplot as plt
from typing import Dict, Any


def black_plotting() -> None:
    """
    Configure matplotlib for dark theme plotting.
    
    Sets up a dark theme with white text and borders on black backgrounds,
    suitable for presentations or dark-themed publications.
    
    Examples:
        >>> import matplotlib.pyplot as plt
        >>> from hcrfish.plotting import black_plotting
        >>> black_plotting()
        >>> plt.plot([1, 2, 3], [1, 4, 9])
        >>> plt.title("Dark Theme Plot")
        >>> plt.show()
    """
    dark_params: Dict[str, Any] = {
        # Background colors
        'figure.facecolor': 'black',
        'axes.facecolor': 'black', 
        'savefig.facecolor': 'black',
        
        # Text and line colors
        'axes.edgecolor': 'white',
        'axes.labelcolor': 'white',
        'xtick.color': 'white',
        'ytick.color': 'white', 
        'axes.titlecolor': 'white',
        'text.color': 'white',
        
        # Grid and legend
        'axes.grid': False,
        'legend.facecolor': 'black',
        'legend.edgecolor': 'white',
        'legend.labelcolor': 'white',
        'legend.fontsize': 'medium',
        
        # Spines
        'axes.spines.left': True,
        'axes.spines.bottom': True,
        'axes.spines.top': False,
        'axes.spines.right': False
    }
    
    plt.rcParams.update(dark_params)


def white_plotting() -> None:
    """
    Configure matplotlib for light theme plotting.
    
    Sets up a clean light theme with black text on white backgrounds,
    suitable for publications and general use.
    
    Examples:
        >>> import matplotlib.pyplot as plt
        >>> from hcrfish.plotting import white_plotting
        >>> white_plotting()
        >>> plt.plot([1, 2, 3], [1, 4, 9])
        >>> plt.title("Light Theme Plot")
        >>> plt.show()
    """
    light_params: Dict[str, Any] = {
        # Background colors
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white',
        
        # Text and line colors  
        'axes.edgecolor': 'black',
        'axes.labelcolor': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black',
        'axes.titlecolor': 'black', 
        'text.color': 'black',
        
        # Grid and legend
        'axes.grid': False,
        'legend.facecolor': 'white',
        'legend.edgecolor': 'black', 
        'legend.labelcolor': 'black',
        'legend.fontsize': 'medium',
        
        # Spines
        'axes.spines.left': True,
        'axes.spines.bottom': True,
        'axes.spines.top': False,
        'axes.spines.right': False
    }
    
    plt.rcParams.update(light_params)


# Alias functions for backwards compatibility
def set_dark_theme() -> None:
    """Alias for black_plotting() for clearer naming."""
    black_plotting()


def set_light_theme() -> None: 
    """Alias for white_plotting() for clearer naming."""
    white_plotting()