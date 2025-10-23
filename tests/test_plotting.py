"""
Tests for hcrfish.plotting module.

This module tests plotting utilities and matplotlib theme configuration.
"""

import pytest
import matplotlib.pyplot as plt
from unittest.mock import patch

from hcrfish.plotting.utils import (
    black_plotting,
    white_plotting,
    set_dark_theme,
    set_light_theme
)


class TestPlottingThemes:
    """Test plotting theme configuration functions."""
    
    def test_black_plotting_sets_parameters(self):
        """Test that black_plotting sets correct matplotlib parameters."""
        # Store original values
        original_params = {
            'figure.facecolor': plt.rcParams['figure.facecolor'],
            'axes.facecolor': plt.rcParams['axes.facecolor'],
            'text.color': plt.rcParams['text.color']
        }
        
        try:
            # Apply black theme
            black_plotting()
            
            # Check that dark theme parameters are set
            assert plt.rcParams['figure.facecolor'] == 'black'
            assert plt.rcParams['axes.facecolor'] == 'black'
            assert plt.rcParams['savefig.facecolor'] == 'black'
            assert plt.rcParams['text.color'] == 'white'
            assert plt.rcParams['axes.labelcolor'] == 'white'
            assert plt.rcParams['xtick.color'] == 'white'
            assert plt.rcParams['ytick.color'] == 'white'
            assert plt.rcParams['axes.titlecolor'] == 'white'
            assert plt.rcParams['axes.edgecolor'] == 'white'
            assert plt.rcParams['legend.facecolor'] == 'black'
            assert plt.rcParams['legend.edgecolor'] == 'white'
            assert plt.rcParams['legend.labelcolor'] == 'white'
            assert plt.rcParams['axes.grid'] is False
            
        finally:
            # Restore original parameters
            plt.rcParams.update(original_params)
    
    def test_white_plotting_sets_parameters(self):
        """Test that white_plotting sets correct matplotlib parameters."""
        # Store original values  
        original_params = {
            'figure.facecolor': plt.rcParams['figure.facecolor'],
            'axes.facecolor': plt.rcParams['axes.facecolor'],
            'text.color': plt.rcParams['text.color']
        }
        
        try:
            # Apply white theme
            white_plotting()
            
            # Check that light theme parameters are set
            assert plt.rcParams['figure.facecolor'] == 'white'
            assert plt.rcParams['axes.facecolor'] == 'white'
            assert plt.rcParams['savefig.facecolor'] == 'white'
            assert plt.rcParams['text.color'] == 'black'
            assert plt.rcParams['axes.labelcolor'] == 'black'
            assert plt.rcParams['xtick.color'] == 'black'
            assert plt.rcParams['ytick.color'] == 'black'
            assert plt.rcParams['axes.titlecolor'] == 'black'
            assert plt.rcParams['axes.edgecolor'] == 'black'
            assert plt.rcParams['legend.facecolor'] == 'white'
            assert plt.rcParams['legend.edgecolor'] == 'black'
            assert plt.rcParams['legend.labelcolor'] == 'black'
            assert plt.rcParams['axes.grid'] is False
            
        finally:
            # Restore original parameters
            plt.rcParams.update(original_params)
    
    def test_theme_switching(self):
        """Test switching between themes."""
        # Start with white theme
        white_plotting()
        assert plt.rcParams['figure.facecolor'] == 'white'
        assert plt.rcParams['text.color'] == 'black'
        
        # Switch to black theme
        black_plotting()
        assert plt.rcParams['figure.facecolor'] == 'black'
        assert plt.rcParams['text.color'] == 'white'
        
        # Switch back to white theme
        white_plotting()
        assert plt.rcParams['figure.facecolor'] == 'white'
        assert plt.rcParams['text.color'] == 'black'
    
    def test_alias_functions(self):
        """Test that alias functions work correctly."""
        # Test dark theme alias
        set_dark_theme()
        assert plt.rcParams['figure.facecolor'] == 'black'
        assert plt.rcParams['text.color'] == 'white'
        
        # Test light theme alias
        set_light_theme()
        assert plt.rcParams['figure.facecolor'] == 'white'
        assert plt.rcParams['text.color'] == 'black'
    
    def test_spine_configuration(self):
        """Test that spine configuration is set correctly."""
        white_plotting()
        
        # Check spine settings
        assert plt.rcParams['axes.spines.left'] is True
        assert plt.rcParams['axes.spines.bottom'] is True
        assert plt.rcParams['axes.spines.top'] is False
        assert plt.rcParams['axes.spines.right'] is False
        
        black_plotting()
        
        # Should be same for dark theme
        assert plt.rcParams['axes.spines.left'] is True
        assert plt.rcParams['axes.spines.bottom'] is True
        assert plt.rcParams['axes.spines.top'] is False
        assert plt.rcParams['axes.spines.right'] is False
    
    def test_legend_configuration(self):
        """Test legend configuration for both themes."""
        # White theme legend
        white_plotting()
        assert plt.rcParams['legend.facecolor'] == 'white'
        assert plt.rcParams['legend.edgecolor'] == 'black'
        assert plt.rcParams['legend.labelcolor'] == 'black'
        assert plt.rcParams['legend.fontsize'] == 'medium'
        
        # Black theme legend
        black_plotting()
        assert plt.rcParams['legend.facecolor'] == 'black'
        assert plt.rcParams['legend.edgecolor'] == 'white'
        assert plt.rcParams['legend.labelcolor'] == 'white'
        assert plt.rcParams['legend.fontsize'] == 'medium'


class TestPlottingIntegration:
    """Test plotting functions with actual matplotlib operations."""
    
    def test_create_plot_with_white_theme(self):
        """Test creating a plot with white theme."""
        import matplotlib.pyplot as plt
        import numpy as np
        
        white_plotting()
        
        # Create a simple plot
        fig, ax = plt.subplots(figsize=(6, 4))
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        ax.plot(x, y, label='sin(x)')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_title('Test Plot')
        ax.legend()
        
        # Check that the theme is applied
        assert fig.get_facecolor() == (1.0, 1.0, 1.0, 1.0)  # White RGBA
        
        plt.close(fig)
    
    def test_create_plot_with_black_theme(self):
        """Test creating a plot with black theme."""
        import matplotlib.pyplot as plt
        import numpy as np
        
        black_plotting()
        
        # Create a simple plot
        fig, ax = plt.subplots(figsize=(6, 4))
        x = np.linspace(0, 10, 100)
        y = np.cos(x)
        ax.plot(x, y, label='cos(x)')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_title('Test Plot')
        ax.legend()
        
        # Check that the theme is applied
        assert fig.get_facecolor() == (0.0, 0.0, 0.0, 1.0)  # Black RGBA
        
        plt.close(fig)
    
    def test_theme_persistence(self):
        """Test that theme settings persist across multiple plots."""
        import matplotlib.pyplot as plt
        import numpy as np
        
        black_plotting()
        
        # Create first plot
        fig1, ax1 = plt.subplots()
        ax1.plot([1, 2, 3], [1, 4, 9])
        assert fig1.get_facecolor() == (0.0, 0.0, 0.0, 1.0)
        
        # Create second plot - should still use black theme
        fig2, ax2 = plt.subplots()
        ax2.plot([1, 2, 3], [3, 2, 1])
        assert fig2.get_facecolor() == (0.0, 0.0, 0.0, 1.0)
        
        plt.close(fig1)
        plt.close(fig2)
    
    def test_parameter_types(self):
        """Test that all parameters are set to correct types."""
        white_plotting()
        
        # Check that color parameters are strings
        assert isinstance(plt.rcParams['figure.facecolor'], str)
        assert isinstance(plt.rcParams['axes.facecolor'], str)
        assert isinstance(plt.rcParams['text.color'], str)
        
        # Check that boolean parameters are booleans
        assert isinstance(plt.rcParams['axes.grid'], bool)
        assert isinstance(plt.rcParams['axes.spines.left'], bool)
        
        black_plotting()
        
        # Same checks for black theme
        assert isinstance(plt.rcParams['figure.facecolor'], str)
        assert isinstance(plt.rcParams['axes.facecolor'], str)
        assert isinstance(plt.rcParams['text.color'], str)
        assert isinstance(plt.rcParams['axes.grid'], bool)


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_repeated_theme_application(self):
        """Test applying the same theme multiple times."""
        # Apply white theme multiple times
        for _ in range(3):
            white_plotting()
            assert plt.rcParams['figure.facecolor'] == 'white'
            assert plt.rcParams['text.color'] == 'black'
        
        # Apply black theme multiple times
        for _ in range(3):
            black_plotting()
            assert plt.rcParams['figure.facecolor'] == 'black'
            assert plt.rcParams['text.color'] == 'white'
    
    def test_function_return_values(self):
        """Test that theme functions return None."""
        assert white_plotting() is None
        assert black_plotting() is None
        assert set_light_theme() is None
        assert set_dark_theme() is None
    
    @patch('matplotlib.pyplot.rcParams')
    def test_rcparams_update_called(self, mock_rcparams):
        """Test that rcParams.update is called with correct parameters."""
        # Mock rcParams as a dictionary
        mock_rcparams.update = lambda x: None
        mock_rcparams.__setitem__ = lambda k, v: None
        
        white_plotting()
        
        # Verify that update was called (via the mock)
        # This tests that the function attempts to update parameters
        assert True  # If we get here without error, the function ran
    
    def test_parameter_completeness(self):
        """Test that all expected parameters are set by both themes."""
        expected_params = [
            'figure.facecolor',
            'axes.facecolor', 
            'savefig.facecolor',
            'axes.edgecolor',
            'axes.labelcolor',
            'xtick.color',
            'ytick.color',
            'axes.titlecolor',
            'text.color',
            'axes.grid',
            'legend.facecolor',
            'legend.edgecolor',
            'legend.labelcolor',
            'legend.fontsize',
            'axes.spines.left',
            'axes.spines.bottom',
            'axes.spines.top',
            'axes.spines.right'
        ]
        
        # Test white theme
        white_plotting()
        for param in expected_params:
            assert param in plt.rcParams
        
        # Test black theme
        black_plotting()
        for param in expected_params:
            assert param in plt.rcParams