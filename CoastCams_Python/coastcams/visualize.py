"""
Visualization tools for CoastCams analysis results.

This module provides plotting functions for timestack images,
wave parameters, bathymetry, and time series analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from typing import Optional, List, Dict
from datetime import datetime
import os


class CoastCamsVisualizer:
    """
    Create visualizations for CoastCams analysis results.

    Provides comprehensive plotting functions for all analysis outputs.
    """

    def __init__(self, output_dir: str = 'Output', save_plots: bool = True):
        """
        Initialize visualizer.

        Parameters
        ----------
        output_dir : str, optional
            Directory to save plots (default: 'Output')
        save_plots : bool, optional
            Automatically save plots (default: True)
        """
        self.output_dir = output_dir
        self.save_plots = save_plots

        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")

        print(f"Visualizer initialized: output_dir={output_dir}, save_plots={save_plots}")

        # Set default style
        plt.style.use('seaborn-v0_8-darkgrid' if 'seaborn-v0_8-darkgrid' in plt.style.available else 'default')

    def plot_timestack(self, timestack: np.ndarray,
                      timestamps: Optional[List[datetime]] = None,
                      cross_shore_positions: Optional[np.ndarray] = None,
                      title: str = 'Timestack Image',
                      filename: str = 'timestack.png'):
        """
        Plot 2D timestack image.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (space x time)
        timestamps : List[datetime], optional
            Time axis labels
        cross_shore_positions : np.ndarray, optional
            Cross-shore position labels
        title : str, optional
            Plot title
        filename : str, optional
            Output filename
        """
        fig, ax = plt.subplots(figsize=(12, 6))

        # Plot timestack
        im = ax.imshow(timestack, aspect='auto', cmap='viridis', origin='lower')

        ax.set_xlabel('Time', fontsize=12)
        ax.set_ylabel('Cross-shore Position (pixels)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Intensity', fontsize=11)

        plt.tight_layout()

        if self.save_plots:
            try:
                save_path = os.path.join(self.output_dir, filename)
                print(f"Saving plot: {save_path}")
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Successfully saved: {filename}")
            except Exception as e:
                print(f"ERROR saving {filename}: {e}")
                import traceback
                traceback.print_exc()

        plt.close()

    def plot_shoreline(self, image: np.ndarray,
                      shoreline: np.ndarray,
                      title: str = 'Shoreline Detection',
                      filename: str = 'shoreline.png'):
        """
        Plot image with detected shoreline overlay.

        Parameters
        ----------
        image : np.ndarray
            Original image
        shoreline : np.ndarray
            Shoreline positions
        title : str, optional
            Plot title
        filename : str, optional
            Output filename
        """
        fig, ax = plt.subplots(figsize=(10, 8))

        # Display image
        ax.imshow(image)

        # Overlay shoreline
        x_coords = np.arange(len(shoreline))
        valid_mask = ~np.isnan(shoreline)

        ax.plot(x_coords[valid_mask], shoreline[valid_mask],
               'r-', linewidth=2, label='Detected Shoreline')

        ax.set_xlabel('Cross-shore (pixels)', fontsize=12)
        ax.set_ylabel('Along-shore (pixels)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(fontsize=11)

        plt.tight_layout()

        if self.save_plots:
            try:
                save_path = os.path.join(self.output_dir, filename)
                print(f"Saving plot: {save_path}")
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Successfully saved: {filename}")
            except Exception as e:
                print(f"ERROR saving {filename}: {e}")
                import traceback
                traceback.print_exc()

        plt.close()

    def plot_wave_parameters(self, timestamps: List[datetime],
                            wave_heights: np.ndarray,
                            wave_periods: np.ndarray,
                            wave_celerities: Optional[np.ndarray] = None,
                            title: str = 'Wave Parameters',
                            filename: str = 'wave_parameters.png'):
        """
        Plot time series of wave parameters.

        Parameters
        ----------
        timestamps : List[datetime]
            Time axis
        wave_heights : np.ndarray
            Significant wave heights
        wave_periods : np.ndarray
            Wave periods
        wave_celerities : np.ndarray, optional
            Wave celerities
        title : str, optional
            Plot title
        filename : str, optional
            Output filename
        """
        n_plots = 3 if wave_celerities is not None else 2

        fig, axes = plt.subplots(n_plots, 1, figsize=(12, 4*n_plots))

        if n_plots == 2:
            axes = [axes[0], axes[1]]
        else:
            axes = [axes[0], axes[1], axes[2]]

        # Wave height
        axes[0].plot(timestamps, wave_heights, 'b-', linewidth=1.5)
        axes[0].set_ylabel('Wave Height (m)', fontsize=12)
        axes[0].set_title(title, fontsize=14, fontweight='bold')
        axes[0].grid(True, alpha=0.3)

        # Wave period
        axes[1].plot(timestamps, wave_periods, 'g-', linewidth=1.5)
        axes[1].set_ylabel('Wave Period (s)', fontsize=12)
        axes[1].grid(True, alpha=0.3)

        # Wave celerity
        if wave_celerities is not None:
            axes[2].plot(timestamps, wave_celerities, 'r-', linewidth=1.5)
            axes[2].set_ylabel('Wave Celerity (m/s)', fontsize=12)
            axes[2].set_xlabel('Time', fontsize=12)
            axes[2].grid(True, alpha=0.3)
        else:
            axes[1].set_xlabel('Time', fontsize=12)

        # Format x-axis
        for ax in axes:
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

        plt.tight_layout()

        if self.save_plots:
            try:
                save_path = os.path.join(self.output_dir, filename)
                print(f"Saving plot: {save_path}")
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Successfully saved: {filename}")
            except Exception as e:
                print(f"ERROR saving {filename}: {e}")
                import traceback
                traceback.print_exc()

        plt.close()

    def plot_bathymetry(self, cross_shore_positions: np.ndarray,
                       depths: np.ndarray,
                       reference_depths: Optional[np.ndarray] = None,
                       title: str = 'Bathymetry Profile',
                       filename: str = 'bathymetry.png'):
        """
        Plot bathymetry (depth) profile.

        Parameters
        ----------
        cross_shore_positions : np.ndarray
            Cross-shore positions (meters)
        depths : np.ndarray
            Estimated depths
        reference_depths : np.ndarray, optional
            Reference depths for comparison
        title : str, optional
            Plot title
        filename : str, optional
            Output filename
        """
        fig, ax = plt.subplots(figsize=(12, 6))

        # Plot estimated bathymetry
        ax.plot(cross_shore_positions, -depths, 'b-', linewidth=2,
               label='Estimated Depth')

        # Plot reference if provided
        if reference_depths is not None:
            ax.plot(cross_shore_positions, -reference_depths, 'r--', linewidth=2,
                   label='Reference Depth')

        # Add water surface
        ax.axhline(y=0, color='cyan', linestyle='--', linewidth=1.5,
                  alpha=0.7, label='Water Surface')

        ax.set_xlabel('Cross-shore Distance (m)', fontsize=12)
        ax.set_ylabel('Elevation (m)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)

        # Invert y-axis so depth increases downward
        ax.invert_yaxis()

        plt.tight_layout()

        if self.save_plots:
            try:
                save_path = os.path.join(self.output_dir, filename)
                print(f"Saving plot: {save_path}")
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Successfully saved: {filename}")
            except Exception as e:
                print(f"ERROR saving {filename}: {e}")
                import traceback
                traceback.print_exc()

        plt.close()

    def plot_comprehensive_analysis(self, results: Dict,
                                   timestamps: List[datetime],
                                   filename: str = 'comprehensive_analysis.png'):
        """
        Create comprehensive multi-panel plot of all results.

        Parameters
        ----------
        results : Dict
            Dictionary containing all analysis results
        timestamps : List[datetime]
            Time axis
        filename : str, optional
            Output filename
        """
        fig = plt.figure(figsize=(16, 12))
        gs = GridSpec(3, 2, figure=fig)

        # 1. Timestack
        if 'timestack' in results:
            ax1 = fig.add_subplot(gs[0, :])
            im = ax1.imshow(results['timestack'], aspect='auto',
                          cmap='viridis', origin='lower')
            ax1.set_title('Timestack Image', fontsize=12, fontweight='bold')
            ax1.set_xlabel('Time')
            ax1.set_ylabel('Cross-shore Position')
            plt.colorbar(im, ax=ax1, label='Intensity')

        # 2. Wave height (use actual time series, not mean)
        if 'wave_heights_timeseries' in results and len(timestamps) > 0:
            ax2 = fig.add_subplot(gs[1, 0])
            wave_heights_ts = results['wave_heights_timeseries']
            ax2.plot(timestamps, wave_heights_ts, 'b-', linewidth=1.5, marker='o')
            ax2.set_title('Significant Wave Height', fontsize=12, fontweight='bold')
            ax2.set_ylabel('Hs (m)')
            ax2.grid(True, alpha=0.3)
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

        # 3. Wave period (use actual time series, not mean)
        if 'wave_periods_timeseries' in results and len(timestamps) > 0:
            ax3 = fig.add_subplot(gs[1, 1])
            wave_periods_ts = results['wave_periods_timeseries']
            ax3.plot(timestamps, wave_periods_ts, 'g-', linewidth=1.5, marker='o')
            ax3.set_title('Wave Period', fontsize=12, fontweight='bold')
            ax3.set_ylabel('T (s)')
            ax3.grid(True, alpha=0.3)
            ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45)

        # 4. Bathymetry
        if 'depths' in results and 'cross_shore_positions' in results:
            ax4 = fig.add_subplot(gs[2, 0])
            ax4.plot(results['cross_shore_positions'], -results['depths'],
                    'b-', linewidth=2)
            ax4.axhline(y=0, color='cyan', linestyle='--', linewidth=1.5, alpha=0.7)
            ax4.set_title('Bathymetry Profile', fontsize=12, fontweight='bold')
            ax4.set_xlabel('Cross-shore Distance (m)')
            ax4.set_ylabel('Elevation (m)')
            ax4.grid(True, alpha=0.3)
            ax4.invert_yaxis()

        # 5. Shoreline position
        if 'shoreline_positions' in results:
            ax5 = fig.add_subplot(gs[2, 1])
            ax5.plot(timestamps, results['shoreline_positions'], 'r-', linewidth=1.5)
            ax5.set_title('Shoreline Position', fontsize=12, fontweight='bold')
            ax5.set_xlabel('Time')
            ax5.set_ylabel('Position (m)')
            ax5.grid(True, alpha=0.3)
            ax5.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            plt.setp(ax5.xaxis.get_majorticklabels(), rotation=45)

        plt.tight_layout()

        if self.save_plots:
            try:
                save_path = os.path.join(self.output_dir, filename)
                print(f"Saving plot: {save_path}")
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Successfully saved: {filename}")
            except Exception as e:
                print(f"ERROR saving {filename}: {e}")
                import traceback
                traceback.print_exc()

        plt.close()

    def plot_spectral_analysis(self, frequencies: np.ndarray,
                              psd: np.ndarray,
                              title: str = 'Wave Spectrum',
                              filename: str = 'spectrum.png'):
        """
        Plot power spectral density.

        Parameters
        ----------
        frequencies : np.ndarray
            Frequency axis
        psd : np.ndarray
            Power spectral density
        title : str, optional
            Plot title
        filename : str, optional
            Output filename
        """
        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(frequencies, psd, 'b-', linewidth=1.5)
        ax.set_xlabel('Frequency (Hz)', fontsize=12)
        ax.set_ylabel('Power Spectral Density', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, max(frequencies))

        plt.tight_layout()

        if self.save_plots:
            try:
                save_path = os.path.join(self.output_dir, filename)
                print(f"Saving plot: {save_path}")
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Successfully saved: {filename}")
            except Exception as e:
                print(f"ERROR saving {filename}: {e}")
                import traceback
                traceback.print_exc()

        plt.close()

    def plot_correlation_matrix(self, correlation: np.ndarray,
                                title: str = '2D Correlation',
                                filename: str = 'correlation_2d.png'):
        """
        Plot 2D correlation matrix.

        Parameters
        ----------
        correlation : np.ndarray
            2D correlation matrix
        title : str, optional
            Plot title
        filename : str, optional
            Output filename
        """
        fig, ax = plt.subplots(figsize=(10, 8))

        im = ax.imshow(correlation, cmap='RdBu_r', origin='lower')
        ax.set_xlabel('Time Lag', fontsize=12)
        ax.set_ylabel('Space Lag', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Correlation', fontsize=11)

        plt.tight_layout()

        if self.save_plots:
            try:
                save_path = os.path.join(self.output_dir, filename)
                print(f"Saving plot: {save_path}")
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Successfully saved: {filename}")
            except Exception as e:
                print(f"ERROR saving {filename}: {e}")
                import traceback
                traceback.print_exc()

        plt.close()

    def create_summary_report(self, results: Dict,
                             config: Optional[object] = None,
                             filename: str = 'analysis_summary.txt'):
        """
        Create text summary report of analysis results.

        Parameters
        ----------
        results : Dict
            Analysis results
        config : object, optional
            Configuration object
        filename : str, optional
            Output filename
        """
        report_path = os.path.join(self.output_dir, filename)

        with open(report_path, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("CoastCams Analysis Summary Report\n")
            f.write("=" * 70 + "\n\n")

            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            if config is not None:
                f.write("Configuration Parameters:\n")
                f.write("-" * 70 + "\n")
                f.write(f"  Camera Height: {config.camera_height:.3f} m\n")
                f.write(f"  Pixel Resolution: {config.pixel_resolution:.3f} m/pixel\n")
                f.write(f"  Shoreline Method: {config.shoreline_method}\n\n")

            f.write("Analysis Results:\n")
            f.write("-" * 70 + "\n")

            if 'mean_Hs' in results:
                f.write(f"  Mean Significant Wave Height: {results['mean_Hs']:.3f} m\n")
            if 'mean_Tm' in results:
                f.write(f"  Mean Wave Period: {results['mean_Tm']:.3f} s\n")
            if 'mean_celerity' in results:
                f.write(f"  Mean Wave Celerity: {results['mean_celerity']:.3f} m/s\n")
            if 'mean_depth' in results:
                f.write(f"  Mean Water Depth: {results['mean_depth']:.3f} m\n")

            f.write("\n" + "=" * 70 + "\n")

        print(f"Summary report saved to {report_path}")
