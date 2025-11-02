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

    def plot_matlab_style_summary(self, timestamps: List[datetime],
                                 average_timestack: np.ndarray,
                                 sla_matrix: Optional[np.ndarray],
                                 wave_heights: np.ndarray,
                                 wave_periods: np.ndarray,
                                 water_levels: Optional[np.ndarray] = None,
                                 breakpoint_locations: Optional[np.ndarray] = None,
                                 shoreline_positions: Optional[np.ndarray] = None,
                                 bathymetry: Optional[np.ndarray] = None,
                                 bathymetry_x: Optional[np.ndarray] = None,
                                 rotation: int = 270,
                                 filename: str = 'coastcams_matlab_summary.png'):
        """
        Create MATLAB-style comprehensive summary plot.

        4 subplots:
        1. Average Timestack (grayscale image) with breakpoint and shoreline overlays
        2. Mean Water Level time series (point data)
        3. Significant Wave Height time series
        4. Peak Wave Period time series

        Parameters
        ----------
        timestamps : List[datetime]
            Timestamps for each image
        average_timestack : np.ndarray
            Average intensity profile for each timestack (num_images × spatial_pixels)
        sla_matrix : np.ndarray, optional
            SLA matrix (num_images × spatial_pixels) - plotted in separate figure
        wave_heights : np.ndarray
            Significant wave heights for each timestamp
        wave_periods : np.ndarray
            Peak wave periods for each timestamp
        water_levels : np.ndarray, optional
            Mean sea level / water depth per timestamp (1D array)
        breakpoint_locations : np.ndarray, optional
            Breakpoint locations for each timestamp (in meters)
        shoreline_positions : np.ndarray, optional
            Shoreline positions for each timestamp (in meters)
        bathymetry : np.ndarray, optional
            Bathymetry profile (1D array) - plotted in separate figure
        bathymetry_x : np.ndarray, optional
            Cross-shore distances for bathymetry (1D array)
        rotation : int
            Rotation angle (default: 270)
        filename : str
            Output filename
        """
        print(f"\nCreating MATLAB-style summary plot...")

        # Create figure
        fig = plt.figure(figsize=(12, 10))
        fig.suptitle('CoastCams Analysis Summary', fontsize=16, fontweight='bold')

        # Create 4 subplots
        ax1 = plt.subplot(4, 1, 1)
        ax2 = plt.subplot(4, 1, 2)
        ax3 = plt.subplot(4, 1, 3)
        ax4 = plt.subplot(4, 1, 4)

        # Convert timestamps to matplotlib dates
        time_nums = mdates.date2num(timestamps)

        # Subplot 1: Average Timestack
        # Stack_av shape: (num_timestacks, num_spatial_positions) e.g., (16, 689)
        # MATLAB: rotated_stack = imrotate(uint8(Stack_av), rotation)
        # Then: imagesc(Time_TS, 1:size(rotated_stack,1), rotated_stack)
        # This rotates the (16 x 689) matrix by rotation degrees
        # With rotation=270, (16 x 689) becomes (689 x 16)
        # Then Time_TS (length 16) goes on x-axis, 1:689 on y-axis

        # Apply rotation to match MATLAB (do NOT transpose first!)
        from scipy.ndimage import rotate
        if rotation != 0:
            # Note: MATLAB's imrotate rotates counterclockwise for positive angles
            # scipy's rotate also rotates counterclockwise for positive angles
            rotated_stack = rotate(average_timestack, rotation, reshape=True, order=1)
        else:
            rotated_stack = average_timestack

        # Display with imshow
        # After 270° rotation, (16, 689) becomes (689, 16)
        # extent=[left, right, bottom, top] = [time_start, time_end, 0, num_pixels]
        # With origin='lower': row 0 displayed at bottom, row 688 at top
        # After rotation: original col 0 (beach) → rotated row 688 (top)
        #                 original col 688 (sea) → rotated row 0 (bottom)
        # So: bottom of plot = offshore/sea, top of plot = onshore/beach
        im1 = ax1.imshow(rotated_stack, aspect='auto', cmap='gray',
                        extent=[time_nums[0], time_nums[-1], 0, rotated_stack.shape[0]],
                        origin='lower')

        # Create secondary y-axis for distance in meters
        pixel_resolution = 0.1  # m/pixel
        max_distance = rotated_stack.shape[0] * pixel_resolution
        ax1_right = ax1.twinx()
        # Invert the axis so 0m (beach) is at top, max (offshore) is at bottom
        # This matches the natural interpretation: top = onshore/beach, bottom = offshore/sea
        ax1_right.set_ylim(max_distance, 0)  # Inverted: max to 0
        ax1_right.set_ylabel('Distance from beach [m]', fontsize=12)

        # Overlay breakpoint locations and shoreline positions
        # With inverted axis, we can plot the values directly:
        # - Small values (shoreline ~10m) appear near top (onshore)
        # - Large values (breakpoint ~67m) appear near bottom (offshore)
        if breakpoint_locations is not None and len(breakpoint_locations) > 0:
            valid_bp = ~np.isnan(breakpoint_locations)
            if np.any(valid_bp):
                ax1_right.plot(np.array(time_nums)[valid_bp], breakpoint_locations[valid_bp],
                        'r.', markersize=8, label='Breakpoint Location', alpha=0.8)

        if shoreline_positions is not None and len(shoreline_positions) > 0:
            valid_sl = ~np.isnan(shoreline_positions)
            if np.any(valid_sl):
                ax1_right.plot(np.array(time_nums)[valid_sl], shoreline_positions[valid_sl],
                        'c.', markersize=8, label='Shoreline Position', alpha=0.8)

        # Add legend if overlays were plotted (on right axis)
        if ((breakpoint_locations is not None and np.any(~np.isnan(breakpoint_locations))) or
            (shoreline_positions is not None and np.any(~np.isnan(shoreline_positions)))):
            ax1_right.legend(loc='upper right', fontsize=10, framealpha=0.7)

        ax1.set_title('Average Timestack', fontsize=14, fontweight='bold')
        ax1.set_ylabel('Cross-shore distance [pixels]', fontsize=12)
        # Remove x-axis label (only on bottom subplot)
        ax1.set_xticklabels([])

        # Subplot 2: Mean Water Level time series (point data)
        if water_levels is not None and len(water_levels) > 0:
            ax2.plot(timestamps, water_levels, 'g.-', linewidth=1.5, markersize=6, label='Water Level')
            ax2.set_title('Mean Water Level', fontsize=14, fontweight='bold')
            ax2.set_ylabel('Water Level [m]', fontsize=12)
            ax2.grid(True, alpha=0.3)
        else:
            ax2.text(0.5, 0.5, 'Water level data not available',
                    transform=ax2.transAxes, ha='center', va='center',
                    fontsize=12, color='red')
            ax2.set_title('Mean Water Level', fontsize=14, fontweight='bold')
            ax2.set_ylabel('Water Level [m]', fontsize=12)
        # Remove x-axis label (only on bottom subplot)
        ax2.set_xticklabels([])

        # Subplot 3: Significant Wave Height
        ax3.plot(timestamps, wave_heights, 'r.-', linewidth=1.5, markersize=6)
        ax3.set_title('Significant Wave Height', fontsize=14, fontweight='bold')
        ax3.set_ylabel('$H_s$ [m]', fontsize=12)
        ax3.grid(True, alpha=0.3)
        # Remove x-axis label (only on bottom subplot)
        ax3.set_xticklabels([])

        # Subplot 4: Peak Wave Period
        ax4.plot(timestamps, wave_periods, 'b.-', linewidth=1.5, markersize=6)
        ax4.set_title('Peak Wave Period', fontsize=14, fontweight='bold')
        ax4.set_ylabel('$T_p$ [s]', fontsize=12)
        ax4.set_xlabel('Time', fontsize=12)
        ax4.grid(True, alpha=0.3)
        ax4.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

        # Synchronize x-axes
        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_xlim(time_nums[0], time_nums[-1])

        # Adjust layout
        plt.tight_layout(rect=[0, 0, 1, 0.97])

        # Save figure
        if self.save_plots:
            output_path = os.path.join(self.output_dir, filename)
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"  Saved: {output_path}")

        plt.show()

        # Create separate SLA figure if data is available
        if sla_matrix is not None and not np.all(np.isnan(sla_matrix)):
            self.plot_sla_matrix(timestamps, sla_matrix, filename='sla_matrix.png')

        # Create separate bathymetry figure if data is available
        if bathymetry is not None and bathymetry_x is not None:
            self.plot_bathymetry(bathymetry_x, bathymetry, filename='bathymetry_profile.png')

        return fig

    def plot_sla_matrix(self, timestamps: List[datetime],
                       sla_matrix: np.ndarray,
                       filename: str = 'sla_matrix.png'):
        """
        Create separate figure showing Sea Level Anomaly (SLA) as 2D raster.

        Parameters
        ----------
        timestamps : List[datetime]
            Timestamps for each image
        sla_matrix : np.ndarray
            SLA matrix (num_images × spatial_pixels)
        filename : str
            Output filename
        """
        print(f"\nCreating SLA matrix plot...")

        fig, ax = plt.subplots(figsize=(12, 6))

        # Convert timestamps to matplotlib dates
        time_nums = mdates.date2num(timestamps)
        time_sla = np.linspace(time_nums[0], time_nums[-1], sla_matrix.shape[0])

        # Plot SLA matrix
        im = ax.imshow(sla_matrix.T, aspect='auto', cmap='jet',
                      extent=[time_sla[0], time_sla[-1], 0, sla_matrix.shape[1]],
                      origin='lower')

        ax.set_title('Sea Level Anomaly (SLA)', fontsize=14, fontweight='bold')
        ax.set_xlabel('Time', fontsize=12)
        ax.set_ylabel('Cross-shore Position [pixels]', fontsize=12)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02)
        cbar.set_label('SLA [m]', fontsize=11)

        plt.tight_layout()

        # Save figure
        if self.save_plots:
            output_path = os.path.join(self.output_dir, filename)
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"  Saved: {output_path}")

        plt.show()

        return fig

    def plot_bathymetry(self, x: np.ndarray,
                       bathymetry: np.ndarray,
                       filename: str = 'bathymetry_profile.png'):
        """
        Create separate figure showing bathymetry profile.

        Parameters
        ----------
        x : np.ndarray
            Cross-shore distance array (meters)
        bathymetry : np.ndarray
            Bathymetry depth values (meters)
        filename : str
            Output filename
        """
        print(f"\nCreating bathymetry profile plot...")

        fig, ax = plt.subplots(figsize=(12, 6))

        # Plot bathymetry profile
        ax.plot(x, bathymetry, 'b-', linewidth=2, label='Mean Bathymetry')
        ax.fill_between(x, bathymetry, alpha=0.3)

        ax.set_title('Mean Bathymetry Profile', fontsize=14, fontweight='bold')
        ax.set_xlabel('Cross-shore distance [m]', fontsize=12)
        ax.set_ylabel('Depth [m]', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=11)

        # Invert y-axis to show depth properly (positive down)
        ax.invert_yaxis()

        plt.tight_layout()

        # Save figure
        if self.save_plots:
            output_path = os.path.join(self.output_dir, filename)
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"  Saved: {output_path}")

        plt.show()

        return fig

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
            2D timestack array (time x space) - will be transposed for display
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

        # CRITICAL: Transpose timestack from (time x space) to (space x time)
        # so that time appears on x-axis (imshow shows columns as x-axis)
        # This matches MATLAB's imagesc(Time_TS, 1:size(stack,1), stack)
        timestack_display = timestack.T  # Transpose: (space x time)

        # Plot timestack
        im = ax.imshow(timestack_display, aspect='auto', cmap='viridis', origin='lower')

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
            # CRITICAL: Transpose from (time x space) to (space x time) for correct display
            timestack_display = results['timestack'].T
            im = ax1.imshow(timestack_display, aspect='auto',
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

        # 4. Bathymetry (averaged profile across timestacks)
        if 'depths_smoothed' in results and 'cross_shore_positions' in results:
            ax4 = fig.add_subplot(gs[2, 0])
            ax4.plot(results['cross_shore_positions'], -results['depths_smoothed'],
                    'b-', linewidth=2, marker='o')
            ax4.axhline(y=0, color='cyan', linestyle='--', linewidth=1.5, alpha=0.7)
            ax4.set_title('Averaged Bathymetry Profile', fontsize=12, fontweight='bold')
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
