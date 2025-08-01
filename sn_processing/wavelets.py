"""
Module for wavelet analysis of supernova spectra.

This module contains functions for performing the à trous wavelet transform
and extracting spectral features from the wavelet coefficients.
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Dict, Tuple
from scipy.ndimage import convolve1d


def atrous_transform(signal: np.ndarray, num_scales: int) -> np.ndarray:
    """
    Performs a 1D à trous wavelet transform on a signal.

    Args:
        signal: The 1D input signal (e.g., the flux array).
        num_scales: The number of wavelet scales to compute.

    Returns:
        A 2D numpy array where rows are the wavelet coefficients for each
        scale, plus the final smoothed residual. Shape is (num_scales + 1, len(signal)).
    """
    # The B3 spline scaling function is used as the convolution kernel.
    kernel = np.array([1/16, 1/4, 3/8, 1/4, 1/16])
    wavelet_coeffs = []
    c_prev = signal.copy()

    for j in range(num_scales):
        # At each scale 'j', the kernel is upsampled by inserting 2**j - 1 zeros
        # between its elements. This is equivalent to dilating the kernel.
        stride = 2**j
        upsampled_kernel = np.zeros((len(kernel) - 1) * stride + 1)
        upsampled_kernel[::stride] = kernel

        # Convolve with the upsampled kernel. 'mirror' mode handles boundaries.
        c_next = convolve1d(c_prev, upsampled_kernel, mode='mirror')

        # The wavelet coefficients for this scale are the difference plane.
        w_j = c_prev - c_next
        wavelet_coeffs.append(w_j)

        c_prev = c_next

    # The final smoothed (residual) plane is also part of the output.
    wavelet_coeffs.append(c_prev)
    return np.array(wavelet_coeffs)


def sum_wavelet_scales(
    wavelet_coeffs: np.ndarray,
    scales_to_sum: List[int]
) -> np.ndarray:
    """
    Sums specific wavelet scales from the coefficient array.

    Args:
        wavelet_coeffs: The 2D array of wavelet coefficients from `atrous_transform`.
        scales_to_sum: A list of 1-based indices of the scales to sum
                       (e.g., [2, 3, 4]).

    Returns:
        A 1D numpy array representing the sum of the selected wavelet scales.
    """
    # Convert 1-based user-facing scales to 0-based numpy indices
    indices_to_sum = [s - 1 for s in scales_to_sum]

    # Validate indices
    num_scales = wavelet_coeffs.shape[0] - 1
    for i in indices_to_sum:
        if not 0 <= i < num_scales:
            raise ValueError(f"Invalid scale index {i+1}. Available scales are 1 to {num_scales}.")

    return np.sum(wavelet_coeffs[indices_to_sum], axis=0)

def calculate_chi_index(
    wavelet_sum: np.ndarray,
    spectral_axis: np.ndarray,
    wl_min: float,
    wl_max: float
) -> float:
    """
    Calculates the 'chi' spectral index for a feature in a given wavelength range.
    The definition is based on the notes for sn1a_template_20110517.pro.
    """
    # Find indices corresponding to the wavelength region
    feature_indices = np.where((spectral_axis >= wl_min) & (spectral_axis <= wl_max))[0]

    if len(feature_indices) == 0:
        return 0.0

    # Extract the wavelet sum in that region
    feature_wavelet_sum = wavelet_sum[feature_indices]

    # Calculate standard deviation in the region
    std_dev = np.std(feature_wavelet_sum)

    if std_dev == 0:
        return 0.0

    # The index is the sum of the absolute values, normalized by std dev.
    chi_index = np.sum(np.abs(feature_wavelet_sum)) / std_dev
    return chi_index


def extract_features(
    wavelet_sum: np.ndarray,
    spectral_axis: np.ndarray,
    feature_definitions: Dict[str, Tuple[float, float]]
) -> Dict[str, float]:
    """
    Calculates the chi index for multiple predefined spectral features.

    Args:
        wavelet_sum: The 1D array of the summed wavelet scales.
        spectral_axis: The corresponding 1D wavelength array.
        feature_definitions: A dictionary mapping feature names to their
                             (min_wl, max_wl) tuples.

    Returns:
        A dictionary mapping feature names to their calculated chi index.
    """
    chi_indices = {}
    for feature_name, (wl_min, wl_max) in feature_definitions.items():
        chi = calculate_chi_index(wavelet_sum, spectral_axis, wl_min, wl_max)
        chi_indices[feature_name] = chi
    return chi_indices


def plot_wavelet_decomposition(
    spectral_axis: np.ndarray,
    original_signal: np.ndarray,
    wavelet_coeffs: np.ndarray,
    output_path: Path,
    title: str = "Wavelet Decomposition"
):
    """
    Visualizes the à trous wavelet decomposition and saves it to a file.

    Creates a multi-panel plot showing the original signal, each wavelet scale,
    and the final smoothed residual. The plot is saved to disk and not shown
    interactively.

    Args:
        spectral_axis: The 1D wavelength array for the x-axis.
        original_signal: The original 1D flux array.
        wavelet_coeffs: The 2D array of wavelet coefficients from `atrous_transform`.
        output_path: The path where the output plot image will be saved.
        title: The main title for the plot.
    """
    num_scales = wavelet_coeffs.shape[0] - 1
    num_plots = num_scales + 2  # Original + all scales + residual

    fig, axes = plt.subplots(
        num_plots, 1, figsize=(12, 2.5 * num_plots), sharex=True, constrained_layout=True
    )
    fig.suptitle(title, fontsize=16)

    # Plot Original Signal
    axes[0].plot(spectral_axis, original_signal, color='black', label='Original Signal')
    axes[0].set_title("Original Signal")
    axes[0].set_ylabel("Flux")
    axes[0].legend(loc='upper right')

    # Plot Wavelet Scales
    for i in range(num_scales):
        ax = axes[i + 1]
        ax.plot(spectral_axis, wavelet_coeffs[i, :], color=f'C{i}', label=f'Scale {i+1}')
        ax.set_title(f"Wavelet Scale {i + 1}")
        ax.axhline(0, color='grey', linestyle='--', linewidth=0.8)
        ax.set_ylabel(f"$w_{i+1}$")
        ax.legend(loc='upper right')

    # Plot Smoothed Residual
    ax = axes[num_scales + 1]
    ax.plot(spectral_axis, wavelet_coeffs[-1, :], color='gray', label='Residual')
    ax.set_title("Smoothed Residual")
    ax.set_ylabel("$c_{res}$")
    ax.set_xlabel("Wavelength (Å)")
    ax.legend(loc='upper right')

    # Ensure output directory exists and save the figure
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150)
    plt.close(fig)  # Close the figure to free memory and prevent interactive display
    print(f"    - Saved wavelet decomposition plot to {output_path}")
