"""
Module for wavelet analysis of supernova spectra.

This module contains functions for performing the à trous wavelet transform
and extracting spectral features from the wavelet coefficients.
"""
import numpy as np
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

