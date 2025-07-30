"""
Data models for supernova processing pipeline.

This module defines Pydantic models for representing astronomical data such as
spectra and supernova properties, ensuring data consistency and handling
physical units using astropy.
"""
from typing import List, Optional

import astropy.units as u
from astropy.units import Quantity
from pydantic import BaseModel, ConfigDict, field_validator, model_validator


class Spectrum(BaseModel):
    """
    Data model for a single astronomical spectrum.

    Attributes:
        wavelength (astropy.units.Quantity): Wavelength values as a 1D array.
        flux (astropy.units.Quantity): Flux values as a 1D array.
        flux_error (Optional[astropy.units.Quantity]): Optional flux errors.
    """
    # Allow custom types like astropy.units.Quantity which wraps a numpy.ndarray
    model_config = ConfigDict(arbitrary_types_allowed=True)

    wavelength: Quantity
    flux: Quantity
    flux_error: Optional[Quantity] = None

    @model_validator(mode='after')
    def validate_spectrum_data(self) -> 'Spectrum':
        """
        Validates the consistency of spectrum data after model initialization.
        - All data are 1D arrays.
        - Wavelength, flux, and flux_error (if present) have the same shape.
        - Flux and flux_error (if present) have the same units.
        """
        if self.wavelength.ndim != 1 or self.flux.ndim != 1:
            raise ValueError("Wavelength and flux must be 1D arrays.")

        if self.wavelength.shape != self.flux.shape:
            raise ValueError("Wavelength and flux must have the same shape.")

        if self.flux_error is not None:
            if self.flux_error.ndim != 1:
                raise ValueError("Flux error must be a 1D array.")
            if self.flux.shape != self.flux_error.shape:
                raise ValueError("Flux and flux_error must have the same shape.")
            if self.flux.unit != self.flux_error.unit:
                raise ValueError("Flux and flux_error must have the same units.")

        return self


class SupernovaProperties(BaseModel):
    """
    Data model for supernova properties.

    Attributes:
        name (str): The name of the supernova.
        redshift (float): The redshift of the supernova.
        ra (astropy.units.Quantity): Right Ascension in angular units.
        dec (astropy.units.Quantity): Declination in angular units.
        spectra (List[Spectrum]): A list of spectra for the supernova.
    """
    model_config = ConfigDict(arbitrary_types_allowed=True)

    name: str
    redshift: float
    ra: Quantity
    dec: Quantity
    spectra: List[Spectrum] = []

    @field_validator('redshift')
    @classmethod
    def redshift_must_be_non_negative(cls, v: float) -> float:
        """Validates that redshift is a non-negative number."""
        if v < 0:
            raise ValueError("Redshift must be non-negative.")
        return v

    @field_validator('ra', 'dec')
    @classmethod
    def coordinates_must_be_angular(cls, v: Quantity) -> Quantity:
        """Validates that RA and Dec have angular units."""
        if not v.unit.is_equivalent(u.deg):
            raise ValueError("Coordinates must have angular units (e.g., degrees, radians).")
        return v

