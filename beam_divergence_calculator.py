#!/usr/bin/env python3
"""
Beam Divergence Calculator for Lateral Shear Interferometry

This script calculates the far field beam divergence from Zernike defocus coefficients
obtained from lateral shear interferometry measurements.

Based on the relationship between defocus and beam divergence in optical systems.
"""

import math
from typing import Dict, Tuple

def calculate_beam_divergence_from_defocus(defocus_coeff_waves: float,
                                         wavelength_nm: float,
                                         aperture_diameter_mm: float) -> Dict[str, float]:
    """
    Calculate far field beam divergence from Zernike defocus coefficient
    
    The relationship between defocus and beam divergence is:
    - Defocus coefficient: C_2^0 = (aperture_radius)²/(2R*λ)
    - Beam divergence: θ = aperture_diameter/(2*R)
    - Therefore: θ = (C_2^0 * λ * aperture_diameter)/(aperture_radius)²
    
    Args:
        defocus_coeff_waves: Zernike defocus coefficient in waves
        wavelength_nm: Wavelength in nanometers
        aperture_diameter_mm: Aperture diameter in millimeters
        
    Returns:
        Dictionary containing divergence angles in various units
    """
    # Convert units
    wavelength = wavelength_nm * 1e-9  # Convert nm to m
    aperture_diameter = aperture_diameter_mm * 1e-3  # Convert mm to m
    aperture_radius = aperture_diameter / 2
    
    # Convert defocus coefficient from waves to meters
    defocus_meters = defocus_coeff_waves * wavelength
    
    # Calculate radius of curvature from defocus coefficient
    # C_2^0 = (aperture_radius)²/(2R*λ)
    # Therefore: R = (aperture_radius)²/(2*C_2^0*λ)
    if abs(defocus_meters) > 1e-15:  # Avoid division by zero
        radius_of_curvature = (aperture_radius**2) / (2 * defocus_meters)
    else:
        radius_of_curvature = float('inf')
    
    # Calculate beam divergence
    # For a spherical wavefront: θ = aperture_diameter/(2*R)
    if abs(radius_of_curvature) < float('inf'):
        divergence_radians = aperture_diameter / (2 * radius_of_curvature)
    else:
        divergence_radians = 0.0
    
    # Convert to different units
    divergence_degrees = math.degrees(divergence_radians)
    divergence_mrad = divergence_radians * 1000  # milliradians
    divergence_arcmin = divergence_degrees * 60  # arcminutes
    divergence_arcsec = divergence_degrees * 3600  # arcseconds
    
    # Calculate beam waist and Rayleigh range (for Gaussian beam approximation)
    # For a Gaussian beam with divergence θ: w_0 = λ/(π*θ)
    if divergence_radians > 1e-10:
        beam_waist_um = wavelength / (math.pi * divergence_radians) * 1e6  # Convert to micrometers
        rayleigh_range_mm = (wavelength / (math.pi * divergence_radians**2)) * 1e3  # Convert to mm
    else:
        beam_waist_um = float('inf')
        rayleigh_range_mm = float('inf')
    
    return {
        'divergence_radians': divergence_radians,
        'divergence_degrees': divergence_degrees,
        'divergence_mrad': divergence_mrad,
        'divergence_arcmin': divergence_arcmin,
        'divergence_arcsec': divergence_arcsec,
        'radius_of_curvature_m': radius_of_curvature,
        'beam_waist_um': beam_waist_um,
        'rayleigh_range_mm': rayleigh_range_mm
    }

def calculate_beam_divergence_from_parameters(shear_mm: float,
                                            fringe_separation_mm: float,
                                            aperture_diameter_mm: float = 50.0,
                                            wavelength_nm: float = 632.0) -> Dict[str, float]:
    """
    Calculate beam divergence directly from lateral shear interferometry parameters
    
    This function combines the defocus calculation with divergence calculation
    in one step.
    
    Args:
        shear_mm: Lateral shear distance in mm
        fringe_separation_mm: Fringe separation in mm
        aperture_diameter_mm: Aperture diameter in mm
        wavelength_nm: Wavelength in nm
        
    Returns:
        Dictionary containing both defocus and divergence information
    """
    # Convert units
    wavelength = wavelength_nm * 1e-9  # Convert nm to m
    shear = shear_mm * 1e-3  # Convert mm to m
    aperture_radius = (aperture_diameter_mm * 1e-3) / 2  # Convert mm to m
    
    # Calculate defocus coefficient (same as in the main analysis)
    radius_of_curvature = (shear * aperture_radius) / wavelength
    defocus_coeff_waves = (aperture_radius**2) / (2 * radius_of_curvature * wavelength)
    
    # Calculate divergence from the defocus coefficient
    divergence_info = calculate_beam_divergence_from_defocus(
        defocus_coeff_waves, wavelength_nm, aperture_diameter_mm)
    
    # Add defocus information to the results
    divergence_info['defocus_coeff_waves'] = defocus_coeff_waves
    divergence_info['defocus_nm'] = defocus_coeff_waves * wavelength_nm
    
    return divergence_info

def print_beam_divergence_results(results: Dict[str, float]):
    """
    Print beam divergence results in a formatted way
    """
    print("Beam Divergence Analysis Results")
    print("=" * 50)
    
    print(f"Defocus Coefficient: {results['defocus_coeff_waves']:.4f} waves ({results['defocus_nm']:.1f} nm)")
    print(f"Radius of Curvature: {results['radius_of_curvature_m']:.3f} m")
    print()
    
    print("Beam Divergence Angles:")
    print("-" * 30)
    print(f"Radians:     {results['divergence_radians']:.6f} rad")
    print(f"Degrees:     {results['divergence_degrees']:.4f}°")
    print(f"Milliradians: {results['divergence_mrad']:.2f} mrad")
    print(f"Arcminutes:  {results['divergence_arcmin']:.2f} arcmin")
    print(f"Arcseconds:  {results['divergence_arcsec']:.2f} arcsec")
    print()
    
    print("Beam Parameters (Gaussian approximation):")
    print("-" * 40)
    if results['beam_waist_um'] < float('inf'):
        print(f"Beam Waist:   {results['beam_waist_um']:.2f} μm")
        print(f"Rayleigh Range: {results['rayleigh_range_mm']:.2f} mm")
    else:
        print("Beam Waist:   ∞ (collimated beam)")
        print("Rayleigh Range: ∞ (collimated beam)")
    
    print()
    print("Interpretation:")
    if results['divergence_mrad'] < 0.1:
        print("- Very low divergence: Excellent beam quality")
    elif results['divergence_mrad'] < 1.0:
        print("- Low divergence: Good beam quality")
    elif results['divergence_mrad'] < 10.0:
        print("- Moderate divergence: Fair beam quality")
    else:
        print("- High divergence: Poor beam quality")

if __name__ == "__main__":
    # Your experimental parameters
    shear = 3.5  # mm
    fringe_separation = 5.5  # mm
    aperture_diameter = 50.0  # mm
    wavelength = 632.0  # nm
    
    print("Beam Divergence Analysis from Lateral Shear Interferometry")
    print("=" * 60)
    print(f"Parameters:")
    print(f"  Lateral shear: {shear} mm")
    print(f"  Fringe separation: {fringe_separation} mm")
    print(f"  Aperture diameter: {aperture_diameter} mm")
    print(f"  Wavelength: {wavelength} nm")
    print()
    
    # Calculate beam divergence
    results = calculate_beam_divergence_from_parameters(
        shear_mm=shear,
        fringe_separation_mm=fringe_separation,
        aperture_diameter_mm=aperture_diameter,
        wavelength_nm=wavelength
    )
    
    # Print results
    print_beam_divergence_results(results)
    
    print()
    print("Note: This analysis assumes the beam has a spherical wavefront.")
    print("For non-spherical wavefronts, the divergence may vary across the beam.")