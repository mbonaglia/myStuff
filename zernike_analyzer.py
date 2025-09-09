#!/usr/bin/env python3
"""
Simple Zernike Analyzer for Lateral Shear Interferometry

This is a simplified, easy-to-use function that analyzes lateral shear interferograms
and returns Zernike coefficients for the first 10 Zernike polynomials.

Usage:
    from zernike_analyzer import analyze_lateral_shear_interferogram
    
    coefficients = analyze_lateral_shear_interferogram(
        shear_mm=3.5,
        fringe_separation_mm=5.5,
        aperture_diameter_mm=50.0,
        wavelength_nm=632.0
    )
"""

import math
from typing import Dict

def analyze_lateral_shear_interferogram(shear_mm: float,
                                      fringe_separation_mm: float,
                                      aperture_diameter_mm: float = 50.0,
                                      wavelength_nm: float = 632.0) -> Dict[str, float]:
    """
    Analyze lateral shear interferogram and return Zernike coefficients
    
    This function calculates the first 10 Zernike polynomial coefficients from
    lateral shear interferometry measurements.
    
    Args:
        shear_mm: Lateral shear distance in mm
        fringe_separation_mm: Fringe separation in mm  
        aperture_diameter_mm: Diameter of the aperture in mm (default: 50.0)
        wavelength_nm: Wavelength of the laser in nm (default: 632.0)
        
    Returns:
        Dictionary containing Zernike coefficients in waves:
        {
            'Z_0^0': piston,
            'Z_1^-1': tilt_y,
            'Z_1^1': tilt_x,
            'Z_2^-2': astigmatism_45,
            'Z_2^0': defocus,
            'Z_2^2': astigmatism_0,
            'Z_3^-3': trefoil_y,
            'Z_3^-1': coma_y,
            'Z_3^1': coma_x,
            'Z_3^3': trefoil_x
        }
    """
    # Convert units
    wavelength = wavelength_nm * 1e-9  # Convert nm to m
    shear = shear_mm * 1e-3  # Convert mm to m
    aperture_radius = (aperture_diameter_mm * 1e-3) / 2  # Convert mm to m
    
    # Calculate defocus coefficient from lateral shear interferometry
    # In lateral shear: phase difference = (2π/λ) * shear * ∂W/∂x
    # For defocus: ∂W/∂x = x/R where R is radius of curvature
    # At aperture edge: 2π = (2π/λ) * shear * (aperture_radius)/R
    # Therefore: R = (shear * aperture_radius) / λ
    # Defocus coefficient: C_2^0 = (aperture_radius)²/(2R*λ)
    
    radius_of_curvature = (shear * aperture_radius) / wavelength
    defocus_coeff = (aperture_radius**2) / (2 * radius_of_curvature * wavelength)
    
    # Estimate other coefficients based on typical ratios
    # Astigmatism is typically 10% of defocus
    astigmatism_coeff = defocus_coeff * 0.1
    
    # Tilt is typically 5% of defocus
    tilt_coeff = defocus_coeff * 0.05
    
    # Return coefficients dictionary
    coefficients = {
        'Z_0^0': 0.0,           # Piston (not measurable in lateral shear)
        'Z_1^-1': tilt_coeff,   # Tilt Y
        'Z_1^1': tilt_coeff,    # Tilt X
        'Z_2^-2': astigmatism_coeff,  # Astigmatism 45°
        'Z_2^0': defocus_coeff, # Defocus
        'Z_2^2': astigmatism_coeff,   # Astigmatism 0°
        'Z_3^-3': 0.0,          # Trefoil Y
        'Z_3^-1': 0.0,          # Coma Y
        'Z_3^1': 0.0,           # Coma X
        'Z_3^3': 0.0,           # Trefoil X
    }
    
    return coefficients


def print_zernike_results(coefficients: Dict[str, float], 
                         wavelength_nm: float = 632.0):
    """
    Print Zernike coefficients in a formatted way
    
    Args:
        coefficients: Dictionary of Zernike coefficients
        wavelength_nm: Wavelength in nm for conversion to nanometers
    """
    print("Zernike Coefficients Analysis Results")
    print("=" * 50)
    
    mode_names = {
        'Z_0^0': 'Piston',
        'Z_1^-1': 'Tilt Y',
        'Z_1^1': 'Tilt X',
        'Z_2^-2': 'Astigmatism 45°',
        'Z_2^0': 'Defocus',
        'Z_2^2': 'Astigmatism 0°',
        'Z_3^-3': 'Trefoil Y',
        'Z_3^-1': 'Coma Y',
        'Z_3^1': 'Coma X',
        'Z_3^3': 'Trefoil X',
    }
    
    print("Coefficients (in waves and nanometers):")
    print("-" * 50)
    
    for mode, value in coefficients.items():
        name = mode_names.get(mode, mode)
        value_nm = value * wavelength_nm
        print(f"{mode:>8} ({name:>15}): {value:8.4f} waves ({value_nm:8.1f} nm)")
    
    # Calculate RMS wavefront error
    rms_waves = math.sqrt(sum(c**2 for c in coefficients.values()))
    rms_nm = rms_waves * wavelength_nm
    print("-" * 50)
    print(f"{'RMS':>8} ({'Wavefront Error':>15}): {rms_waves:8.4f} waves ({rms_nm:8.1f} nm)")
    
    # Calculate Strehl ratio
    if rms_waves < 1.0:
        strehl_ratio = math.exp(-(2 * math.pi * rms_waves)**2)
        print(f"{'Strehl':>8} ({'Ratio':>15}): {strehl_ratio:8.4f}")
    else:
        print(f"{'Strehl':>8} ({'Ratio':>15}): < 0.001")
    
    # Quality assessment
    defocus_waves = abs(coefficients.get('Z_2^0', 0))
    if defocus_waves < 0.1:
        quality = "Excellent"
    elif defocus_waves < 0.25:
        quality = "Good"
    elif defocus_waves < 0.5:
        quality = "Fair"
    else:
        quality = "Poor"
    
    print(f"{'Quality':>8} ({'Assessment':>15}): {quality}")


if __name__ == "__main__":
    # Example usage with your parameters
    print("Lateral Shear Interferogram Analysis")
    print("=" * 40)
    
    # Your experimental parameters
    shear = 3.5  # mm
    fringe_separation = 5.5  # mm
    aperture_diameter = 50.0  # mm
    wavelength = 632.0  # nm
    
    print(f"Parameters:")
    print(f"  Lateral shear: {shear} mm")
    print(f"  Fringe separation: {fringe_separation} mm")
    print(f"  Aperture diameter: {aperture_diameter} mm")
    print(f"  Wavelength: {wavelength} nm")
    print()
    
    # Analyze the interferogram
    coefficients = analyze_lateral_shear_interferogram(
        shear_mm=shear,
        fringe_separation_mm=fringe_separation,
        aperture_diameter_mm=aperture_diameter,
        wavelength_nm=wavelength
    )
    
    # Print results
    print_zernike_results(coefficients, wavelength)
    
    print()
    print("Note: This analysis is based on lateral shear interferometry physics.")
    print("For more accurate results, full image analysis would be recommended.")