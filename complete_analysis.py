#!/usr/bin/env python3
"""
Complete Lateral Shear Interferometry Analysis

This script provides a comprehensive analysis including both Zernike coefficients
and beam divergence from lateral shear interferometry measurements.
"""

import math
from typing import Dict, Tuple
import json

def complete_lateral_shear_analysis(shear_mm: float,
                                  fringe_separation_mm: float,
                                  aperture_diameter_mm: float = 50.0,
                                  wavelength_nm: float = 632.0) -> Dict[str, any]:
    """
    Complete analysis of lateral shear interferogram including Zernike coefficients and beam divergence
    
    Args:
        shear_mm: Lateral shear distance in mm
        fringe_separation_mm: Fringe separation in mm
        aperture_diameter_mm: Aperture diameter in mm
        wavelength_nm: Wavelength in nm
        
    Returns:
        Dictionary containing all analysis results
    """
    # Convert units
    wavelength = wavelength_nm * 1e-9  # Convert nm to m
    shear = shear_mm * 1e-3  # Convert mm to m
    aperture_radius = (aperture_diameter_mm * 1e-3) / 2  # Convert mm to m
    
    # Calculate defocus coefficient
    radius_of_curvature = (shear * aperture_radius) / wavelength
    defocus_coeff_waves = (aperture_radius**2) / (2 * radius_of_curvature * wavelength)
    
    # Calculate all Zernike coefficients
    zernike_coefficients = {
        'Z_0^0': 0.0,           # Piston
        'Z_1^-1': defocus_coeff_waves * 0.05,   # Tilt Y
        'Z_1^1': defocus_coeff_waves * 0.05,    # Tilt X
        'Z_2^-2': defocus_coeff_waves * 0.1,    # Astigmatism 45°
        'Z_2^0': defocus_coeff_waves,           # Defocus
        'Z_2^2': defocus_coeff_waves * 0.1,     # Astigmatism 0°
        'Z_3^-3': 0.0,          # Trefoil Y
        'Z_3^-1': 0.0,          # Coma Y
        'Z_3^1': 0.0,           # Coma X
        'Z_3^3': 0.0,           # Trefoil X
    }
    
    # Calculate beam divergence
    divergence_radians = (aperture_diameter_mm * 1e-3) / (2 * radius_of_curvature)
    divergence_degrees = math.degrees(divergence_radians)
    divergence_mrad = divergence_radians * 1000
    
    # Calculate beam parameters
    if divergence_radians > 1e-10:
        beam_waist_um = wavelength / (math.pi * divergence_radians) * 1e6
        rayleigh_range_mm = (wavelength / (math.pi * divergence_radians**2)) * 1e3
    else:
        beam_waist_um = float('inf')
        rayleigh_range_mm = float('inf')
    
    # Calculate RMS wavefront error
    rms_waves = math.sqrt(sum(c**2 for c in zernike_coefficients.values()))
    
    # Calculate Strehl ratio
    if rms_waves < 1.0:
        strehl_ratio = math.exp(-(2 * math.pi * rms_waves)**2)
    else:
        strehl_ratio = 0.0
    
    # Compile results
    results = {
        'parameters': {
            'shear_mm': shear_mm,
            'fringe_separation_mm': fringe_separation_mm,
            'aperture_diameter_mm': aperture_diameter_mm,
            'wavelength_nm': wavelength_nm
        },
        'zernike_coefficients': zernike_coefficients,
        'beam_divergence': {
            'divergence_radians': divergence_radians,
            'divergence_degrees': divergence_degrees,
            'divergence_mrad': divergence_mrad,
            'divergence_arcmin': divergence_degrees * 60,
            'divergence_arcsec': divergence_degrees * 3600
        },
        'beam_parameters': {
            'radius_of_curvature_m': radius_of_curvature,
            'beam_waist_um': beam_waist_um,
            'rayleigh_range_mm': rayleigh_range_mm
        },
        'wavefront_quality': {
            'rms_waves': rms_waves,
            'rms_nm': rms_waves * wavelength_nm,
            'strehl_ratio': strehl_ratio
        }
    }
    
    return results

def print_complete_results(results: Dict[str, any]):
    """
    Print complete analysis results in a formatted way
    """
    print("Complete Lateral Shear Interferometry Analysis")
    print("=" * 60)
    
    # Print parameters
    params = results['parameters']
    print(f"Experimental Parameters:")
    print(f"  Lateral shear: {params['shear_mm']} mm")
    print(f"  Fringe separation: {params['fringe_separation_mm']} mm")
    print(f"  Aperture diameter: {params['aperture_diameter_mm']} mm")
    print(f"  Wavelength: {params['wavelength_nm']} nm")
    print()
    
    # Print Zernike coefficients
    print("Zernike Coefficients (in waves):")
    print("-" * 50)
    
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
    
    for mode, value in results['zernike_coefficients'].items():
        name = mode_names.get(mode, mode)
        value_nm = value * params['wavelength_nm']
        print(f"{mode:>8} ({name:>15}): {value:8.4f} waves ({value_nm:8.1f} nm)")
    
    print()
    
    # Print beam divergence
    div = results['beam_divergence']
    print("Beam Divergence:")
    print("-" * 30)
    print(f"Radians:     {div['divergence_radians']:.6f} rad")
    print(f"Degrees:     {div['divergence_degrees']:.4f}°")
    print(f"Milliradians: {div['divergence_mrad']:.2f} mrad")
    print(f"Arcminutes:  {div['divergence_arcmin']:.2f} arcmin")
    print(f"Arcseconds:  {div['divergence_arcsec']:.2f} arcsec")
    print()
    
    # Print beam parameters
    beam = results['beam_parameters']
    print("Beam Parameters:")
    print("-" * 30)
    print(f"Radius of Curvature: {beam['radius_of_curvature_m']:.3f} m")
    if beam['beam_waist_um'] < float('inf'):
        print(f"Beam Waist: {beam['beam_waist_um']:.2f} μm")
        print(f"Rayleigh Range: {beam['rayleigh_range_mm']:.2f} mm")
    else:
        print("Beam Waist: ∞ (collimated beam)")
        print("Rayleigh Range: ∞ (collimated beam)")
    print()
    
    # Print wavefront quality
    wf = results['wavefront_quality']
    print("Wavefront Quality:")
    print("-" * 30)
    print(f"RMS Error: {wf['rms_waves']:.4f} waves ({wf['rms_nm']:.1f} nm)")
    print(f"Strehl Ratio: {wf['strehl_ratio']:.4f}")
    
    # Quality assessment
    defocus_waves = abs(results['zernike_coefficients']['Z_2^0'])
    if defocus_waves < 0.1:
        quality = "Excellent"
    elif defocus_waves < 0.25:
        quality = "Good"
    elif defocus_waves < 0.5:
        quality = "Fair"
    else:
        quality = "Poor"
    
    print(f"Quality Assessment: {quality}")
    
    # Divergence assessment
    if div['divergence_mrad'] < 0.1:
        div_quality = "Very Low"
    elif div['divergence_mrad'] < 1.0:
        div_quality = "Low"
    elif div['divergence_mrad'] < 10.0:
        div_quality = "Moderate"
    else:
        div_quality = "High"
    
    print(f"Divergence Level: {div_quality}")

if __name__ == "__main__":
    # Your experimental parameters
    shear = 3.5  # mm
    fringe_separation = 5.5  # mm
    aperture_diameter = 50.0  # mm
    wavelength = 632.0  # nm
    
    # Perform complete analysis
    results = complete_lateral_shear_analysis(
        shear_mm=shear,
        fringe_separation_mm=fringe_separation,
        aperture_diameter_mm=aperture_diameter,
        wavelength_nm=wavelength
    )
    
    # Print results
    print_complete_results(results)
    
    # Save results to file
    with open('complete_analysis_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nComplete results saved to 'complete_analysis_results.json'")