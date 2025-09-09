#!/usr/bin/env python3
"""
Refined Lateral Shear Interferogram Analysis for Zernike Polynomial Extraction

This script provides a more accurate analysis of lateral shear interferograms
based on the actual physics of lateral shear interferometry.

Based on Malacara's "Optical Shop Testing" Chapter 4 - Lateral Shear Interferometers
"""

import math
from typing import Dict, List, Tuple
import json

class RefinedLateralShearAnalyzer:
    """
    Refined analyzer for lateral shear interferograms with better physics modeling
    """
    
    def __init__(self, wavelength: float = 632e-9, aperture_diameter: float = 50e-3):
        """
        Initialize the analyzer with experimental parameters
        """
        self.wavelength = wavelength
        self.aperture_diameter = aperture_diameter
        self.aperture_radius = aperture_diameter / 2
        
        # Zernike polynomial definitions (OSA/ANSI standard)
        self.zernike_modes = [
            (0, 0, "Piston"),
            (1, -1, "Tilt Y"),
            (1, 1, "Tilt X"),
            (2, -2, "Astigmatism 45°"),
            (2, 0, "Defocus"),
            (2, 2, "Astigmatism 0°"),
            (3, -3, "Trefoil Y"),
            (3, -1, "Coma Y"),
            (3, 1, "Coma X"),
            (3, 3, "Trefoil X"),
        ]
    
    def calculate_wavefront_gradient_from_fringes(self, 
                                                fringe_separation_mm: float,
                                                shear_mm: float) -> float:
        """
        Calculate wavefront gradient from fringe separation in lateral shear interferometry
        
        In lateral shear interferometry:
        - The phase difference between sheared wavefronts is: Δφ = (2π/λ) * shear * ∂W/∂x
        - Fringe separation corresponds to 2π phase difference
        - Therefore: ∂W/∂x = λ / shear (for 2π phase difference)
        
        Args:
            fringe_separation_mm: Fringe separation in mm
            shear_mm: Lateral shear in mm
            
        Returns:
            Wavefront gradient in meters/meter
        """
        # Convert to meters
        fringe_separation = fringe_separation_mm * 1e-3
        shear = shear_mm * 1e-3
        
        # For lateral shear, the wavefront gradient is:
        # ∂W/∂x = (λ * fringe_separation) / (2π * shear)
        # But more accurately, for a 2π phase difference:
        wavefront_gradient = self.wavelength / shear
        
        return wavefront_gradient
    
    def estimate_zernike_coefficients_from_gradient(self, 
                                                  wavefront_gradient: float,
                                                  aperture_radius: float) -> Dict[str, float]:
        """
        Estimate Zernike coefficients from wavefront gradient
        
        This uses a more realistic model based on the relationship between
        wavefront gradients and Zernike coefficients in lateral shear interferometry.
        
        Args:
            wavefront_gradient: Wavefront gradient in meters/meter
            aperture_radius: Aperture radius in meters
            
        Returns:
            Dictionary of estimated Zernike coefficients in waves
        """
        coefficients = {}
        
        # For lateral shear interferometry, the most significant contributions are:
        
        # 1. Defocus (Z_2^0) - causes straight, parallel fringes
        # The defocus coefficient is related to the average gradient across the aperture
        # For a circular aperture with defocus, the gradient varies linearly
        defocus_coeff = wavefront_gradient * aperture_radius / (4 * math.pi)
        coefficients['Z_2^0'] = defocus_coeff / self.wavelength
        
        # 2. Astigmatism (Z_2^±2) - causes curved fringes
        # Astigmatism creates hyperbolic fringe patterns
        # The coefficient is typically smaller than defocus
        astigmatism_coeff = wavefront_gradient * aperture_radius / (8 * math.pi)
        coefficients['Z_2^-2'] = astigmatism_coeff / self.wavelength  # Astigmatism 45°
        coefficients['Z_2^2'] = astigmatism_coeff / self.wavelength   # Astigmatism 0°
        
        # 3. Tilt coefficients (Z_1^±1) - cause fringe rotation
        # These are typically small in lateral shear unless there's systematic tilt
        # Estimate based on potential misalignment
        tilt_coeff = wavefront_gradient * aperture_radius / (20 * math.pi)  # Small contribution
        coefficients['Z_1^-1'] = tilt_coeff / self.wavelength  # Tilt Y
        coefficients['Z_1^1'] = tilt_coeff / self.wavelength   # Tilt X
        
        # 4. Higher order terms are typically much smaller
        coefficients['Z_0^0'] = 0.0   # Piston (not measurable in lateral shear)
        coefficients['Z_3^-3'] = 0.0  # Trefoil Y
        coefficients['Z_3^-1'] = 0.0  # Coma Y
        coefficients['Z_3^1'] = 0.0   # Coma X
        coefficients['Z_3^3'] = 0.0   # Trefoil X
        
        return coefficients
    
    def calculate_geometric_optics_estimate(self, 
                                          fringe_separation_mm: float,
                                          shear_mm: float,
                                          aperture_diameter_mm: float) -> Dict[str, float]:
        """
        Calculate Zernike coefficients using geometric optics approximation
        
        This method uses the relationship between fringe patterns and wavefront
        curvature in lateral shear interferometry.
        
        Args:
            fringe_separation_mm: Fringe separation in mm
            shear_mm: Lateral shear in mm
            aperture_diameter_mm: Aperture diameter in mm
            
        Returns:
            Dictionary of Zernike coefficients in waves
        """
        # Convert units
        aperture_radius = (aperture_diameter_mm * 1e-3) / 2  # Convert to meters
        shear = shear_mm * 1e-3  # Convert to meters
        fringe_separation = fringe_separation_mm * 1e-3  # Convert to meters
        
        # Calculate wavefront curvature from fringe separation
        # In lateral shear, the fringe separation is related to the wavefront curvature
        # For a spherical wavefront: fringe_separation = λ * R / (2 * shear)
        # where R is the radius of curvature
        
        # Estimate radius of curvature from fringe separation
        if fringe_separation > 0:
            radius_of_curvature = (2 * shear * fringe_separation) / self.wavelength
        else:
            radius_of_curvature = float('inf')  # Plane wave
        
        # Calculate defocus coefficient from radius of curvature
        # For a spherical wavefront: W(r) = r²/(2R)
        # The defocus coefficient is related to the curvature
        if abs(radius_of_curvature) < float('inf'):
            defocus_coeff = aperture_radius**2 / (2 * radius_of_curvature)
        else:
            defocus_coeff = 0.0
        
        coefficients = {}
        coefficients['Z_2^0'] = defocus_coeff / self.wavelength  # Defocus
        
        # Other coefficients are typically much smaller
        coefficients['Z_0^0'] = 0.0   # Piston
        coefficients['Z_1^-1'] = 0.0  # Tilt Y
        coefficients['Z_1^1'] = 0.0   # Tilt X
        coefficients['Z_2^-2'] = 0.0  # Astigmatism 45°
        coefficients['Z_2^2'] = 0.0   # Astigmatism 0°
        coefficients['Z_3^-3'] = 0.0  # Trefoil Y
        coefficients['Z_3^-1'] = 0.0  # Coma Y
        coefficients['Z_3^1'] = 0.0   # Coma X
        coefficients['Z_3^3'] = 0.0   # Trefoil X
        
        return coefficients
    
    def analyze_interferogram_parameters(self, 
                                       shear_mm: float,
                                       fringe_separation_mm: float,
                                       aperture_diameter_mm: float = 50.0) -> Dict[str, float]:
        """
        Complete analysis based on interferogram parameters using refined methods
        
        Args:
            shear_mm: Lateral shear in mm
            fringe_separation_mm: Fringe separation in mm
            aperture_diameter_mm: Aperture diameter in mm
            
        Returns:
            Dictionary of Zernike coefficients in waves
        """
        # Use geometric optics method for more accurate results
        coefficients = self.calculate_geometric_optics_estimate(
            fringe_separation_mm, shear_mm, aperture_diameter_mm)
        
        return coefficients
    
    def print_results(self, coefficients: Dict[str, float]):
        """
        Print analysis results in a formatted way
        """
        print("Refined Lateral Shear Interferogram Analysis Results")
        print("=" * 60)
        print(f"Wavelength: {self.wavelength * 1e9:.1f} nm")
        print(f"Aperture diameter: {self.aperture_diameter * 1e3:.1f} mm")
        print()
        
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
        
        for mode, value in coefficients.items():
            name = mode_names.get(mode, mode)
            value_nm = value * self.wavelength * 1e9
            print(f"{mode:>8} ({name:>15}): {value:8.4f} waves ({value_nm:8.1f} nm)")
        
        # Calculate RMS wavefront error
        rms_waves = math.sqrt(sum(c**2 for c in coefficients.values()))
        rms_nm = rms_waves * self.wavelength * 1e9
        print("-" * 50)
        print(f"{'RMS':>8} ({'Wavefront Error':>15}): {rms_waves:8.4f} waves ({rms_nm:8.1f} nm)")
        
        # Calculate Strehl ratio (approximation)
        if rms_waves < 1.0:  # Only calculate if reasonable
            strehl_ratio = math.exp(-(2 * math.pi * rms_waves)**2)
            print(f"{'Strehl':>8} ({'Ratio':>15}): {strehl_ratio:8.4f}")
        else:
            print(f"{'Strehl':>8} ({'Ratio':>15}): < 0.001")
        
        print()
        print("Analysis Notes:")
        print("- This analysis uses geometric optics approximation")
        print("- Results are based on fringe separation and lateral shear")
        print("- Defocus is the primary aberration detected")
        print("- Higher-order aberrations require full image analysis")
        
        # Provide interpretation
        defocus_waves = abs(coefficients.get('Z_2^0', 0))
        if defocus_waves < 0.1:
            print("- Wavefront quality: Excellent (defocus < 0.1 waves)")
        elif defocus_waves < 0.25:
            print("- Wavefront quality: Good (defocus < 0.25 waves)")
        elif defocus_waves < 0.5:
            print("- Wavefront quality: Fair (defocus < 0.5 waves)")
        else:
            print("- Wavefront quality: Poor (defocus > 0.5 waves)")


def analyze_lateral_shear_interferogram(shear_mm: float,
                                      fringe_separation_mm: float,
                                      aperture_diameter_mm: float = 50.0,
                                      wavelength_nm: float = 632.0) -> Dict[str, float]:
    """
    Main function to analyze lateral shear interferogram parameters and return Zernike coefficients
    
    Args:
        shear_mm: Lateral shear distance in mm
        fringe_separation_mm: Fringe separation in mm
        aperture_diameter_mm: Diameter of the aperture in mm
        wavelength_nm: Wavelength of the laser in nm
        
    Returns:
        Dictionary containing Zernike coefficients in waves
    """
    # Convert units
    wavelength = wavelength_nm * 1e-9  # Convert nm to m
    aperture_diameter = aperture_diameter_mm * 1e-3  # Convert mm to m
    
    # Create analyzer
    analyzer = RefinedLateralShearAnalyzer(wavelength=wavelength, 
                                         aperture_diameter=aperture_diameter)
    
    # Analyze interferogram
    coefficients = analyzer.analyze_interferogram_parameters(
        shear_mm, fringe_separation_mm, aperture_diameter_mm)
    
    # Print results
    analyzer.print_results(coefficients)
    
    return coefficients


if __name__ == "__main__":
    # Your experimental parameters
    aperture_diameter = 50.0  # mm
    wavelength = 632.0  # nm
    shear = 3.5  # mm
    fringe_separation = 5.5  # mm
    
    print("Analyzing lateral shear interferogram with refined method...")
    print(f"Parameters:")
    print(f"  Aperture diameter: {aperture_diameter} mm")
    print(f"  Wavelength: {wavelength} nm")
    print(f"  Lateral shear: {shear} mm")
    print(f"  Fringe separation: {fringe_separation} mm")
    print()
    
    try:
        coefficients = analyze_lateral_shear_interferogram(
            shear_mm=shear,
            fringe_separation_mm=fringe_separation,
            aperture_diameter_mm=aperture_diameter,
            wavelength_nm=wavelength
        )
        
        # Save results to file
        results = {
            'parameters': {
                'aperture_diameter_mm': aperture_diameter,
                'wavelength_nm': wavelength,
                'shear_mm': shear,
                'fringe_separation_mm': fringe_separation
            },
            'zernike_coefficients': coefficients,
            'rms_waves': math.sqrt(sum(c**2 for c in coefficients.values())),
            'rms_nm': math.sqrt(sum(c**2 for c in coefficients.values())) * wavelength
        }
        
        with open('refined_lateral_shear_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\nResults saved to 'refined_lateral_shear_results.json'")
        
    except Exception as e:
        print(f"Error analyzing interferogram: {e}")
        import traceback
        traceback.print_exc()