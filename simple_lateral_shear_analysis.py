#!/usr/bin/env python3
"""
Simplified Lateral Shear Interferogram Analysis for Zernike Polynomial Extraction

This script provides a basic analysis of lateral shear interferograms to extract 
wavefront aberrations in terms of the first 10 Zernike polynomials.

Based on Malacara's "Optical Shop Testing" Chapter 4 - Lateral Shear Interferometers
"""

import math
from typing import Dict, List, Tuple
import json

class SimpleLateralShearAnalyzer:
    """
    Simplified analyzer for lateral shear interferograms
    """
    
    def __init__(self, wavelength: float = 632e-9, aperture_diameter: float = 50e-3):
        """
        Initialize the analyzer with experimental parameters
        
        Args:
            wavelength: Laser wavelength in meters (default: 632 nm)
            aperture_diameter: Aperture diameter in meters (default: 50 mm)
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
    
    def factorial(self, n: int) -> int:
        """Calculate factorial"""
        if n <= 1:
            return 1
        return n * self.factorial(n - 1)
    
    def zernike_radial_polynomial(self, n: int, m: int, rho: float) -> float:
        """
        Calculate Zernike radial polynomial R_n^m
        
        Args:
            n: Radial order
            m: Azimuthal frequency
            rho: Normalized radial coordinate (0 to 1)
            
        Returns:
            Radial polynomial value
        """
        R = 0.0
        for k in range((n - abs(m)) // 2 + 1):
            try:
                coeff = ((-1)**k * self.factorial(n - k)) / (
                    self.factorial(k) * 
                    self.factorial((n + abs(m)) // 2 - k) * 
                    self.factorial((n - abs(m)) // 2 - k)
                )
                R += coeff * (rho**(n - 2*k))
            except:
                continue
        return R
    
    def zernike_polynomial(self, n: int, m: int, rho: float, theta: float) -> float:
        """
        Calculate Zernike polynomial Z_n^m
        
        Args:
            n: Radial order
            m: Azimuthal frequency
            rho: Normalized radial coordinate (0 to 1)
            theta: Azimuthal coordinate (0 to 2π)
            
        Returns:
            Zernike polynomial value
        """
        # Normalization factor
        if m == 0:
            N = math.sqrt(n + 1)
        else:
            N = math.sqrt(2 * (n + 1))
        
        # Radial polynomial
        R = self.zernike_radial_polynomial(n, m, rho)
        
        # Angular part
        if m >= 0:
            angular = math.cos(m * theta)
        else:
            angular = math.sin(abs(m) * theta)
        
        return N * R * angular
    
    def estimate_fringe_phase_from_separation(self, fringe_separation_mm: float, 
                                            shear_mm: float) -> float:
        """
        Estimate phase from fringe separation in lateral shear interferometry
        
        In lateral shear interferometry, the phase difference is related to the
        wavefront gradient by: Δφ = (2π/λ) * shear * ∂W/∂x
        
        Args:
            fringe_separation_mm: Fringe separation in mm
            shear_mm: Lateral shear in mm
            
        Returns:
            Estimated phase gradient in radians/mm
        """
        # Convert to meters
        fringe_separation = fringe_separation_mm * 1e-3
        shear = shear_mm * 1e-3
        
        # Phase difference between fringes is 2π
        phase_gradient = 2 * math.pi / fringe_separation
        
        return phase_gradient
    
    def calculate_zernike_coefficients_from_parameters(self, 
                                                     shear_mm: float,
                                                     fringe_separation_mm: float,
                                                     aperture_diameter_mm: float) -> Dict[str, float]:
        """
        Calculate Zernike coefficients from experimental parameters
        
        This is a simplified approach that estimates coefficients based on
        the fringe pattern characteristics.
        
        Args:
            shear_mm: Lateral shear in mm
            fringe_separation_mm: Fringe separation in mm
            aperture_diameter_mm: Aperture diameter in mm
            
        Returns:
            Dictionary of estimated Zernike coefficients in waves
        """
        # Convert units
        shear = shear_mm * 1e-3  # Convert to meters
        fringe_separation = fringe_separation_mm * 1e-3  # Convert to meters
        aperture_radius = (aperture_diameter_mm * 1e-3) / 2  # Convert to meters
        
        # Estimate phase gradient from fringe separation
        phase_gradient = self.estimate_fringe_phase_from_separation(
            fringe_separation_mm, shear_mm)
        
        # Convert to wavefront gradient (in meters/meter)
        wavefront_gradient = phase_gradient * self.wavelength / (2 * math.pi * shear)
        
        # Estimate coefficients based on typical lateral shear patterns
        coefficients = {}
        
        # For a circular aperture with lateral shear, we can estimate some coefficients
        # This is a simplified model - in practice, you'd need full image analysis
        
        # Estimate defocus (Z_2^0) - most common in lateral shear
        # Defocus creates straight, parallel fringes
        if fringe_separation_mm > 0:
            # Rough estimation: defocus coefficient proportional to fringe density
            defocus_estimate = wavefront_gradient * aperture_radius / (2 * math.pi)
            coefficients['Z_2^0'] = defocus_estimate / self.wavelength
        
        # Estimate tilt coefficients (Z_1^±1) - cause fringe rotation
        # These are typically small in lateral shear unless there's systematic tilt
        coefficients['Z_1^-1'] = 0.0  # Tilt Y
        coefficients['Z_1^1'] = 0.0   # Tilt X
        
        # Estimate astigmatism (Z_2^±2) - causes curved fringes
        # Astigmatism creates hyperbolic fringe patterns
        astigmatism_estimate = wavefront_gradient * aperture_radius / (4 * math.pi)
        coefficients['Z_2^-2'] = astigmatism_estimate / self.wavelength  # Astigmatism 45°
        coefficients['Z_2^2'] = astigmatism_estimate / self.wavelength   # Astigmatism 0°
        
        # Higher order terms are typically smaller
        coefficients['Z_0^0'] = 0.0   # Piston
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
        Complete analysis based on interferogram parameters
        
        Args:
            shear_mm: Lateral shear in mm
            fringe_separation_mm: Fringe separation in mm
            aperture_diameter_mm: Aperture diameter in mm
            
        Returns:
            Dictionary of Zernike coefficients in waves
        """
        coefficients = self.calculate_zernike_coefficients_from_parameters(
            shear_mm, fringe_separation_mm, aperture_diameter_mm)
        
        return coefficients
    
    def print_results(self, coefficients: Dict[str, float]):
        """
        Print analysis results in a formatted way
        """
        print("Lateral Shear Interferogram Analysis Results")
        print("=" * 50)
        print(f"Wavelength: {self.wavelength * 1e9:.1f} nm")
        print(f"Aperture diameter: {self.aperture_diameter * 1e3:.1f} mm")
        print()
        
        print("Zernike Coefficients (in waves):")
        print("-" * 40)
        
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
        print("-" * 40)
        print(f"{'RMS':>8} ({'Wavefront Error':>15}): {rms_waves:8.4f} waves ({rms_nm:8.1f} nm)")
        
        # Calculate Strehl ratio (approximation)
        strehl_ratio = math.exp(-(2 * math.pi * rms_waves)**2)
        print(f"{'Strehl':>8} ({'Ratio':>15}): {strehl_ratio:8.4f}")
        
        print()
        print("Analysis Notes:")
        print("- This is a simplified analysis based on fringe separation")
        print("- For accurate results, full image analysis is recommended")
        print("- The coefficients are estimates based on typical lateral shear patterns")
        print("- Higher-order aberrations may not be accurately captured")


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
    analyzer = SimpleLateralShearAnalyzer(wavelength=wavelength, 
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
    
    print("Analyzing lateral shear interferogram...")
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
        
        with open('lateral_shear_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\nResults saved to 'lateral_shear_results.json'")
        
    except Exception as e:
        print(f"Error analyzing interferogram: {e}")
        import traceback
        traceback.print_exc()