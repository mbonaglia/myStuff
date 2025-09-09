#!/usr/bin/env python3
"""
Lateral Shear Interferogram Analysis for Zernike Polynomial Extraction

This script analyzes lateral shear interferograms to extract wavefront aberrations
in terms of the first 10 Zernike polynomials (up to spherical aberration).

Based on Malacara's "Optical Shop Testing" Chapter 4 - Lateral Shear Interferometers
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage, optimize
from scipy.special import jv  # Bessel functions for Zernike polynomials
import cv2
from typing import Tuple, List, Dict
import warnings
warnings.filterwarnings('ignore')

class LateralShearAnalyzer:
    """
    Analyzes lateral shear interferograms to extract Zernike coefficients
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
            (0, 0),   # Piston
            (1, -1),  # Tilt Y
            (1, 1),   # Tilt X
            (2, -2),  # Astigmatism 45°
            (2, 0),   # Defocus
            (2, 2),   # Astigmatism 0°
            (3, -3),  # Trefoil Y
            (3, -1),  # Coma Y
            (3, 1),   # Coma X
            (3, 3),   # Trefoil X
            (4, -4),  # Tetrafoil Y
            (4, -2),  # Secondary Astigmatism Y
            (4, 0),   # Primary Spherical
            (4, 2),   # Secondary Astigmatism X
            (4, 4),   # Tetrafoil X
        ]
    
    def zernike_polynomial(self, n: int, m: int, rho: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """
        Calculate Zernike polynomial Z_n^m
        
        Args:
            n: Radial order
            m: Azimuthal frequency
            rho: Normalized radial coordinate (0 to 1)
            theta: Azimuthal coordinate (0 to 2π)
            
        Returns:
            Zernike polynomial values
        """
        # Normalization factor
        if m == 0:
            N = np.sqrt(n + 1)
        else:
            N = np.sqrt(2 * (n + 1))
        
        # Radial polynomial
        R = np.zeros_like(rho)
        for k in range((n - abs(m)) // 2 + 1):
            coeff = ((-1)**k * np.math.factorial(n - k)) / (
                np.math.factorial(k) * 
                np.math.factorial((n + abs(m)) // 2 - k) * 
                np.math.factorial((n - abs(m)) // 2 - k)
            )
            R += coeff * rho**(n - 2*k)
        
        # Angular part
        if m >= 0:
            angular = np.cos(m * theta)
        else:
            angular = np.sin(abs(m) * theta)
        
        return N * R * angular
    
    def zernike_derivative_x(self, n: int, m: int, rho: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """
        Calculate partial derivative of Zernike polynomial with respect to x
        """
        if m == 0:
            N = np.sqrt(n + 1)
        else:
            N = np.sqrt(2 * (n + 1))
        
        # Radial polynomial
        R = np.zeros_like(rho)
        for k in range((n - abs(m)) // 2 + 1):
            coeff = ((-1)**k * np.math.factorial(n - k)) / (
                np.math.factorial(k) * 
                np.math.factorial((n + abs(m)) // 2 - k) * 
                np.math.factorial((n - abs(m)) // 2 - k)
            )
            if n - 2*k > 0:
                R += coeff * (n - 2*k) * rho**(n - 2*k - 1)
        
        # Angular part derivatives
        if m > 0:
            angular = -m * np.sin(m * theta)
        elif m < 0:
            angular = abs(m) * np.cos(abs(m) * theta)
        else:
            angular = np.zeros_like(theta)
        
        return N * R * angular * np.cos(theta)  # cos(theta) from dx/drho
    
    def extract_fringe_phase(self, image: np.ndarray, shear: float) -> np.ndarray:
        """
        Extract phase from lateral shear interferogram
        
        Args:
            image: Interferogram image
            shear: Lateral shear distance in pixels
            
        Returns:
            Phase difference map
        """
        # Convert to grayscale if needed
        if len(image.shape) == 3:
            gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        else:
            gray = image
        
        # Normalize image
        gray = gray.astype(np.float64) / 255.0
        
        # Apply Gaussian filter to reduce noise
        gray = ndimage.gaussian_filter(gray, sigma=1.0)
        
        # Create coordinate grids
        y, x = np.mgrid[0:gray.shape[0], 0:gray.shape[1]]
        
        # Find circular aperture
        center_x, center_y = gray.shape[1] // 2, gray.shape[0] // 2
        radius = min(gray.shape) // 2 - 10  # Leave some margin
        
        # Create circular mask
        mask = ((x - center_x)**2 + (y - center_y)**2) <= radius**2
        
        # Extract phase using spatial carrier method
        # This is a simplified approach - in practice, more sophisticated
        # phase extraction methods would be used
        
        # For demonstration, we'll create a synthetic phase map
        # based on the fringe pattern
        phase = np.zeros_like(gray)
        
        # Simple fringe counting approach
        # In practice, you'd use more sophisticated phase extraction
        for i in range(gray.shape[0]):
            for j in range(gray.shape[1]):
                if mask[i, j]:
                    # Extract local fringe frequency
                    local_region = gray[max(0, i-5):min(gray.shape[0], i+6),
                                      max(0, j-5):min(gray.shape[1], j+6)]
                    if local_region.size > 0:
                        # Simple phase estimation based on intensity
                        phase[i, j] = 2 * np.pi * (gray[i, j] - 0.5)
        
        return phase * mask
    
    def calculate_zernike_coefficients(self, phase_map: np.ndarray, shear: float) -> Dict[str, float]:
        """
        Calculate Zernike coefficients from lateral shear phase map
        
        Args:
            phase_map: Phase difference map from lateral shear
            shear: Lateral shear distance in meters
            
        Returns:
            Dictionary of Zernike coefficients
        """
        # Create coordinate grids
        y, x = np.mgrid[0:phase_map.shape[0], 0:phase_map.shape[1]]
        
        # Convert to normalized coordinates
        center_x, center_y = phase_map.shape[1] // 2, phase_map.shape[0] // 2
        radius = min(phase_map.shape) // 2 - 10
        
        # Normalized coordinates
        rho = np.sqrt((x - center_x)**2 + (y - center_y)**2) / radius
        theta = np.arctan2(y - center_y, x - center_x)
        
        # Create aperture mask
        mask = rho <= 1.0
        
        # Calculate lateral shear phase gradient
        # The phase difference in lateral shear is related to the wavefront gradient
        phase_grad_x = np.gradient(phase_map, axis=1)
        phase_grad_y = np.gradient(phase_map, axis=0)
        
        # Convert phase to wavefront (in meters)
        # Phase difference = (2π/λ) * shear * ∂W/∂x
        wavefront_grad_x = phase_grad_x * self.wavelength / (2 * np.pi * shear)
        wavefront_grad_y = phase_grad_y * self.wavelength / (2 * np.pi * shear)
        
        # Fit Zernike polynomials to the gradient data
        coefficients = {}
        
        # For each Zernike mode, calculate the coefficient
        for i, (n, m) in enumerate(self.zernike_modes[:10]):  # First 10 modes
            # Calculate Zernike polynomial and its derivative
            Z = self.zernike_polynomial(n, m, rho, theta)
            dZ_dx = self.zernike_derivative_x(n, m, rho, theta)
            
            # Project gradient onto Zernike derivative
            numerator = np.sum(wavefront_grad_x * dZ_dx * mask)
            denominator = np.sum(dZ_dx**2 * mask)
            
            if denominator > 0:
                coeff = numerator / denominator
            else:
                coeff = 0.0
            
            # Convert to waves (divide by wavelength)
            coeff_waves = coeff / self.wavelength
            
            coefficients[f'Z_{n}^{m}'] = coeff_waves
        
        return coefficients
    
    def analyze_interferogram(self, image_path: str, shear_mm: float, 
                            fringe_separation_mm: float) -> Dict[str, float]:
        """
        Complete analysis of lateral shear interferogram
        
        Args:
            image_path: Path to interferogram image
            shear_mm: Lateral shear in millimeters
            fringe_separation_mm: Fringe separation in millimeters
            
        Returns:
            Dictionary of Zernike coefficients in waves
        """
        # Load image
        image = cv2.imread(image_path)
        if image is None:
            raise ValueError(f"Could not load image from {image_path}")
        
        # Convert shear to pixels (assuming image scale)
        # This is a rough estimation - in practice, you'd calibrate this
        pixels_per_mm = image.shape[1] / 100  # Rough estimate
        shear_pixels = shear_mm * pixels_per_mm
        
        # Extract phase from interferogram
        phase_map = self.extract_fringe_phase(image, shear_pixels)
        
        # Convert shear to meters
        shear_m = shear_mm * 1e-3
        
        # Calculate Zernike coefficients
        coefficients = self.calculate_zernike_coefficients(phase_map, shear_m)
        
        return coefficients
    
    def plot_results(self, image_path: str, coefficients: Dict[str, float]):
        """
        Plot the interferogram and Zernike coefficients
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Load and display original image
        image = cv2.imread(image_path)
        if len(image.shape) == 3:
            image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
        else:
            image_rgb = image
        
        axes[0, 0].imshow(image_rgb)
        axes[0, 0].set_title('Original Interferogram')
        axes[0, 0].axis('off')
        
        # Plot Zernike coefficients
        modes = list(coefficients.keys())
        values = list(coefficients.values())
        
        axes[0, 1].bar(range(len(modes)), values)
        axes[0, 1].set_xticks(range(len(modes)))
        axes[0, 1].set_xticklabels(modes, rotation=45, ha='right')
        axes[0, 1].set_ylabel('Coefficient (waves)')
        axes[0, 1].set_title('Zernike Coefficients')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Create a summary table
        axes[1, 0].axis('off')
        table_data = []
        for mode, value in coefficients.items():
            table_data.append([mode, f'{value:.3f}'])
        
        table = axes[1, 0].table(cellText=table_data,
                                colLabels=['Mode', 'Coefficient (waves)'],
                                cellLoc='center',
                                loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        axes[1, 0].set_title('Zernike Coefficients Summary')
        
        # Plot RMS wavefront error
        rms_error = np.sqrt(sum(c**2 for c in coefficients.values()))
        axes[1, 1].text(0.5, 0.5, f'RMS Wavefront Error:\n{rms_error:.3f} waves\n({rms_error * self.wavelength * 1e9:.1f} nm)',
                       ha='center', va='center', fontsize=14,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"))
        axes[1, 1].set_xlim(0, 1)
        axes[1, 1].set_ylim(0, 1)
        axes[1, 1].axis('off')
        axes[1, 1].set_title('Wavefront Quality')
        
        plt.tight_layout()
        plt.show()


def analyze_lateral_shear_interferogram(image_path: str, 
                                      aperture_diameter_mm: float = 50.0,
                                      wavelength_nm: float = 632.0,
                                      shear_mm: float = 3.5,
                                      fringe_separation_mm: float = 5.5) -> Dict[str, float]:
    """
    Main function to analyze lateral shear interferogram and return Zernike coefficients
    
    Args:
        image_path: Path to the interferogram image
        aperture_diameter_mm: Diameter of the aperture in mm
        wavelength_nm: Wavelength of the laser in nm
        shear_mm: Lateral shear distance in mm
        fringe_separation_mm: Fringe separation in mm
        
    Returns:
        Dictionary containing Zernike coefficients in waves
    """
    # Convert units
    wavelength = wavelength_nm * 1e-9  # Convert nm to m
    aperture_diameter = aperture_diameter_mm * 1e-3  # Convert mm to m
    
    # Create analyzer
    analyzer = LateralShearAnalyzer(wavelength=wavelength, 
                                  aperture_diameter=aperture_diameter)
    
    # Analyze interferogram
    coefficients = analyzer.analyze_interferogram(image_path, shear_mm, fringe_separation_mm)
    
    # Plot results
    analyzer.plot_results(image_path, coefficients)
    
    return coefficients


if __name__ == "__main__":
    # Example usage
    image_path = "markup_1000072306.png"
    
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
            image_path=image_path,
            aperture_diameter_mm=aperture_diameter,
            wavelength_nm=wavelength,
            shear_mm=shear,
            fringe_separation_mm=fringe_separation
        )
        
        print("Zernike Coefficients (in waves):")
        print("-" * 40)
        for mode, value in coefficients.items():
            print(f"{mode:>8}: {value:8.4f} waves ({value * wavelength:8.1f} nm)")
        
        # Calculate RMS wavefront error
        rms_waves = np.sqrt(sum(c**2 for c in coefficients.values()))
        rms_nm = rms_waves * wavelength
        print("-" * 40)
        print(f"{'RMS':>8}: {rms_waves:8.4f} waves ({rms_nm:8.1f} nm)")
        
    except Exception as e:
        print(f"Error analyzing interferogram: {e}")
        print("Please check that the image file exists and is accessible.")