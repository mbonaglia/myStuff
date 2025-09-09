# Lateral Shear Interferometry Analysis

This repository contains a comprehensive analysis system for lateral shear interferograms, based on Malacara's "Optical Shop Testing" Chapter 4. The system extracts Zernike polynomial coefficients and calculates beam divergence from experimental parameters.

## Overview

Lateral shear interferometry is a powerful technique for evaluating optical wavefront quality by comparing a wavefront with a laterally shifted version of itself. This self-referencing method is particularly useful for testing optical components without requiring an external reference wavefront.

## Features

- **Zernike Polynomial Analysis**: Extracts the first 10 Zernike coefficients (up to spherical aberration)
- **Beam Divergence Calculation**: Computes far-field beam divergence from defocus measurements
- **Multiple Analysis Methods**: Provides both simplified and refined analysis approaches
- **Comprehensive Output**: Includes RMS wavefront error, Strehl ratio, and quality assessments
- **Easy-to-Use Interface**: Simple function calls with experimental parameters

## Files Description

### Main Analysis Scripts

1. **`zernike_analyzer.py`** - Main, easy-to-use function for Zernike coefficient extraction
2. **`beam_divergence_calculator.py`** - Dedicated beam divergence analysis
3. **`complete_analysis.py`** - Comprehensive analysis combining both Zernike and divergence calculations
4. **`final_lateral_shear_analysis.py`** - Most accurate analysis method with detailed physics
5. **`simple_lateral_shear_analysis.py`** - Simplified analysis for basic understanding
6. **`lateral_shear_analysis.py`** - Full-featured analysis with image processing capabilities

### Data Files

- **`markup_1000072306.png`** - Sample lateral shear interferogram
- **`Malacara - Optical shop testing (1).pdf`** - Reference textbook
- **`*_results.json`** - Analysis results in JSON format

## Quick Start

### Basic Usage

```python
from zernike_analyzer import analyze_lateral_shear_interferogram

# Your experimental parameters
coefficients = analyze_lateral_shear_interferogram(
    shear_mm=3.5,
    fringe_separation_mm=5.5,
    aperture_diameter_mm=50.0,
    wavelength_nm=632.0
)

print(coefficients)
```

### Complete Analysis

```python
from complete_analysis import complete_lateral_shear_analysis

# Get comprehensive results including beam divergence
results = complete_lateral_shear_analysis(
    shear_mm=3.5,
    fringe_separation_mm=5.5,
    aperture_diameter_mm=50.0,
    wavelength_nm=632.0
)

print(results['beam_divergence'])
print(results['zernike_coefficients'])
```

## Example Results

Based on the sample interferogram with parameters:
- **Aperture diameter**: 50 mm
- **Wavelength**: 632 nm
- **Lateral shear**: 3.5 mm
- **Fringe separation**: 5.5 mm

### Zernike Coefficients (in waves):
- **Z₂⁰ (Defocus)**: 3.57 waves (2,257 nm)
- **Z₂⁻² (Astigmatism 45°)**: 0.36 waves (226 nm)
- **Z₂² (Astigmatism 0°)**: 0.36 waves (226 nm)
- **Z₁⁻¹ (Tilt Y)**: 0.18 waves (113 nm)
- **Z₁¹ (Tilt X)**: 0.18 waves (113 nm)

### Beam Divergence:
- **0.18 mrad** (milliradians)
- **0.0103°** (degrees)
- **37.25 arcsec** (arcseconds)

### Wavefront Quality:
- **RMS Error**: 3.62 waves (2,285 nm)
- **Strehl Ratio**: < 0.001
- **Quality Assessment**: Poor (high defocus)
- **Divergence Level**: Low (good beam propagation)

## Theory

### Lateral Shear Interferometry

In lateral shear interferometry, the phase difference between sheared wavefronts is:

```
Δφ = (2π/λ) × shear × ∂W/∂x
```

Where:
- `Δφ` is the phase difference
- `λ` is the wavelength
- `shear` is the lateral shear distance
- `∂W/∂x` is the wavefront gradient

### Zernike Polynomials

The system uses the OSA/ANSI standard Zernike polynomials to describe wavefront aberrations:

- **Z₀⁰**: Piston
- **Z₁⁻¹, Z₁¹**: Tilt (Y, X)
- **Z₂⁻², Z₂⁰, Z₂²**: Astigmatism (45°), Defocus, Astigmatism (0°)
- **Z₃⁻³, Z₃⁻¹, Z₃¹, Z₃³**: Trefoil (Y), Coma (Y), Coma (X), Trefoil (X)

### Beam Divergence Calculation

The beam divergence is calculated from the defocus coefficient:

```
θ = aperture_diameter / (2 × radius_of_curvature)
```

Where the radius of curvature is derived from the lateral shear measurements.

## Requirements

- Python 3.6+
- Standard library modules (math, json, typing)
- Optional: numpy, matplotlib, opencv-python, scipy (for advanced analysis)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/mbonaglia/myStuff.git
cd myStuff
```

2. Run the analysis:
```bash
python3 zernike_analyzer.py
python3 complete_analysis.py
```

## Usage Examples

### Command Line Usage

```bash
# Basic Zernike analysis
python3 zernike_analyzer.py

# Complete analysis with beam divergence
python3 complete_analysis.py

# Beam divergence only
python3 beam_divergence_calculator.py
```

### Programmatic Usage

```python
# Import the analysis functions
from zernike_analyzer import analyze_lateral_shear_interferogram, print_zernike_results
from complete_analysis import complete_lateral_shear_analysis

# Analyze your interferogram
coefficients = analyze_lateral_shear_interferogram(
    shear_mm=3.5,
    fringe_separation_mm=5.5,
    aperture_diameter_mm=50.0,
    wavelength_nm=632.0
)

# Print formatted results
print_zernike_results(coefficients)

# Get complete analysis
results = complete_lateral_shear_analysis(
    shear_mm=3.5,
    fringe_separation_mm=5.5,
    aperture_diameter_mm=50.0,
    wavelength_nm=632.0
)

# Access specific results
print(f"Beam divergence: {results['beam_divergence']['divergence_mrad']:.2f} mrad")
print(f"RMS wavefront error: {results['wavefront_quality']['rms_waves']:.3f} waves")
```

## Output Files

The analysis generates several output files:

- **`*_results.json`**: Analysis results in JSON format
- **Console output**: Formatted results with interpretations
- **Quality assessments**: Wavefront quality and beam divergence evaluations

## References

1. Malacara, D. "Optical Shop Testing" - Chapter 4: Lateral Shear Interferometers
2. Zernike polynomial definitions (OSA/ANSI standard)
3. Lateral shear interferometry theory and applications

## Contributing

Feel free to submit issues, feature requests, or pull requests to improve the analysis capabilities.

## License

This project is part of the myStuff repository. Please refer to the main repository license.

## Contact

For questions about the lateral shear interferometry analysis, please open an issue in the repository.

---

**Note**: This analysis is based on the physics of lateral shear interferometry as described in Malacara's textbook. For the most accurate results, full image analysis of the interferogram would be recommended, but this parameter-based approach provides good estimates for most applications.