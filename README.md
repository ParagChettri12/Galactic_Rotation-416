# Santi Parag Isaac - 416 Project - An Analysis of Milkey Way's Rotation Curve

This Pluto.jl notebook presents an analysis of the Milky Way's rotation curve using data from Gaia Data Release 3 (DR3) compared with DR2. The project was developed by Parag Chettri, Santiago Berumen, and Isaac Whitson as part of a computational astrophysics study. Our goal was to use Gaia's data to effectively find the models that describe the rotation curve of our galaxy and the parameters that determine the results we were finding.

Key Features:

 - Interactive data exploration with Pluto.jl
 - Galactic coordinate transformations
 - Kinematic modeling of stellar motions
 - Comparison of simple and complex galactic rotation models
 - Visualization of rotation curves and residual analyses

Data Sources:
   - Primary dataset: Gaia DR3 (training data)
   - Validation dataset: Gaia DR2 (test data)

Methodology:
 - Data querying and filtering of Gaia sources near the galactic plane
 - Coordinate transformations from equatorial to galactocentric frames
 - Velocity calculations including proper motion corrections
 - Keplerian and dark matter-corrected rotation curve modeling
 - Nonlinear regression fitting to determine galactic mass components

Notebook Structure
  The notebook is organized into several key sections:
  - Data loading and preprocessing
  - Coordinate transformations
  - Velocity calculations
  - Rotation curve modeling
  - Model fitting and validation
  - Visualization of results
