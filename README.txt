Calibration and Demodulation Instructions

#. Repository Structure

The repository contains MATLAB scripts, calibration files, and experimental spectral datasets for three multimode–no-core–multimode (MNM) fiber-optic sensors:

- **MNM-1:** Experimental datasets identified by the suffix `960`.
- **MNM-2:** Experimental datasets identified by the suffix `965`.
- **MNM-3:** Experimental datasets identified by the suffix `977`.

The repository is organized into three main directories:

1. MATLAB_CODE

Contains the MATLAB scripts:

- `svdsystem.m`: SVD-based spectral calibration and validation.
- `demodula.m`: Spectral demodulation and measurand identification.
- `simulatedmnm.m`: Numerical simulation of the MNM spectral response.

2. CALIBRATION_DATA

Contains the calibration models generated for MNM3:

- `calibracioncurvatura_977.mat`
- `calibracionstrain_977.mat`
- `calibraciontemperatura_977.mat`

3. EXPERIMENTAL_DATA

Contains the experimental spectra organized by measurand and sensor:

```text
4. EXPERIMENTAL_DATA/
│
├── Curvature_response/
│   ├── CURV960/       (MNM-1)
│   ├── CURV965/       (MNM-2)
│   └── CURV977/       (MNM-3)
│
├── Strain_response/
│   ├── STRAIN960/     (MNM-1)
│   ├── STRAIN965/     (MNM-2)
│   └── STRAIN977/     (MNM-3)
│
└── Temperature_response/
    ├── TEMP960/       (MNM-1)
    ├── TEMP965/       (MNM-2)
    └── TEMP977/       (MNM-3)
```

Each experimental dataset contains the spectral measurements associated with the corresponding sensor and measurand.

##. SVD-Based Calibration (`svdsystem.m`)

To reproduce the calibration procedure, open MATLAB, set the working directory to the repository's root folder, and run `svdsystem.m`. No user input is required.

The following procedure describes the calibration of **MNM-3 (977)**. The corresponding experimental spectra are located in:

- `4. EXPERIMENTAL_DATA/Curvature_response/CURV977/`
- `4. EXPERIMENTAL_DATA/Strain_response/STRAIN977/`
- `4. EXPERIMENTAL_DATA/Temperature_response/TEMP977/`

Each dataset contains five measurement rounds, organized into subfolders from `1ra medicion` to `5ta medicion`. The applied curvature, strain, or temperature values are automatically extracted from the filenames.

Using the first three measurement rounds (T1–T3), the script constructs a rank-four singular value decomposition (SVD) basis and sixth-order polynomial calibration trajectories over the wavelength range of 840–1000 nm. The remaining rounds (T4 and T5) are excluded from calibration and used for independent validation.

For each measurand, the script generates ten figures (101–110 for curvature, 201–210 for strain, and 301–310 for temperature), illustrating:

- Calibration spectra in linear and logarithmic (dB) scales.
- The four spectral components obtained from SVD.
- Component amplitudes.
- Normalized coefficients (A2/A1, A3/A1, and A4/A1) and their sixth-order polynomial fits.
- Three-dimensional calibration trajectories.
- Estimated versus applied values for the independent validation rounds (T4 and T5).

The script also displays the root mean square error (RMSE) for the validation rounds, together with the coefficient of determination (R²) and standard deviation (σ) obtained from linear regression of the calibration data.

The resulting MNM3 calibration models are saved as:

- `calibracioncurvatura_977.mat`
- `calibracionstrain_977.mat`
- `calibraciontemperatura_977.mat`

These calibration files are also provided in the `3. CALIBRATION_DATA/` directory.

**Calibration of MNM-1 and MNM-2**

The same calibration methodology can be applied to the other two sensors by selecting their corresponding experimental datasets:

| Sensor | Curvature | Strain      | Temperature |
--------------------------------------------------
| MNM1   | `CURV960` | `STRAIN960` | `TEMP960`   |
| MNM2   | `CURV965` | `STRAIN965` | `TEMP965`   |
| MNM3   | `CURV977` | `STRAIN977` | `TEMP977`   |

The sensor-specific datasets must be selected in the script before running the corresponding calibration. The provided MNM3 procedure serves as the reference example for reproducing the SVD-based calibration and validation methodology.

###. Spectral Demodulation (`demodula.m`)

After completing the calibration procedure, run `demodula.m` to analyze experimental spectra.

The script prompts the user to select the three previously generated calibration files in the following order:

1. Curvature calibration (`calibracioncurvatura_977.mat`)
2. Strain calibration (`calibracionstrain_977.mat`)
3. Temperature calibration (`calibraciontemperatura_977.mat`)

For MNM3, these files are available in the `3. CALIBRATION_DATA/` directory.

The user can then select one or multiple experimental spectra for analysis.

For each spectrum, the algorithm identifies the most likely measurand by determining which calibration model produces the lowest normalized residual. It then estimates the corresponding measurand value.

When the applied value can be extracted from the filename, the script also reports the reference value and estimation error. All results are compiled into a summary table.

Additionally, the last analyzed spectrum is visualized relative to the three calibration trajectories, together with the normalized residual associated with each calibration model.

For MNM1 and MNM2, the corresponding sensor-specific calibration files should be generated and selected before demodulation.