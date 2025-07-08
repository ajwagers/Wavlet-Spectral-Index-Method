## Important scripts ##
### \wavelet\DISTRIB_V4\idl\pro\dmplot.pro ###
#### What dmplot.pro Does ###
This is a comprehensive scientific analysis script written in IDL. Its primary purpose is to analyze observational data from a set of Type Ia supernovae. It loads data, performs various calculations and statistical fits, and generates a large, multi-page PostScript file (dm15plot.ps) containing dozens of plots that explore the relationships between different supernova properties.

The script performs the following main tasks:

Data Loading and Preparation: It starts by defining two lists of supernova names (snname_new, snname_hamuy). It then iterates through these, calling a helper function getsninfo to retrieve photometric and spectroscopic data for each one. It also loads a pre-saved IDL data file (XWave.sav) which contains processed spectral data (like POW6150MAT, SIG6355Mat, etc.).

Data Filtering and Reshaping: The script filters the supernovae based on their redshift (z0) and other quality flags to create a clean sample for analysis. It uses the REFORM function extensively to reshape multi-dimensional data arrays.

Per-Supernova Plotting: It enters a large loop that iterates through each unique supernova. Inside this loop, it generates a series of plots for that specific object, showing how various spectral features (e.g., P!l6150!N, which likely refers to the pseudo-equivalent width of a spectral line at 6150 Å) evolve over time ('Day'). It also performs linear and spline fits to the data on these plots.

Group and Correlation Analysis: After the main loop, the script performs a broader analysis across the entire sample.

It groups supernovae into different classes (sn91bglike, sn91tlike, normal) and creates comparison plots of their velocity evolution.
It creates numerous correlation plots, comparing intrinsic properties like the decline rate (dm15) and absolute magnitude (M!lB!N) against various measured spectral features.
It uses functions like fitexy to perform robust linear fits to data that has errors in both the x and y dimensions.
In summary, this is a classic data analysis pipeline from an astronomical research project, designed to turn raw observational data into a series of publication-quality plots that investigate the physics of supernovae.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\xwhub20120421.pro ###
#### Purpose ####
This IDL procedure (pro xwhub20120421) plots a Hubble diagram for Type Ia supernovae (SNe Ia), comparing observed data to a theoretical cosmological model.

#### Step-by-step Breakdown ####
##### Restore Data #####

Loads data from a saved IDL structure/file named sn1a_hubble.sav.
The file contains arrays:
z: redshift values
mu: observed distance modulus
dmu: error in distance modulus
Set Cosmological Parameters

h0: Hubble constant (km/s/Mpc)
om: Omega matter (matter density parameter)
ol: Omega lambda (dark energy density parameter)

##### Calculate Theoretical Distance Modulus #####
For each redshift value, computes the theoretical distance modulus using a function distmod (not shown, but presumably computes the cosmological distance modulus for given parameters).

##### Plot Observed Data #####
Plots the observed data as points (psym=4 is a filled circle).
X-axis: redshift (z)
Y-axis: distance modulus (mu)
Sets axis ranges and labels.

##### Overplot Theoretical Curve #####
Plots the theoretical curve (distance modulus vs. redshift) in red over the data.

##### Add Error Bars #####
For each data point, draws a vertical blue line representing the error bar in distance modulus.

#### Summary ####
Loads SNe Ia Hubble diagram data (redshift, distance modulus, errors).
Computes the theoretical distance modulus for each redshift using cosmological parameters.
Plots the observed data with error bars.
Overplots the theoretical cosmological model for comparison.
In short:
This code visualizes how well the observed SNe Ia data fit a standard cosmological model on a Hubble diagram.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\waveletpca_20120405.pro ###

#### General Overview ####
The file is written in IDL (Interactive Data Language), commonly used for data analysis, visualization, and image processing, especially in astronomy and remote sensing.

The filename suggests it performs PCA (Principal Component Analysis) on wavelet-transformed data. The date in the filename (20120405) likely indicates its creation or last modification date.

#### Step-by-Step Explanation ####
1. Procedure Definition
The file likely starts with:

data: Input data array (probably 2D or 3D).
ncomp: Number of principal components to retain.
result: Output variable for the PCA result.
/PLOT: Optional keyword to plot results.
/HELP: Optional keyword to print usage/help.
2. Help Option
If /HELP is set, the procedure prints usage instructions and returns.

3. Wavelet Transform
The code likely applies a wavelet transform to the input data. This could be done using IDL’s WAVELET or similar routines, decomposing the data into different frequency components.

4. Data Reshaping
After the wavelet transform, the data is probably reshaped into a 2D array (samples × features) suitable for PCA.

5. Principal Component Analysis (PCA)
The code performs PCA, which involves:

Mean-centering the data.
Computing the covariance matrix.
Eigen-decomposition of the covariance matrix.
Sorting eigenvalues/eigenvectors to select the top ncomp components.
Projecting the data onto the principal components.
IDL has built-in functions like PCOMP for PCA.

6. Output
The result (principal components, projections, or reconstructed data) is stored in the result variable.

7. Optional Plotting
If /PLOT is set, the code generates plots of the principal components, eigenvalues, or projections for visualization.

#### Summary ####
waveletpca_20120405.pro is an IDL procedure that:

Takes input data.
Applies a wavelet transform to extract features at multiple scales.
Performs PCA on the transformed data to reduce dimensionality or extract main modes of variation.
Optionally plots the results.
Returns the PCA output.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\XWMC_20090128.pro ###

#### Purpose ####
This IDL file is a toolkit for wavelet analysis of supernova spectra, specifically for Type Ia supernovae. It automates the process of reading in spectral data, performing wavelet transforms, extracting statistical features (like variance, skewness, kurtosis, and integrated "power" in specific wavelength regions), and saving the results for further analysis or Monte Carlo simulations.

#### Key Components ####
1. MERGE_ARRAY Function
Purpose: Merges two arrays of possibly different sizes along their last dimension.
Use: Useful for combining simulation results or data blocks.
2. INT_ATROUS Function
Purpose: Integrates the wavelet power spectrum over a specified wavelength region, with options for error estimation.
Use: Quantifies the "power" (variance) in a spectral feature, e.g., Si II 6355 Å.
3. MR1D_ATROU Procedure
Purpose: Performs a 1D wavelet transform using the à trous algorithm.
Use: Decomposes a spectrum into multiple scales (resolutions) for feature extraction.
4. snatrous Procedure
Purpose: Reads a spectrum, pads it, applies the wavelet transform, and plots the results.
Use: For visual inspection and debugging of the wavelet decomposition.
5. extractdateinfo Procedure
Purpose: Extracts supernova name, observation date, and other info from a filename.
Use: Automates metadata extraction from file naming conventions.
6. xwavemc Procedure
Purpose: The core analysis routine.
Reads a spectrum file.
Adds simulated noise.
Performs wavelet transforms at multiple redshifts.
Extracts features (variance, skewness, kurtosis, integrated power) in key wavelength regions (e.g., 6355, 5485, 5150, 4570 Å).
Optionally plots results.
Writes results to a file.
Use: Used in Monte Carlo simulations to study how noise and other effects impact feature extraction from supernova spectra.
7. snxwave Procedure
Purpose: Batch driver for running xwavemc on many spectra and storing results.
Finds all relevant data files.
Loops over files and noise realizations.
Calls xwavemc and collects results into large arrays.
Saves the results as an IDL .sav file for later use.
Use: Automates large-scale simulations and feature extraction for statistical analysis.
#### Typical Workflow ####
Prepare a list of supernova spectra files.
For each file and for multiple noise realizations:
Add noise to the spectrum.
Perform wavelet decomposition.
Extract statistical features from specific wavelength regions.
Save the results.
Aggregate results for all files and realizations.
Save the aggregated results for further analysis (e.g., PCA, cosmological studies).
#### Scientific Context ####
The code is designed for quantitative analysis of supernova spectra using wavelet techniques.
It focuses on extracting robust, noise-insensitive features (variance, skewness, kurtosis, integrated power) from spectral lines that are important for supernova classification and cosmology.
The Monte Carlo approach (adding noise and repeating the analysis) helps assess the reliability of feature extraction under realistic observational conditions.
#### Summary ####
This file is a comprehensive IDL toolkit for simulating, analyzing, and extracting features from supernova spectra using wavelet transforms, with a focus on robust statistical analysis for large datasets and Monte Carlo studies.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\dm15sav.pro ###

#### Purpose ####
This script reads a data file containing information about Type Ia supernovae (SNe Ia), extracts relevant columns, and saves them into an IDL .sav file for later use in analysis.

#### Step-by-Step Explanation ####
1. Set Plot Device
```
SET_PLOT, 'WIN'
```
Sets the graphics device to 'WIN' (Windows display). This is not essential for the data processing, but may be a default setting.

2. Open and Read Data File
Opens the file dm15.dat for reading.
Reads the data as ASCII, using a template to parse the columns (comma-separated).
Closes the file.
```
dm15FILE = 'C:\TAMU_research\dm15.dat'
OPENR, 19, dm15FILE
dm15_data = READ_ASCII(DATA_START=0, DELIMITER=',', TEMPLATE=ASCII_TEMPLATE(dm15FILE))
CLOSE, 19
```

3. Extract Columns
Assigns each column of the data file to a variable with a descriptive name (e.g., supernova name, decline rate, redshift, photometric data, host galaxy info, etc.).
```
SN_name = dm15_data.field01
SN_dm15 = dm15_data.field02
SN_dm15source = dm15_data.field03
SN_datasource = dm15_data.field04
SN_z = dm15_data.field05
...
SN_cerr = dm15_data.field25
```

4. Get Number of Supernovae
```
SN_num = SIZE(SN_name)
```
Determines the number of supernovae (rows) in the data.

5. Save Data to .sav File
```
SAVE, SN_name, SN_dm15, SN_z, SN_edm15, SN_Bmax, SN_eBmax, SN_Bsource, SN_Vmax, SN_eVmax, SN_Vsource, SN_datasource, SN_dm15source, SN_A, SN_Asource, SN_EBV, SN_hEBV, $
      SN_hEBVerr, SN_hEBVsource, SN_zerr, SN_x1, SN_x1err, SN_c, SN_cerr, FILENAME = '\TAMU_research\dm15.sav', /VERBOSE
```
Saves all the extracted variables into an IDL .sav file for efficient future loading and analysis.

6. Print Completion Message
```
print, 'DONE!'
```

#### In Short ####
dm15sav.pro reads a CSV file of supernova properties, extracts the columns into variables, and saves them in an IDL .sav file for later use in scientific analysis.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\error_sav.pro ###

#### Purpose ####
This script processes a set of supernova spectral data files (for SN 2001el), bins the spectra, pads them for wavelet analysis, and saves the processed data into .sav files for each epoch (observation date).

#### Step-by-Step Explanation ####
1. Set Up File Paths
```
fileDIR = '\TAMU_research\SNSPEC\'
datFileName = fileDIR + '2001el\sn2001el*'
FileList = FINDFILE(datFileName)
NUMfiles = SIZE(FileList)
```
Sets the directory and file pattern for SN 2001el spectra.
Finds all matching files and counts them.

2. Loop Over Each File (Epoch)
```
FOR epochNUM = 0,NUMfiles(1)-1 DO BEGIN
```
Processes each file (epoch) one by one.

3. Extract Epoch Information from Filename
Uses string operations to find the epoch (days since maximum light) from the filename, handling both positive (p) and negative (m) epochs.

4. Read Spectrum Data
If FITS file:
Reads the spectrum and wavelength calibration from FITS headers.
If ASCII file:
Reads wavelength and flux columns from the file.

5. Bin the Data
Bins the spectrum into 5 Ångström intervals using a simple average, to standardize the data for further analysis.

6. Pad Data for Wavelet Analysis
Pads the binned data to the next power of two (4096 points), filling extra points with a decaying exponential to minimize edge effects in wavelet transforms.

7. Save the Processed Data
```
SAVE, Nlam, NbinData, lam_range, sn1a_lammax, sn1a_lammin, nnnBIN, delta_lam, FILENAME='/TAMU_research/error_sav/2001el.'+STRING(epoch,FORMAT='(I03)')+'.sav'
```
Saves the wavelength array, binned flux, and metadata into a .sav file named for the epoch.

8. Completion Message
```
print, 'DONE!'
```

#### In Short ####
error_sav.pro automates the preparation of supernova spectra for wavelet analysis by:

Reading each spectrum,
Binning and padding the data,
Saving the result for each epoch in a standardized format.
This is a typical preprocessing step for further time-series or spectral feature analysis.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\wimspec.pro ###
#### Purpose ####
This script is designed to process, organize, and plot supernova spectra from various sources and categories. It reads in spectral data, organizes it by type and epoch, performs some wavelet analysis, and generates a series of publication-quality plots (PostScript files) for different supernova groups (e.g., Type Ia, non-Ia, early, late, etc.).

#### Step-by-Step Overview ####
1. File Discovery and Setup
Finds all relevant spectral data files in specified directories (e.g., Hsiao templates, observed spectra, early/late/very late/NOT Ia/garbage spectra).
Initializes arrays to store filenames, dates, and other metadata.

2. Hsiao Template Processing
Reads in Hsiao template spectra (standardized SN Ia templates).
Extracts epoch information from filenames.
Sorts templates by epoch.
Reads and stores wavelength and flux data for each template.

3. Observed Spectra Processing
Reads in observed spectra filenames and extracts SN names and observation dates.
Matches each spectrum to its maximum-light date using a lookup file (maxdates.txt).

4. Classification
Classifies spectra into different groups:
Type Ia (normal SNe Ia)
NOT Ia (other types: Ib, Ic, II, etc.)
Early (pre-maximum)
Late (post-maximum)
Very Late
Garbage (problematic or unclassified spectra)

5. Wavelet Analysis and Cross-Correlation
For each Type Ia spectrum:
Reads and bins the spectrum.
Pads the data for wavelet analysis.
Performs a wavelet transform (à trous algorithm).
Sums selected wavelet scales to emphasize certain features.
Cross-correlates the observed spectrum with Hsiao templates over a range of redshifts to find the best match.
Plots the results, including the best-matching template and the observed spectrum.

6. Plotting
For each group (Ia, NOT Ia, Early, Late, Very Late, Garbage):
Generates a stacked plot of all spectra in the group, normalized and offset for clarity.
Annotates each spectrum with its name, date, and epoch.
Saves each plot as a PostScript file in the appropriate directory.

7. Completion
Prints 'DONE!' when finished.

#### In Short ####
wimspec.pro is a comprehensive script for:
Organizing and classifying supernova spectra,
Performing wavelet-based feature extraction and template matching,
Generating detailed, annotated plots for different supernova groups.
It is a key tool for visualizing and comparing supernova spectral data in a research context.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\xwavenoise.pro ###

#### Purpose ####
This script simulates the effect of noise on a supernova spectrum and visualizes the result. It adds Poisson noise to the flux values in a specified wavelength region and plots both the original and noisy spectra for comparison.

#### Step-by-Step Explanation ####
1. Procedure Definition
```
PRO XWAVENOISE, lamval, fluxval, snname, snepoch, q
```
Inputs:
lamval: Array of wavelength values.
fluxval: Array of flux values (spectrum).
snname: Supernova name (string).
snepoch: Epoch (integer, e.g., days since max light).
q: Noise scaling factor.

2. Set Plotting to PostScript
```
SET_PLOT, 'PS'
DEVICE, FILENAME = '\TAMU_research\results\error_out' + STRING(snname) + STRING(snepoch,FORMAT='(I03)') + '.ps'
```
Sets the output to a PostScript file named with the supernova and epoch.

3. Noise Calculation and Addition
```
nfluxval = FLTARR(N_ELEMENTS(fluxval))
w6355 = WHERE(lamval GE 5500 AND lamval LE 6500)
IF w6355 NE [-1] THEN BEGIN
  sig3 = STDDEV(fluxval(w6355))
  sig1 = STDDEV(fluxval(w6355))
  rho13 = sig3/sig1
  numlamval = N_ELEMENTS(lamval)
  rnum = RANDOMN(42,numlamval,POISSON=sig3, /DOUBLE)
  nfluxval = fluxval + rnum*q
ENDIF
```
Identifies the wavelength region 5500–6500 Å (often associated with the Si II 6355 feature in SN Ia spectra).
Calculates the standard deviation (sig3) of the flux in this region.
Generates Poisson-distributed random noise (rnum) with the same standard deviation.
Adds the scaled noise (q is the scaling factor) to the original flux to create a noisy spectrum (nfluxval).

4. Plotting
```
PLOT, lamval, fluxval, TITLE = snname + ' ' + STRING(snepoch,FORMAT='(I03)') + ' noise test'
OPLOT, lamval, nfluxval, COLOR = 50
DEVICE, /CLOSE
```
Plots the original spectrum.
Overplots the noisy spectrum in a different color.
Closes the PostScript device.

5. Signal-to-Noise Calculation
```
STNR = TOTAL(fluxval)/TOTAL(rnum)
```

6. Reset Plotting to Screen
```
SET_PLOT, 'X'
```
#### In Short ####
xwavenoise.pro takes a supernova spectrum, adds simulated Poisson noise to it, and plots both the original and noisy spectra for visual comparison, saving the result as a PostScript file.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\hubpcafit20110604.pro ###

#### Purpose ####
This script performs a Principal Component Analysis (PCA) on Type Ia supernova data, focusing on the relationship between Hubble residuals and various supernova parameters (such as color, light-curve shape, and spectral features). It then fits and visualizes how these parameters (and their principal components) relate to the observed Hubble diagram, aiming to improve the standardization of SNe Ia as distance indicators.

#### Key Components ####
1. Utility Functions
- *COSMOFUNC, SIMPCOSMO, SIMPCOSMO2:*
Functions to compute theoretical distance modulus as a function of redshift and cosmological parameters.
- *MYLINE, MYQUAD:*
Simple linear and quadratic model functions.
- *PCAnalysis:*
A custom function to perform PCA on a 2D data array, with options for standardization and error estimation via bootstrapping.

2. Main Procedure: HUBPCAfit
- Loads supernova names and groups (CN, BL, CL, SS, KP, UK), which are different subclasses of SNe Ia.
- Loads supernova data from a .sav file (e.g., redshift, B-band maximum, color, DM15, etc.).
- Loads spectral feature measurements (chi values at different wavelengths and epochs) for each supernova from individual .sav files.
- Computes mean and standard deviation of these features over a specified epoch range for each SN.
- Filters the sample to select SNe Ia with good data and within certain parameter ranges (e.g., redshift, DM15, etc.).
- Computes Hubble residuals (difference between observed and expected B-band magnitude).
- Fits a simple cosmological model to the data using MPFITFUN.
- Performs PCA on a matrix of Hubble residuals and SN parameters (color, DM15, spectral features).
- Extracts and prints the relationships (fit parameters) between Hubble residuals and the principal components, color, and DM15.
- Plots:
    - The Hubble diagram and residuals.
    - The effect of correcting the Hubble residuals using PC1, color, DM15, and combinations.
    - The relationships between standardized parameters and the Hubble residuals.
- Calculates and prints the scatter (standard deviation) of the Hubble residuals before and after corrections.

#### Scientific Context ####
*Goal:*
To improve the precision of Type Ia supernovae as standard candles by identifying and correcting for additional sources of scatter in the Hubble diagram, using both traditional parameters (color, DM15) and spectral features via PCA.
*Approach:*
By projecting the data onto principal components, the script identifies the combinations of parameters that best explain the variance in Hubble residuals, potentially leading to better standardization and cosmological constraints.

#### In Short ####
hubpcafit20110604.pro is a comprehensive analysis script for:

Loading and organizing SN Ia data,
Extracting and averaging spectral features,
Computing Hubble residuals,
Fitting cosmological and empirical models,
Performing PCA to find optimal corrections,
Visualizing and quantifying the improvement in SN Ia standardization.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\bailey20110615.pro ###

#### Purpose ####
This script analyzes and visualizes the relationship between supernova spectral features and their peak brightness (Bmax). It processes spectra for a set of Type Ia supernovae, bins the spectra into velocity-based wavelength intervals, normalizes them, and computes correlations between binned fluxes and Bmax across the sample. The results are visualized as plots and correlation matrices.

#### Key Components ####
1. Function: binspec_bl\
**Purpose:**\
Bins a spectrum into intervals defined by a constant velocity width (rather than constant wavelength width).\
**How:**\
- Starts at 3500 Å and increments by a velocity-based step (using vlam).
- Extends the input wavelength and flux arrays if needed to cover the output range.
- Interpolates the input flux onto the new wavelength grid using spline interpolation.
- Integrates the flux within each bin to produce the binned spectrum.
Returns the binned flux array.
2. Procedure: bailey\
**Purpose:**\
Main analysis routine for the Bailey et al. (2009) style spectral binning and correlation analysis.\
**Steps:**\
    1. Define Supernova Groups:
Lists of SN names grouped by spectral subclass (e.g., "PC", "HV", "FD", etc.).
    2. Create Output Wavelength Grid:
Sets up a velocity-binned wavelength grid from 3500 Å to 8500 Å.
    3. Initialize Arrays:
Prepares arrays to hold binned wavelengths, fluxes, and other SN properties.
    4. Load SN Properties:
Loads SN parameters (redshift, dm15, Bmax, etc.) from a .sav file.
    5. Match SN Names to Properties:
For each SN in the list, finds and stores its properties from the loaded data.
    6. (Optional) Redshift Adjustment and Plotting:
If the ZADJUST keyword is set, for each SN:
        - Finds the spectrum closest to maximum light.
        - Loads and bins the spectrum, as well as template spectra.
        - Performs wavelet analysis and cross-correlation with templates over a range of redshifts.
        - Plots the normalized spectra and saves the results.
        - Stores the binned, normalized spectra for later analysis.
    7. Correlation Analysis:
        - For each wavelength bin, computes the correlation between the binned flux and Bmax across all SNe.
        - Plots the correlation as a function of wavelength.
        - Computes and visualizes a 2D correlation matrix between all pairs of bins and Bmax using images and contour plots.

#### In Short ####
bailey20110615.pro implements a spectral binning and correlation analysis for Type Ia supernovae:

Bins spectra by velocity intervals,
Normalizes and stores the binned spectra,
Computes how each spectral bin correlates with peak brightness (Bmax),
Visualizes these correlations as plots and 2D matrices.
This approach helps identify which spectral regions are most predictive of SN Ia luminosity, aiding in standardization and cosmological studies.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\sn1a_template_20110517.pro ###

#### Purpose ####
This script is a comprehensive tool for analyzing Type Ia supernova (SN Ia) spectra using à trous wavelet decomposition. It processes a large set of SN Ia spectra, performs wavelet analysis, extracts spectral features, normalizes and bins the data, and saves the results for further analysis. It also generates plots and computes "spectral indices" (feature strengths) as a function of time (epoch) and compares them to light-curve parameters like Dm15.

#### Key Steps and Components ####
1. Setup and Initialization
Sets up plotting parameters and line thicknesses for publication-quality plots.
Defines lists of supernova names to process, including various subgroups and special cases.
Loads supporting data files, such as Dm15 values and redshifts.

2. Main Loop Over Supernovae
For each supernova in the list:
Finds all spectral files for that SN (across different epochs).
Loops over each epoch (spectrum) for the SN.

3. Spectrum Reading and Preprocessing
Reads the spectrum from either FITS or ASCII files.
Bins the spectrum into 5 Å bins.
Pads the binned data to the next power of two (for wavelet analysis), filling extra points with a decaying exponential to minimize edge effects.

4. Redshift Correction
Adjusts the wavelength scale for redshift, using either a precomputed value or by cross-correlation with template spectra (SN 1998aq, 1998bu, 2001el).
Stores the best-fit redshift for each epoch and SN.

5. Wavelet Decomposition
Performs à trous wavelet decomposition on the binned spectrum and on versions with different extinction corrections.
Sums selected wavelet scales to isolate features of interest.

6. Feature Extraction and Normalization
Calculates the standard deviation of the wavelet sum in specific wavelength regions (e.g., around 6355, 5485, 5150, 4570 Å) to quantify feature strengths.
Normalizes the wavelet coefficients.

7. Saving Results
Saves the processed data, wavelet coefficients, and feature strengths to .sav files for each SN and epoch.

8. Plotting
Generates and saves plots of the original spectrum, wavelet scales, and wavelet sums for each SN and epoch.
Plots the time evolution of spectral indices for each feature and SN.
Plots the relationship between spectral indices at maximum light and Dm15 (a light-curve decline rate parameter).

9. Spectral Index Calculation
For each feature and epoch, calculates a "spectral index" (CHI) using the wavelet sum, normalized by the standard deviation in the feature region.
Saves these indices for further analysis.

10. Summary Plots
Plots the time evolution of spectral indices for all SNe and for selected SNe.
Plots spectral index vs. Dm15 for each feature.

#### Scientific Context #####
- Wavelet analysis allows for robust, scale-dependent feature extraction from noisy spectra.
- Spectral indices quantify the strength of key features (e.g., Si II 6355 Å) as a function of time and SN properties.
- Comparison to Dm15 helps relate spectral features to light-curve shape, which is important for SN Ia standardization in cosmology.
Summary Table
Step	Purpose
Setup	Define SN list, load parameters, set plotting
Spectrum reading	Read and bin spectra, pad for wavelet analysis
Redshift correction	Adjust wavelength scale using templates or saved values
Wavelet analysis	Decompose spectra, sum scales for features
Feature extraction	Measure feature strengths (spectral indices)
Saving	Store results for each SN/epoch/feature
Plotting	Visualize spectra, wavelets, feature evolution, correlations

####In Short ####
sn1a_template_20110517.pro is a full pipeline for:

Reading and preprocessing SN Ia spectra,
Performing wavelet-based feature extraction,
Quantifying and saving feature strengths,
Plotting the evolution and correlations of these features with SN properties.

### C:\Users\ajwag\Documents\Code\TAMU\TAMU_research_old\IDL_Archive\evofit3_20100709.pro ###

#### High-Level Summary ####
This script, evofit3_20100709.pro, is designed to analyze the properties of spectral features in a large sample of supernovae. Its main goals are:

1. *Process and Correct Data:* It loads pre-processed supernova data and applies a sophisticated, custom noise-correction model.
2. *Measure Feature Evolution:* It fits mathematical models to the time evolution of several key spectral features for each supernova.
3. *Classify Supernovae:* It uses the measured feature properties to perform cluster analysis, objectively grouping the supernovae based on their similarities.
4. *Generate Publication-Quality Outputs:* It creates a large number of plots (scatter plots, evolution plots, dendrograms) and data tables intended for a scientific paper. The file paths like paper2/noiseplots3_b.ps strongly suggest this. 

In essence, it's a tool to take raw measurements, clean them up, extract meaningful physical parameters, and explore the relationships between those parameters to understand the diversity of supernovae.

#### Detailed Breakdown ####
The script can be divided into several logical parts:

1. Function Definitions (Lines 1-20)
idl
```
 Show full code block 
FUNCTION MYFUNC, X, P
  RETURN, P[0]/(P[3]+P[5]*EXP((P[1]-X)/P[2]))+P[4]
END

FUNCTION SIIFUNC, X, P
  RETURN, P[0]*(1/(1+EXP((P[1]-X)/P[2]))-1)
END
...
```
The script starts by defining several mathematical functions (MYFUNC, SIIFUNC, EFFUNC, etc.). These are various forms of logistic or sigmoid functions, which are commonly used to model processes that show an "S-shaped" curve—like the changing strength of a spectral feature over time. These functions are used later in the script to fit the observational data.

2. Clustering Procedures (Lines 22-66)
idl
``` 
PRO clustering, X, names, colors, lblnames, lfnds, VERBOSE=verbose, LINKAVE=linkave, FAR=far;, CLUSRES
...
END

PRO clusterzoom, X, names
...
END
```
These are helper procedures that perform and visualize hierarchical cluster analysis.

The clustering procedure takes a dataset X, calculates the "distance" between each item, groups them into a cluster tree, and then creates a dendrogram (a tree-like diagram) to show the groupings.
It's a key statistical method used later in the script to find natural groupings within the supernova sample based on their measured properties.

3. Main Program: evofit3 (Lines 68-end)
This is the heart of the script where the main analysis happens.

    A. Data Loading and Preparation (Lines 68-168)

        - Loading Data: It uses the RESTORE command to load data from several .sav files, such as dm15.sav and XWMC.sav. These files likely contain:
            - dm15: A key supernova parameter (the decline in brightness 15 days after maximum).
            - Wavelet power (POW...) and noise (SIG...) measurements for different spectral features, identified by their wavelengths (e.g., 6355 for the Si II 6355Å line, 5485 for the Si II 5485Å line, etc.).
            - Data for NS supernovae, with NMC Monte Carlo simulations for each to estimate errors.
        - Noise Correction: It performs a series of calculations to model and correct for noise in the data. It does linear fits (linfit) between noise measurements at different scales to create a robust correction factor (c1sm). This correction is crucial for comparing features between different observations and supernovae.

    B. Analysis and Plot Generation (Multiple large loops)

    The rest of the script is a series of very large FOR loops that generate the plots and data products for the paper.

    - Noise Plots (Lines 170-305): The first two sets of plots investigate the noise characteristics of the data and validate the correction model derived in the previous step.

    - Feature Correlation Plots (Lines 307-2411): This is the script's largest section, containing three nearly identical, massive, nested loops.

        - The loops iterate through every possible pair of spectral features (e.g., 5800Å vs. 6100Å).
        - For each supernova, it:
            1. Loads the data for the two features.
            2. Applies the noise correction: chi3a_sav1 = chi3a_sav(epsav_sort)/SQRT(1-c1sm*rhoval6355(epsav_sort)^2).
            3. Fits the time evolution of the corrected feature strength (chi3a_sav) using POLY_FIT or LINFIT.
            4. Determines the feature's strength at maximum light (epoch = 0), either from a direct measurement or the fit.
            5. Plots the strength of feature 1 vs. the strength of feature 2 on a scatter plot. Each point is one supernova.
        - The points are color-coded based on pre-defined supernova types (snNAMES_CN, snNAMES_BL, etc.), which correspond to known classes like Core-Normal, Broad-Line, Cool, and Shallow-Silicon. This allows the author to see how these classes relate to the measured feature strengths.

    - Clustering and Dendrograms (Lines 2413-2479):

        - After calculating the feature strengths for all supernovae, it takes a subset of these parameters (the strengths of the 5800Å and 6100Å features, and Dm15).
        - It standardizes this data (subtracts the mean and divides by the standard deviation).
        - It then calls the clustering procedure to group the supernovae and plots the resulting dendrograms. This is a data-driven way to see which supernovae are most similar to each other.

    - Final Evolution Plots and Data Saving (Lines 2481-end):

        - The final sections generate more plots, this time showing the fitted evolution curves for each feature, grouped by supernova type.
        - It uses MPFITFUN with the custom functions from the top of the file to perform more complex, non-linear fits to the data.
        - Finally, it saves the results of all these fits into .sav files and a LaTeX-formatted table using PRINTF, ready for the final publication.

In summary, this is a powerful script that automates a very complex astronomical data analysis workflow, from raw measurements to final, publication-ready scientific results.
