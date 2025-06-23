import numpy as np
from typing import Tuple
from scipy.signal import hilbert
from scipy.special import gamma, psi, factorial, kv

from .transforms import compute_hanning_power_spec
from .decorators import output_vars, required_kwargs, dependencies
from ...data_objs.analysis_config import RfAnalysisConfig
from ...data_objs.analysis import Window
from ...data_objs.image import UltrasoundRfImage

@output_vars("f", "nps", "r_ps", "ps")
def compute_power_spectra(scan_rf_window: np.ndarray, phantom_rf_window: np.ndarray, 
                    window: Window, config: RfAnalysisConfig, 
                    image_data: UltrasoundRfImage, **kwargs) -> None:
    """Compute power spectra for a single window.
    """
    f, ps = compute_hanning_power_spec(
        scan_rf_window, config.transducer_freq_band[0],
        config.transducer_freq_band[1], config.sampling_frequency
    )
    ps = 20 * np.log10(ps)
    f, rPs = compute_hanning_power_spec(
        phantom_rf_window, config.transducer_freq_band[0],
        config.transducer_freq_band[1], config.sampling_frequency
    )
    rPs = 20 * np.log10(rPs)
    nps = np.asarray(ps) - np.asarray(rPs)
    # Fill in attributes defined in ResultsClass above
    window.results.nps = nps # dB
    window.results.f = f # Hz
    window.results.ps = ps # dB
    window.results.r_ps = rPs # dB

@output_vars("mbf", "ss", "si")
@dependencies("compute_power_spectra")
def lizzi_feleppa(scan_rf_window: np.ndarray, phantom_rf_window: np.ndarray, 
                    window: Window, config: RfAnalysisConfig, 
                    image_data: UltrasoundRfImage, **kwargs) -> None:
    """Compute spectral analysis values for a single window.
    
    Args:
        scan_rf_window (np.ndarray): RF data of the window in the scan image.
        phantom_rf_window (np.ndarray): RF data of the window in the phantom image.
        window (Window): Window object to store results.
        config (RfAnalysisConfig): Configuration object for analysis.
    """
    # Accessible from above function call
    nps = window.results.nps
    f = window.results.f
    
    def compute_spectral_params(nps: np.ndarray, f: np.ndarray, 
                               low_f: int, high_f: int) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
        """Perform spectral analysis on the normalized power spectrum.
        source: Lizzi et al. https://doi.org/10.1016/j.ultrasmedbio.2006.09.002
        
        Args:
            nps (np.ndarray): normalized power spectrum.
            f (np.ndarray): frequency array (Hz).
            low_f (int): lower bound of the frequency window for analysis (Hz).
            high_f (int): upper bound of the frequency window for analysis (Hz).
            
        Returns:
            Tuple: midband fit, frequency range, linear fit, and linear regression coefficients.
        """
        # 1. in one scan / run-through of data file's f array, find the data points on
        # the frequency axis closest to reference file's analysis window's LOWER bound and UPPER bounds
        smallest_diff_low_f = 999999999
        smallest_diff_high_f = 999999999

        for i in range(len(f)):
            current_diff_low_f = abs(low_f - f[i])
            current_diff_high_f = abs(high_f - f[i])

            if current_diff_low_f < smallest_diff_low_f:
                smallest_diff_low_f = current_diff_low_f
                smallest_diff_index_low_f = i

            if current_diff_high_f < smallest_diff_high_f:
                smallest_diff_high_f = current_diff_high_f
                smallest_diff_index_high_f = i

        # 2. compute linear regression within the analysis window
        f
        f_band = f[
            smallest_diff_index_low_f:smallest_diff_index_high_f
        ]  # transpose row vector f in order for it to have same dimensions as column vector nps
        p = np.polyfit(
            f_band, nps[smallest_diff_index_low_f:smallest_diff_index_high_f], 1
        )
        nps_linfit = np.polyval(p, f_band)  # y_linfit is a column vecotr

        mbfit = p[0] * f_band[round(f_band.shape[0] / 2)] + p[1]

        return mbfit, f_band, nps_linfit, p
    
    mbf, _, _, p = compute_spectral_params(nps, f, config.analysis_freq_band[0], config.analysis_freq_band[1])
    
    # Fill in attributes defined in "output_vars" decorator
    window.results.mbf = mbf # dB
    window.results.ss = p[0]*1e6 # dB/MHz
    window.results.si = p[1] # dB
    
@output_vars("att_coef")
@required_kwargs("ref_attenuation")
@dependencies("compute_power_spectra")
def attenuation_coef(scan_rf_window: np.ndarray, phantom_rf_window: np.ndarray, 
                    window: Window, config: RfAnalysisConfig, 
                    image_data: UltrasoundRfImage, **kwargs) -> None:
    """Compute the local attenuation coefficient of the ROI using the Spectral Difference
    Method for Local Attenuation Estimation. This method computes the attenuation coefficient
    for multiple frequencies and returns the slope of the attenuation as a function of frequency.

    Args:
        rf_data (np.ndarray): RF data of the ROI (n lines x m samples).
        ref_rf_data (np.ndarray): RF data of the phantom (n lines x m samples).
        overlap (float): Overlap percentage for analysis windows.
        window_depth (int): Depth of each window in samples.
    Updated and verified : Feb 2025 - IR
    """
    overlap = 50
    window_depth = min(100, scan_rf_window.shape[0] // 3)
    ref_attenuation = kwargs['ref_attenuation']  # Reference attenuation coefficient (dB/cm/MHz)
    
    sampling_frequency = config.sampling_frequency
    start_frequency = config.analysis_freq_band[0]
    end_frequency = config.analysis_freq_band[1]

    # Initialize arrays for storing intensities (log of power spectrum for each frequency)
    ps_sample = []  # ROI power spectra
    ps_ref = []     # Phantom power spectra

    start_idx = 0
    end_idx = window_depth
    window_center_indices = []
    counter = 0

    # Loop through the windows in the RF data
    while end_idx < scan_rf_window.shape[0]:
        sub_window_rf = scan_rf_window[start_idx:end_idx]
        f, ps = compute_hanning_power_spec(sub_window_rf, start_frequency, end_frequency, sampling_frequency)
        ps_sample.append(20 * np.log10(ps))  # Log scale intensity for the ROI

        ref_sub_window_rf = phantom_rf_window[start_idx:end_idx]
        ref_f, ref_ps = compute_hanning_power_spec(ref_sub_window_rf, start_frequency, end_frequency, sampling_frequency)
        ps_ref.append(20 * np.log10(ref_ps))  # Log scale intensity for the phantom

        window_center_indices.append((start_idx + end_idx) // 2)

        start_idx += int(window_depth * (1 - (overlap / 100)))
        end_idx = start_idx + window_depth
        counter += 1

    # Convert window depths to cm
    axial_res_cm = image_data.axial_res / 10
    window_depths_cm = np.array(window_center_indices) * axial_res_cm

    attenuation_coefficients = []  # One coefficient for each frequency

    f = f / 1e6
    ps_sample = np.array(ps_sample)
    ps_ref = np.array(ps_ref)

    mid_idx = f.shape[0] // 2
    start_idx = max(0, mid_idx - 25)
    end_idx = min(f.shape[0], mid_idx + 25)

    # Compute attenuation for each frequency
    for f_idx in range(start_idx, end_idx):
        normalized_intensities = np.subtract(ps_sample[:, f_idx], ps_ref[:, f_idx])
        p = np.polyfit(window_depths_cm, normalized_intensities, 1)
        local_attenuation = ref_attenuation * f[f_idx] - (1 / 4) * p[0]  # dB/cm
        attenuation_coefficients.append(local_attenuation / f[f_idx])  # dB/cm/MHz

    attenuation_coef = np.mean(attenuation_coefficients)
    
    # Fill in attributes defined in "output_vars" decorator
    window.results.att_coef = attenuation_coef # dB/cm/MHz
    
@output_vars("bsc")
@required_kwargs("ref_bsc", "ref_attenuation")
@dependencies("compute_power_spectra", "attenuation_coef")
def bsc(scan_rf_window: np.ndarray, phantom_rf_window: np.ndarray, 
            window: Window, config: RfAnalysisConfig, 
            image_data: UltrasoundRfImage, **kwargs) -> None:
    """Compute the backscatter coefficient of the ROI using the reference phantom method.
    Assumes instrumentation and beam terms have the same effect on the signal from both 
    image and phantom. 

    Source: Yao et al. (1990): https://doi.org/10.1177/016173469001200105. PMID: 2184569

    Args:
        freq_arr (np.ndarray): Frequency array of power spectra (Hz).
        scan_ps (np.ndarray): Power spectrum of the analyzed scan at the current region.
        ref_ps (np.ndarray): Power spectrum of the reference phantom at the current region.
        att_coef (float): Attenuation coefficient of the current region (dB/cm/MHz).
        frequency (int): Frequency on which to compute backscatter coefficient (Hz).
        roi_depth (int): Depth of the start of the ROI in samples.
        
    Returns:
        float: Backscatter coefficient of the ROI for the central frequency (1/cm-sr).
        Updated and verified : Feb 2025 - IR
    """
    freq_arr = window.results.f
    scan_ps = window.results.ps
    ref_ps = window.results.r_ps
    att_coef = window.results.att_coef
    roi_depth = scan_rf_window.shape[0]
    
    # Required kwargs
    ref_attenuation = kwargs['ref_attenuation']
    ref_backscatter_coef = kwargs['ref_bsc']
    
    # Optional kwarg
    frequency = kwargs.get('bsc_freq', config.center_frequency)  # Frequency for backscatter coefficient calculation
    
    index = np.argmin(np.abs(freq_arr - frequency))
    ps_sample = scan_ps[index]
    ps_ref = ref_ps[index]
    s_ratio = ps_sample / ps_ref

    np_conversion_factor = np.log(10) / 20 
    converted_att_coef = att_coef * np_conversion_factor  # dB/cm/MHz -> Np/cm/MHz
    converted_ref_att_coef = ref_attenuation * np_conversion_factor  # dB/cm/MHz -> Np/cm/MHz

    window_depth_cm = roi_depth * image_data.axial_res / 10  # cm
    converted_att_coef *= frequency / 1e6  # Np/cm
    converted_ref_att_coef *= frequency / 1e6  # Np/cm        

    att_comp = np.exp(4 * window_depth_cm * (converted_att_coef - converted_ref_att_coef)) 
    bsc = s_ratio * ref_backscatter_coef * att_comp

    # Fill in attributes defined in "output_vars" decorator
    window.results.bsc = bsc # 1/cm-sr

@output_vars("nak_w", "nak_u")
def nakagami_params(scan_rf_window: np.ndarray, phantom_rf_window: np.ndarray, 
            window: Window, config: RfAnalysisConfig, 
            image_data: UltrasoundRfImage, **kwargs) -> None:
    """Compute Nakagami parameters for the ROI.

    Source: Tsui, P. H., Wan, Y. L., Huang, C. C. & Wang, M. C. 
    Effect of adaptive threshold filtering on ultrasonic Nakagami 
    parameter to detect variation in scatterer concentration. Ultrason. 
    Imaging 32, 229–242 (2010). https://doi.org/10.1177%2F016173461003200403

    Args:
        rf_data (np.ndarray): RF data of the ROI (n lines x m samples).
        
    Returns:
        Tuple: Nakagami parameters (w, u) for the ROI.
    """
    r = np.abs(hilbert(scan_rf_window, axis=1))
    w = np.nanmean(r ** 2, axis=1)
    u = (w ** 2) / np.var(r ** 2, axis=1)

    # Averaging to get single parameter values
    w = np.nanmean(w)
    u = np.nanmean(u)
    
    # Fill in attributes defined in "output_vars" decorator
    window.results.nak_w = w
    window.results.nak_u = u
    
@output_vars("kappa", "mu", "hk", "s", "omega", "alpha")
def hkd_params(scan_rf_window: np.ndarray, phantom_rf_window: np.ndarray, 
            window: Window, config: RfAnalysisConfig, 
            image_data: UltrasoundRfImage, **kwargs) -> None:
    """
    Implemented: IR, April 2025
    Compute Homodyned K distribution parameters for the ROI.
    Based on: 
        [1] Destrempes, F., et al. (2013). "Estimation Method of the Homodyned K-Distribution Based on the Mean Intensity and Two Log-Moments." SIAM J Imaging Sci 6(3): 1499-1530.
        [2] Cristea, A., et al. (2020). "Quantitative assessment of media concentration using the Homodyned K distribution." Ultrasonics 101: 105986.
        Args:
    env: envelope of RF data
    Returns:
    alpha: scatterer clustering parameter
    kappa: structure parameter
    *****************************************************
    Generally following notation in Destrempes 2013, the parameters are:
    alpha = scatterer clustering parameter
    kappa = structure parameter
    mu = mean intensity
    hk = ratio of coherent to diffuse signal
    epsilon = coherent signal
    sigma = related to incoherent signal
    omega = related to sigma
    """
    def simplif_hyper_term1(alpha, gama):
        value = 0.0
        oldvalue = value
        k = 0
        vect = np.arange(1, k + 1) - alpha
        prod_vect = np.prod(vect) if len(vect) > 0 else 1
        value += (1 / ((k + 1) * factorial(k + 1))) * (1 / prod_vect) * (gama ** k)
        while abs(value - oldvalue) > 0.001:
            oldvalue = value
            k += 1
            vect = np.arange(1, k + 1) - alpha
            prod_vect = np.prod(vect)
            value += (1 / ((k + 1) * factorial(k + 1))) * (1 / prod_vect) * (gama ** k)
        return value

    def simplif_hyper_term2(alpha, gama):
        value = 0.0
        oldvalue = value
        k = 0
        vect = np.append(np.arange(1, k + 1), k) + alpha
        prod_vect = np.prod(vect) if len(vect) > 0 else 1
        value += alpha / factorial(k) * (1 / prod_vect) * (gama ** k)
        while abs(value - oldvalue) > 0.001:
            oldvalue = value
            k += 1
            vect = np.append(np.arange(1, k + 1), k) + alpha
            prod_vect = np.prod(vect)
            value += alpha / factorial(k) * (1 / prod_vect) * (gama ** k)
        return value

    def simplif_hyper_term3(alpha, gama):
        value = 0.001
        oldvalue = value
        k = 0
        vect = np.append(np.arange(2, k + 2), k + 1) + alpha
        prod_vect = np.prod(vect) if len(vect) > 0 else 1
        value += (alpha + 1) / factorial(k) * (1 / prod_vect) * (gama ** k)
        while abs(value - oldvalue) > 0.001:
            oldvalue = value
            k += 1
            vect = np.append(np.arange(2, k + 2), k + 1) + alpha
            prod_vect = np.prod(vect)
            value += (alpha + 1) / factorial(k) * (1 / prod_vect) * (gama ** k)
        return value

    def simplif_hyper_term4(alpha, gama):
        value = 0.0
        oldvalue = value
        k = 0
        vect = np.arange(2, k + 2) - alpha
        prod_vect = np.prod(vect) if vect.size > 0 else 1
        value += 1 / ((k + 1) * factorial(k + 1)) * (1 / prod_vect) * (gama ** k)
        while abs(value - oldvalue) > 0.001:
            oldvalue = value
            k += 1
            vect = np.arange(2, k + 2) - alpha
            prod_vect = np.prod(vect) if vect.size > 0 else 1
            value += 1 / ((k + 1) * factorial(k + 1)) * (1 / prod_vect) * (gama ** k)
        return value

    def UHK(gama, alpha):
        eulergamma = 0.5772156649
        term1 = -eulergamma - np.log(gama + alpha) + psi(alpha)
        term2 = (np.pi / np.sin(np.pi * alpha) * gama ** alpha *
                    simplif_hyper_term2(alpha, gama) /
                    (alpha * gamma(alpha) * gamma(alpha + 1)))
        if alpha < 1:
            term3 = (-np.pi / np.sin(np.pi * (1 - alpha)) *
                        gama / (gamma(alpha) * gamma(2 - alpha)) *
                        simplif_hyper_term4(alpha, gama))
        else:
            term3 = gama / (alpha - 1) * simplif_hyper_term4(alpha, gama)
        return term1 + term2 + term3

    def UHK_loop(gama, alpha):
        if round(alpha) == alpha:
            Uhk1 = UHK(gama, alpha - 1e-7)
            Uhk2 = UHK(gama, alpha + 1e-7)
            return Uhk1 + (Uhk2 - Uhk1) / 2
        else:
            return UHK(gama, alpha)

    def XHK(gama, alpha):
        term0 = ((1 + 2 * alpha) / (gama + alpha) -
                    2 * gama ** (alpha / 2 + 0.5) /
                    ((gama + alpha) * gamma(alpha)) *
                    kv(alpha + 1, 2 * np.sqrt(gama)))
        term1 = simplif_hyper_term1(alpha, gama)
        term2 = (-np.pi / np.sin(np.pi * alpha) *
                    gama ** (alpha - 1) * simplif_hyper_term2(alpha, gama) /
                    (gamma(alpha) * gamma(alpha + 1)))
        term3 = (np.pi / np.sin(np.pi * (alpha + 1)) *
                    gama ** alpha * simplif_hyper_term3(alpha, gama) /
                    (gamma(alpha) * gamma(alpha + 2) * (1 + alpha)))
        if alpha < 1:
            term4 = (np.pi / np.sin(np.pi * (1 - alpha)) *
                        alpha / (gamma(alpha) * gamma(2 - alpha)) *
                        simplif_hyper_term4(alpha, gama))
        else:
            term4 = -alpha / (alpha - 1) * simplif_hyper_term4(alpha, gama)
        return term0 + gama / (gama + alpha) * (term1 + term2 + term3 + term4)

    def XHK_loop(gama, alpha):
        if round(alpha) == alpha:
            Xhk1 = XHK(gama, alpha - 1e-7)
            Xhk2 = XHK(gama, alpha + 1e-7)
            return Xhk1 + (Xhk2 - Xhk1) / 2
        else:
            return XHK(gama, alpha)

    def compute_gama(alpha: float, X: float) -> float:
        gama1 = 0
        gama2 = 1
        while XHK_loop(gama2, alpha) > X:
            gama2 += 1
        tol = 1e-4
        while abs(gama1 - gama2) > tol:
            gama = (gama1 + gama2) / 2
            if XHK_loop(gama, alpha) <= X:
                gama2 = gama
            else:
                gama1 = gama
        return gama

    def compute_gama_alpha(X, U):
        alpha0 = 1 / (X - 1)
        alpha1 = 0
        if X <= 1:
            alpha2 = 1
            gama = compute_gama(alpha2, X)
            while UHK_loop(gama, alpha2) > U:
                alpha2 += 1
                gama = compute_gama(alpha2, X)
        else:
            alpha2 = min(alpha0, 80)
        tol = 1e-2
        while abs(alpha1 - alpha2) > tol:
            alpha = (alpha1 + alpha2) / 2
            gama = compute_gama(alpha, X)
            if UHK_loop(gama, alpha) <= U:
                alpha2 = alpha
            else:
                alpha1 = alpha
        return gama, alpha

    # --- Envelope detection and moment statistics ---
    env = np.abs(hilbert(scan_rf_window, axis=1))
    #env = rfData()
    data = env.flatten()
    I = data ** 2
    I = I.astype(float)
    I[I == 0] = np.finfo(float).eps

    MEAN = np.mean(I)
    X = np.mean(I * np.log(I)) / MEAN - np.mean(np.log(I))
    U = np.mean(np.log(I)) - np.log(MEAN)

    if X <= 1:
        alpha_star = 80.0
    else:
        alpha0 = 1.0 / (X - 1)
        alpha_star = min(alpha0, 80.0)

    computed_gamma = compute_gama(alpha_star, X)

    if UHK_loop(computed_gamma, alpha_star) <= U:
        computed_gamma, alpha = compute_gama_alpha(X, U)
    else:
        alpha = alpha_star

    epsilon = np.sqrt(MEAN * computed_gamma / (computed_gamma + alpha))
    sigma = np.sqrt(MEAN / (2 * (computed_gamma + alpha)))
    kappa = epsilon**2 / (2 * sigma**2 * alpha)
    mu = np.mean(env**2)
    hk = epsilon / (sigma * np.sqrt(alpha))
    omega = mu / (hk**2 + 2)
    s = hk * np.sqrt(omega)
    
    window.results.kappa = kappa
    window.results.mu = mu
    window.results.hk = hk
    window.results.s = s
    window.results.omega = omega
    window.results.alpha = alpha



# experimenting with adding functions
@output_vars("mean_rf_amplitude")
def compute_mean_rf(scan_rf_window: np.ndarray, phantom_rf_window: np.ndarray, 
                    window: Window, config: RfAnalysisConfig, 
                    image_data: UltrasoundRfImage, **kwargs) -> None:
    """
    Compute the mean amplitude of the RF window.

    Args:
        scan_rf_window (np.ndarray): RF data of the window in the scan image.
        window (Window): Window object to store results.
        seg_data: Unused placeholder for compatibility.
        analysis_meta: Unused placeholder for compatibility.
        analysis_config: Unused placeholder for compatibility.
    """
    mean_val = np.mean(np.abs(scan_rf_window))
    window.results.mean_rf_amplitude = mean_val