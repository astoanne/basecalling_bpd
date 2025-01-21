import streamlit as st
from matplotlib import pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

def analyze_signal_gradients(original_signal, start, end, sigma_value=0.75, threshold_factor=0.5):
    # Slice the signal according to the start and end indices
    raw_signal = original_signal[start:end]
    
    # Step 1: Apply Gaussian smoothing
    smoothed_signal = gaussian_filter1d(raw_signal, sigma=sigma_value)

    # Step 2: Compute the gradient of the smoothed signal and its absolute value
    gradient = np.gradient(smoothed_signal)
    absolute_gradient = np.abs(gradient)

    # Step 3: Define a threshold and identify significant gradients
    threshold = np.mean(absolute_gradient) + threshold_factor * np.std(absolute_gradient)
    significant_gradients = np.where(absolute_gradient > threshold)[0]

    # Additional: Using peak finding to refine the selection of gradient points
    peaks, _ = find_peaks(absolute_gradient, height=threshold)
    
    with st.expander("Intermediate Steps: Absolute Gradient and Peak Detection"):
        # Plotting each step using Streamlit
        st.subheader("Raw Signal")
        fig1, ax1 = plt.subplots(figsize=(10, 5))
        ax1.plot(raw_signal, label='Raw Signal')
        ax1.set_title('Raw Signal')
        ax1.set_xlabel('Sample Index')
        ax1.set_ylabel('Signal Amplitude')
        ax1.legend()
        st.pyplot(fig1)

        st.subheader("Smoothed Signal with Identified Points")
        fig2, ax2 = plt.subplots(figsize=(10, 5))
        ax2.plot(smoothed_signal, label='Smoothed Signal', color='orange')
        ax2.plot(significant_gradients, smoothed_signal[significant_gradients], 'rx', label='Identified Points')
        ax2.set_title('Smoothed Signal with Identified Points')
        ax2.set_xlabel('Sample Index')
        ax2.set_ylabel('Signal Amplitude')
        ax2.legend()
        st.pyplot(fig2)

        st.subheader("Absolute Gradient with Peaks")
        fig3, ax3 = plt.subplots(figsize=(10, 5))
        ax3.plot(absolute_gradient, label='Absolute Gradient', color='grey')
        ax3.plot(peaks, absolute_gradient[peaks], 'ro', label='Peaks in Gradient')
        ax3.set_title('Absolute Smoothed Signal Gradient')
        ax3.set_xlabel('Sample Index')
        ax3.set_ylabel('Gradient Amplitude')
        ax3.legend()
        st.pyplot(fig3)

    # Adjust peaks to match original signal indices
    adjusted_peaks = peaks + start

    return adjusted_peaks
def analyze_signal_gradients_no_plot(original_signal, start, end, sigma_value=0.75, threshold_factor=0.5):
    # Slice the signal according to the start and end indices
    raw_signal = original_signal[start:end]
    
    # Step 1: Apply Gaussian smoothing
    smoothed_signal = gaussian_filter1d(raw_signal, sigma=sigma_value)

    # Step 2: Compute the gradient of the smoothed signal and its absolute value
    gradient = np.gradient(smoothed_signal)
    absolute_gradient = np.abs(gradient)

    # Step 3: Define a threshold and identify significant gradients
    threshold = np.mean(absolute_gradient) + threshold_factor * np.std(absolute_gradient)

    # Additional: Using peak finding to refine the selection of gradient points
    peaks, _ = find_peaks(absolute_gradient, height=threshold)

    # Adjust peaks to match original signal indices
    adjusted_peaks = peaks + start

    return adjusted_peaks