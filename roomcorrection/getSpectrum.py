import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import welch


# Read the eigen frequencies from the file
with open('roomcorrection/eigfreq.txt', 'r') as f:
    eigen_frequencies = [float(line.strip()) for line in f]

# Read the wav file
sample_rate, data = wavfile.read('roomcorrection/whitenoise.wav')

# Calculate the segment length
nperseg = int(sample_rate / 3)

# Perform pwelch to calculate the power spectral density with the new segment length
frequencies, psd = welch(data, sample_rate, nperseg=nperseg)

# Plot the power spectral density
plt.figure(figsize=(10, 6))
plt.semilogy(frequencies, psd)
plt.xlim((0, 200))
plt.title('Power Spectral Density of White Noise')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power/Frequency [dB/Hz]')
plt.grid()

# Plot vertical lines for each eigen frequency
for freq in eigen_frequencies:
    plt.axvline(x=freq, color='r', linestyle='--')

plt.show()
