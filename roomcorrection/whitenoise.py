import pyaudio
import numpy as np

# Parameters
duration = 65  # Duration in seconds
sample_rate = 44100  # Sample rate in Hz
volume = 0.5  # Volume (0.0 to 1.0)

# Generate white noise
samples = (np.random.rand(int(duration * sample_rate)) * 2 - 1).astype(np.float32)

# Initialize PyAudio
p = pyaudio.PyAudio()

# Open stream
stream = p.open(format=pyaudio.paFloat32,
                channels=1,
                rate=sample_rate,
                output=True)

# Play white noise
stream.write(volume * samples)

# Stop and close the stream
stream.stop_stream()
stream.close()

# Terminate PyAudio
p.terminate()
