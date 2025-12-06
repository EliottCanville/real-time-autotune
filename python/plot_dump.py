import numpy as np
from matplotlib import pyplot as plt
from scipy.io.wavfile import write
import sys

# Question 10: Adapted script to read and plot dump files

# Choose which file to process (can be changed via command line)
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    # Default filename
    filename = 'dump_input'

# Sampling rate - change if necessary
fs = 44100

print(f"Reading file: {filename}")

# Read binary file (double format)
try:
    x = np.fromfile(filename + '.bin', dtype=np.double)
    print(f"Loaded {len(x)} samples")
except FileNotFoundError:
    print(f"File {filename}.bin not found!")
    sys.exit(1)

if len(x) == 0:
    print("File is empty!")
    sys.exit(1)

# Write wav file (for audio signals only, not for f0)
if 'f0' not in filename:
    write(filename + '.wav', fs, x)
    print(f"Wrote {filename}.wav")

# Plot
plt.figure(figsize=(12, 6))

# For f0 data, plot vs time with appropriate labels
if 'f0' in filename:
    # Convert sample index to time (each f0 value per buffer)
    # Assuming buffer size of 512 samples
    buffer_size = 2048
    time = np.arange(len(x)) * buffer_size / fs
    plt.plot(time, x)
    plt.xlabel('Time (s)')
    plt.ylabel('Fundamental Frequency f0 (Hz)')
    plt.title('f0 estimation over time')
    plt.ylim([80, 500])  # Typical vocal range
else:
    # For audio signals, plot vs sample index or time
    time = np.arange(len(x)) / fs
    plt.plot(time, x)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title(f'Signal: {filename}')

plt.grid()
plt.tight_layout()
plt.show()

print("Plot displayed")
