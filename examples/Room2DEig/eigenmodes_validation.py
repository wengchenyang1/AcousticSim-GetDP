import numpy as np


def compute_room_acoustic_modes(Lx, Ly, m_values=[0, 1, 2, 3], n_values=[0, 1, 2, 3]):
    modes = []
    for m in m_values:
        for n in n_values:
            k_squared = (m * np.pi / Lx) ** 2 + (n * np.pi / Ly) ** 2
            modes.append(np.sqrt(k_squared))
    return modes


# Define the dimensions of the room
Lx = 10.0  # Length in x-direction
Ly = 8.0  # Length in y-direction

# Compute room acoustic modes
modes = compute_room_acoustic_modes(Lx, Ly)

# Print the computed modes
c0 = 343
print("Room Acoustic Modes:")
for i, mode in enumerate(modes):
    print(f"Mode {i+1} (k): {mode}, with f = {mode*c0/2/np.pi} Hz")
