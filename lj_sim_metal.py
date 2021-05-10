# Lennard-Jones fluid simulation

# Same as v2 but using metal units for argon


import os

import numpy as np


OUT_DIR = "data/raw/"

# Lennard-Jones parameters
EPS = 117.7 * 8.617e-5
SIG = 3.504
R_CUT = 5.0 * SIG

# Simulation params
L = 10 * SIG  # simulation box dimensions
nparts = 20
nframes = 1000


def wrapped_cutoff_dists(r, R):
    absdiff = np.abs(R - r)
    wrapped_absdiff = L - absdiff
    pbc_dists = np.linalg.norm(np.minimum(absdiff, wrapped_absdiff), axis=1)
    pbc_dists = pbc_dists[pbc_dists <= L/2]
    return pbc_dists


def wrap_diffs(diff):
    wrapped_diff = diff.copy()

    neg_to_wrap = (diff < 0) & (-diff > L/2)
    wrapped_diff[neg_to_wrap] += L

    pos_to_wrap = (diff >= 0) & (diff > L/2)
    wrapped_diff[pos_to_wrap] -= L

    return wrapped_diff


def ljpotential(coords):
    u = 0
    for i, r in enumerate(coords):
        pbc_dists = wrapped_cutoff_dists(r, coords[i+1:])
        u += (4*EPS * ((SIG / pbc_dists)**12 - (SIG / pbc_dists)**6)).sum()
    return u


def ljforces(coords):
    forces = np.zeros(coords.shape)
    for i, pos in enumerate(coords):
        mask = np.ones(coords.shape[0], dtype=bool)
        mask[i] = False
        others = coords[mask]
        r = wrap_diffs(others - pos)
        r_mag = np.linalg.norm(r, axis=1).reshape((-1, 1))
        r_mag[r_mag >= R_CUT] = np.inf
        forces[i] = np.sum(-24*EPS/SIG * (2 * (SIG/r_mag)**13 - (SIG/r_mag)**7) * r/r_mag, axis=0)
    return forces


frame_coords = []
frame_energies = []
frame_forces = []
frame_boxes = []

for i in range(nframes):
    energy = 1
    while energy > 0:
        coords = np.random.uniform(0, L, size=(nparts, 3))
        energy = ljpotential(coords)
    frame_coords.append(np.concatenate(coords))
    frame_energies.append(energy)
    frame_forces.append(np.concatenate(ljforces(coords)))
    frame_boxes.append([L, 0, 0, 0, L, 0, 0, 0, L])
    if i % 100 == 99:
        print(f"Finished {i+1}")

np.savetxt(os.path.join(OUT_DIR, "box.raw"), np.array(frame_boxes), fmt="%.6f")
np.savetxt(os.path.join(OUT_DIR, "coord.raw"), np.array(frame_coords), fmt="%.6f")
np.savetxt(os.path.join(OUT_DIR, "energy.raw"), np.array(frame_energies), fmt="%.6f")
np.savetxt(os.path.join(OUT_DIR, "force.raw"), np.array(frame_forces), fmt="%.6f")

with open(os.path.join(OUT_DIR, "type.raw"), "w") as outfile:
    outfile.write(" ".join(["0"]*nparts))
