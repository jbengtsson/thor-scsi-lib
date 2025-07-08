import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Biot-Savart kernel (regularized to avoid singularity)
def biot_savart(x, xp, w, delta):
    r = x - xp
    r2 = np.sum(r**2)
    if r2 < 1e-12:
        return np.zeros(3)
    factor = (np.cross(w, r)) / (r2 + delta**2)**(3/2)
    return factor / (4 * np.pi)

# Initialize vortex ring as a toroidal loop of vortex particles
def init_vortex_ring(R=1.0, a=0.1, N_ring=100):
    thetas = np.linspace(0, 2*np.pi, N_ring, endpoint=False)
    positions = np.zeros((N_ring, 3))
    vorticities = np.zeros((N_ring, 3))

    for i, theta in enumerate(thetas):
        x = R * np.cos(theta)
        y = R * np.sin(theta)
        z = 0
        positions[i] = [x, y, z]

        # Vorticity direction tangential to ring
        vx = -np.sin(theta)
        vy = np.cos(theta)
        vz = 0
        vorticities[i] = a * np.array([vx, vy, vz])  # a scales strength

    return positions, vorticities

# Parameters
N = 200            # number of particles in ring
delta = 0.05       # core smoothing radius
dt = 0.02          # time step
T = 2.0            # total time
steps = int(T/dt)

# Initialize vortex ring
pos, vort = init_vortex_ring(R=1.0, a=1.0, N_ring=N)

# Time loop
for step in range(steps):
    vel = np.zeros_like(pos)
    for i in range(N):
        for j in range(N):
            if i != j:
                vel[i] += biot_savart(pos[i], pos[j], vort[j], delta)

    pos += dt * vel  # Euler time integration

    # Visualization
    if step % 10 == 0:
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], 'bo-', markersize=3)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ax.set_zlim(-1, 1)
        ax.set_title(f"Step {step}, Time {step*dt:.2f}")
        plt.tight_layout()
        plt.pause(0.01)
        plt.clf()

plt.show()
