import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from scipy.interpolate import splprep, splev

# Biot–Savart kernel (regularized)
def biot_savart(x, xp, w, delta):
    r = x - xp
    r2 = np.sum(r**2)
    if r2 < 1e-12:
        return np.zeros(3)
    factor = np.cross(w, r) / (r2 + delta**2)**1.5
    return factor / (4 * np.pi)

# Initialize vortex ring
def init_vortex_ring(R=1.0, a=0.1, N_ring=300):
    thetas = np.linspace(0, 2 * np.pi, N_ring, endpoint=False)
    positions = np.zeros((N_ring, 3))
    vorticities = np.zeros((N_ring, 3))
    for i, theta in enumerate(thetas):
        x = R * np.cos(theta)
        y = R * np.sin(theta)
        z = 0
        positions[i] = [x, y, z]
        vx = -np.sin(theta)
        vy = np.cos(theta)
        vorticities[i] = a * np.array([vx, vy, 0])
    return positions, vorticities

# Simulation parameters
N = 300
delta = 0.05
dt = 0.01
T = 2.0
steps = int(T / dt)
pos, vort = init_vortex_ring(R=1.0, a=1.0, N_ring=N)

# Prepare figure
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
line, = ax.plot([], [], [], 'b-', linewidth=2)

# Set up axis limits (adjust dynamically in animation)
def init():
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_zlim(-2, 2)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("3D Vortex Ring")
    return line,

# Animation step
def update(frame):
    global pos
    vel = np.zeros_like(pos)
    for i in range(N):
        for j in range(N):
            if i != j:
                vel[i] += biot_savart(pos[i], pos[j], vort[j], delta)
    pos += dt * vel

    try:
        tck, _ = splprep(pos.T, s=1.0, per=True)
        t_fine = np.linspace(0, 1, 1000)
        smooth_pos = np.array(splev(t_fine, tck)).T
        line.set_data(smooth_pos[:, 0], smooth_pos[:, 1])
        line.set_3d_properties(smooth_pos[:, 2])
    except:
        line.set_data(pos[:, 0], pos[:, 1])
        line.set_3d_properties(pos[:, 2])

    # Centering
    center = np.mean(pos, axis=0)
    r = 1.5
    ax.set_xlim(center[0] - r, center[0] + r)
    ax.set_ylim(center[1] - r, center[1] + r)
    ax.set_zlim(center[2] - r, center[2] + r)

    return line,

# Create animation
ani = FuncAnimation(fig, update, frames=steps, init_func=init, blit=False, interval=30, repeat=False)

plt.tight_layout()
plt.show()
