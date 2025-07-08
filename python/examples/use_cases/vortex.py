import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, ifft2, fftfreq

# Simulation parameters
N = 128              # Grid size
L = 2 * np.pi        # Domain size
dt = 0.01            # Time step
nt = 500             # Number of time steps
vis_plot = 50        # Plot every n steps

x = np.linspace(0, L, N, endpoint=False)
y = np.linspace(0, L, N, endpoint=False)
X, Y = np.meshgrid(x, y, indexing='ij')

dx = dy = L / N

# Initial vorticity: Gaussian vortex
def init_vortex(x0, y0, strength, width):
    return strength * np.exp(-((X - x0)**2 + (Y - y0)**2) / (2 * width**2))

omega = init_vortex(L/2, L/2, 5.0, 0.2)

# Wave numbers for FFT
kx = fftfreq(N, d=dx) * 2 * np.pi
ky = fftfreq(N, d=dy) * 2 * np.pi
KX, KY = np.meshgrid(kx, ky, indexing='ij')
K2 = KX**2 + KY**2
K2[0, 0] = 1.0  # Avoid division by zero

# Time stepping loop
for t in range(nt):
    # Solve Poisson equation: ∇²ψ = -ω → ψ̂ = -ω̂ / k²
    omega_hat = fft2(omega)
    psi_hat = -omega_hat / K2
    psi = np.real(ifft2(psi_hat))

    # Compute velocity: u = -dψ/dy, v = dψ/dx
    u = -np.real(ifft2(1j * KY * psi_hat))
    v =  np.real(ifft2(1j * KX * psi_hat))

    # Advect vorticity using 2nd-order upwind scheme
    dωdx = (np.roll(omega, -1, axis=0) - np.roll(omega, 1, axis=0)) / (2 * dx)
    dωdy = (np.roll(omega, -1, axis=1) - np.roll(omega, 1, axis=1)) / (2 * dy)
    omega -= dt * (u * dωdx + v * dωdy)

    # Visualization
    if t % vis_plot == 0:
        plt.clf()
        plt.contourf(X, Y, omega, levels=50, cmap='RdBu')
        plt.colorbar(label='Vorticity')
        plt.title(f"Time = {t * dt:.2f}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.axis('equal')
        plt.pause(0.01)

plt.show()
