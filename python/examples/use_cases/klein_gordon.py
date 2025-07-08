import numpy as np
import matplotlib.pyplot as plt

# Parameters
r_max = 50.0
N = 500
dr = r_max / N
dt = 0.01
c = 1.0  # speed of light (units)
t_max = 20.0

# Radial grid
r = np.linspace(dr, r_max, N)  # avoid r=0 singularity

# Initial condition: Gaussian pulse
phi = np.exp(-((r - 10)**2) / 4)
phi_old = phi.copy()
phi_new = np.zeros_like(phi)

# Nonlinear potential derivative: V'(phi) = dV/dphi (e.g., quartic)
def V_prime(phi):
    return phi**3

# Time evolution loop using leapfrog method
t = 0.0
while t < t_max:
    # Apply finite difference for second spatial derivative (1D radial Laplacian with spherical symmetry)
    laplacian = (np.roll(phi, -1) - 2*phi + np.roll(phi, 1)) / dr**2
    laplacian[0] = (phi[1] - 2*phi[0]) / dr**2  # boundary at r=dr

    # Radial term: (2/r) d_phi/dr
    dphi_dr = (np.roll(phi, -1) - np.roll(phi, 1)) / (2*dr)
    radial_term = (2 / r) * dphi_dr

    # Wave eq: d2phi/dt2 = c^2 (laplacian + radial_term) - V'(phi)
    phi_new = 2*phi - phi_old + dt**2 * (c**2 * (laplacian + radial_term) - V_prime(phi))

    # Update steps
    phi_old = phi.copy()
    phi = phi_new.copy()

    t += dt

# Plot final profile
plt.plot(r, phi)
plt.xlabel('r')
plt.ylabel('phi(r,t_max)')
plt.title('Nonlinear scalar field profile at t = {:.2f}'.format(t_max))
plt.show()
