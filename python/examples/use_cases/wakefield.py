import os
import math
import numpy as np
from   numpy.polynomial import hermite_e
import matplotlib.pyplot as plt


def rd_csv_file(file_name, k):
    data = np.loadtxt(file_name, delimiter=" ", dtype=float)
    bunch = np.array([data[0,:], data[k, :]])
    return bunch


def compute_moment(n, bunch):
    # Time step.
    dt = bunch[0, 1] - bunch[0, 0]
    return np.sum(bunch[0, :]**n*bunch[1, :])*dt


def normalise_distr(bunch):
    # Normalise distribution.
    bunch[1] /= compute_moment(0, bunch)
    # Remove mean.
    bunch[0] -= compute_moment(1, bunch)
    # Normalise by sigma.
    bunch[0] /= np.sqrt(compute_moment(2, bunch))

    # Renormalise distribution.
    bunch[1] /= compute_moment(0, bunch)

    m_1, m_2 = compute_moment(1, bunch), compute_moment(2, bunch)
    print("\n[mean, sigma, sigma^2] =\n  [{:10.3e}, {:9.3e}, {:9.3e}]".
          format(m_1, np.sqrt(m_2), m_2))
    return bunch


def Gaussian(m_1, m_2, t):
   return np.exp(-(t-m_1)**2/(2e0*m_2))/np.sqrt(2e0*np.pi*m_2)


def compute_residual(m_1, m_2, bunch):
    bunch[1, :] /= Gaussian(m_1, m_2, bunch[0, :])
    return bunch


def wrt_csv_file(file_name, bunch):
    np.savetxt(file_name, bunch, fmt="%.18e", delimiter=" ", newline='\n')


def plot_Hermite():

    x = np.arange(-2.4, 3.5, 0.01)
    plt.plot(x, hermite.hermval(x, 1),                  c = "r",
             label = "n = 0")
    plt.plot(x, hermite.hermval(x, [0, 1]),             c = "g",
             label = "n = 1")
    plt.plot(x, hermite.hermval(x, [0, 0, 1]),          c = "b",
             label = "n = 2")
    plt.plot(x, hermite.hermval(x, [0, 0, 0, 1]),       c = "c",
             label = "n = 3")
    plt.plot(x, hermite.hermval(x, [0, 0, 0, 0, 1]),    c = "g",
             label = "n = 4")
    plt.plot(x, hermite.hermval(x, [0, 0, 0, 0, 0, 1]), c = "m",
             label = "n = 5")
    plt.legend()
    plt.title("Hermite Polynomials (physicists)")

    ax = plt.gca()
    ax.set_xlim([-2.4, 3.5])
    ax.set_ylim([-50, 50])

    plt.show()


def plot_bunch(title, bunch, bunch_fit, with_Gaussian):
    fig, gr = plt.subplots(1)
    gr.set_title(title)
    gr.set_xlabel("t [nsec]")
    gr.set_ylabel(r"$\rho$ []")
    gr.plot(bunch[0], bunch[1],  "blue", label="Bunch Profile")
    gr.plot(bunch[0], bunch_fit, "cyan", label="Fitted")
    gr.legend()

    fig.tight_layout()

    if not with_Gaussian:
        ax = plt.gca()
        ax.set_xlim([-3.5, 3.0])
        ax.set_ylim([-0.1, 2.0])
    plt.show()


def compute_c_n(n, bunch):
    # Time step.
    dt = bunch[0][1] - bunch[0][0]
    h_n = np.zeros(n+1)
    h_n[n] = 1
    coeff = \
        np.sum(hermite_e.hermeval(bunch[0, :], h_n)**2 \
               *Gaussian(0e0, 1e0, bunch[0, :]))*dt
    c_n = np.sum(bunch[1, :]*hermite_e.hermeval(bunch[0, :], h_n))*dt/coeff

    return c_n


if False:
    plot_Hermite()

file_dir = os.path.join(
    os.environ["HOME"], "/Volumes/Ext. HD (Dropbox)", "Dropbox",
    "collective-effects-tracking", "TO", "examples-profiles")
if True:
    file_name = os.path.join(file_dir, "profiles_complexFF_173bunches.csv")
    bunch_first  = rd_csv_file(file_name, 173)
    bunch_middle = rd_csv_file(file_name, 88)
    bunch_last   = rd_csv_file(file_name, 1)
else:
    file_name = os.path.join(file_dir, "profiles_complexFF_164bunches.csv")
    bunch_first  = rd_csv_file(file_name, 164)
    bunch_middle = rd_csv_file(file_name, 82)
    bunch_last   = rd_csv_file(file_name, 1)

bunch = bunch_first;

np.set_printoptions(precision=3, linewidth=132)

# Change time unit to [nsec].
bunch[0] = 1e9*bunch[0]
bunch = normalise_distr(bunch)

if not False:
    c_max = 6

    with_Gaussian = True

    c = np.zeros(c_max+1)
    for j in range(0, c_max+1):
        c[j] = compute_c_n(j, bunch)
    print("c =\n ", c)

    if not with_Gaussian:
        bunch = compute_residual(0e0, 1e0, bunch)

    bunch_fit = np.zeros(len(bunch[0]))
    bunch_fit[:] = hermite_e.hermeval(bunch[0, :], c)
    if with_Gaussian:
        bunch_fit[:] *= Gaussian(0e0, 1e0, bunch[0, :])

    plot_bunch("Bunch Profile - First", bunch, bunch_fit, with_Gaussian)


if False:
    if not True:
        c_max = 6
    else:
        c_max = 16

    bunch = compute_residual(0e0, 1e0, bunch)
    c = hermite_e.hermefit(bunch[0], bunch[1], c_max)
    print("c =\n ", c)

    bunch_fit = np.zeros(len(bunch[0]))
    bunch_fit[:] = hermite_e.hermeval(bunch[0, :], c)

    plot_bunch("Bunch Profile - First", bunch, bunch_fit, True)


if False:
    # Convert to polynomial.
    p = hermite_e.herme2poly(c)
    bunch_fit = np.zeros(len(bunch[0]))
    bunch_fit[:] = np.polynomial.polynomial.polyval(bunch[0, :], p)
    if with_Gaussian:
        bunch_fit[:] *= Gaussian(0e0, 1e0, bunch[0, :])

    plot_bunch("Bunch Profile - First", bunch, bunch_fit, with_Gaussian)
