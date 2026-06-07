#!/usr/bin/env python3
"""
current_strip_poisson_solver.py

Finite-difference Poisson solver and least-squares optimizer for local
current-strip nonlinear corrections in insertion devices.

Model
-----
For long, straight current strips carrying current along the local s/z axis,
solve the 2-D magnetostatic Poisson equation

    (d_x^2 + d_y^2) A_s(x,y) = -mu0 * J_s(x,y)

on a transverse rectangular box with Dirichlet boundary A_s = 0.
The box should be chosen large compared with the beam aperture / strip array.

Right-handed field convention used internally:

    B_x =  dA_s/dy
    B_y = -dA_s/dx

For a local left-handed accelerator coordinate convention, flip signs in
`fields_from_A` or in the target first integrals consistently.

The current-strip correction is then obtained by building the linear response
matrix Q from unit currents in the strips and solving

    Q @ I_strip ~= target_first_integral

using least squares, matching the approach described in the NSLS-II current-strip
papers.

Dependencies
------------
numpy, scipy, matplotlib

Typical use
-----------
1. Replace the synthetic `target_Iy_Tm` in `main()` by a target first-integral
   curve derived from the ID second-order kick map.
2. Adjust strip geometry, gap, length, and domain.
3. Run:

       python current_strip_poisson_solver.py

Outputs:
    strip_currents.csv
    current_strip_fit.png

Units
-----
SI throughout. A first field integral is in T*m.
1 G*cm = 1e-6 T*m.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
from numpy.typing import NDArray

try:
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    from scipy.interpolate import RegularGridInterpolator
    from scipy.optimize import lsq_linear
except ImportError as exc:
    raise SystemExit(
        "This script requires scipy. Install with: pip install scipy"
    ) from exc

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "This script requires matplotlib. Install with: pip install matplotlib"
    ) from exc


MU0 = 4.0e-7 * np.pi


@dataclass(frozen=True)
class Strip:
    """Rectangular current strip cross-section.

    x0, y0      center position [m]
    width       horizontal strip width [m]
    thickness   vertical strip thickness [m]
    length      effective magnetic length [m]
    name        optional label
    """
    x0: float
    y0: float
    width: float
    thickness: float
    length: float
    name: str = ""


def gcm_to_t_m(value_gcm: NDArray[np.float64] | float) -> NDArray[np.float64] | float:
    """Convert gauss*cm to tesla*m."""
    return np.asarray(value_gcm) * 1.0e-6


def t_m_to_gcm(value_tm: NDArray[np.float64] | float) -> NDArray[np.float64] | float:
    """Convert tesla*m to gauss*cm."""
    return np.asarray(value_tm) * 1.0e6


def make_two_strip_rows(
    *,
    n_per_row: int = 20,
    pitch: float = 3.0e-3,
    width: float = 2.0e-3,
    thickness: float = 0.3e-3,
    gap: float = 10.7e-3,
    length: float = 6.5,
) -> list[Strip]:
    """Create symmetric upper/lower strip rows.

    The defaults are close to the published current-strip dimensions:
    2 mm wide, 0.3 mm thick, 1 mm spacing, 6.5 m long, with an example
    upper/lower gap of 10.7 mm.

    Positive current is along +s/+z for both rows. The fitted currents may
    naturally choose opposite signs for the two rows.
    """
    x_centers = (np.arange(n_per_row) - 0.5 * (n_per_row - 1)) * pitch
    strips: list[Strip] = []
    for row_name, y0 in (("top", +gap / 2.0), ("bottom", -gap / 2.0)):
        for i, x0 in enumerate(x_centers):
            strips.append(
                Strip(
                    x0=float(x0),
                    y0=float(y0),
                    width=width,
                    thickness=thickness,
                    length=length,
                    name=f"{row_name}_{i:02d}",
                )
            )
    return strips


class PoissonStripSolver:
    """Finite-difference Poisson solver for transverse current strips."""

    def __init__(
        self,
        *,
        xlim: tuple[float, float] = (-0.070, 0.070),
        ylim: tuple[float, float] = (-0.035, 0.035),
        nx: int = 241,
        ny: int = 141,
        mu0: float = MU0,
    ) -> None:
        if nx < 5 or ny < 5:
            raise ValueError("nx and ny must be at least 5.")
        self.x = np.linspace(xlim[0], xlim[1], nx)
        self.y = np.linspace(ylim[0], ylim[1], ny)
        self.dx = float(self.x[1] - self.x[0])
        self.dy = float(self.y[1] - self.y[0])
        self.nx = nx
        self.ny = ny
        self.mu0 = mu0

        self._lap = self._build_laplacian_dirichlet()
        # Factorization is reused for every unit-current strip.
        self._solve = spla.factorized(self._lap.tocsc())

    def _build_laplacian_dirichlet(self) -> sp.csr_matrix:
        """Sparse Laplacian on interior points with A=0 boundary."""
        nxi = self.nx - 2
        nyi = self.ny - 2
        n = nxi * nyi

        main = np.full(n, -2.0 / self.dx**2 - 2.0 / self.dy**2)
        east_west = np.full(n - 1, 1.0 / self.dx**2)
        north_south = np.full(n - nxi, 1.0 / self.dy**2)

        # Break east/west coupling at row boundaries.
        for j in range(1, nyi):
            east_west[j * nxi - 1] = 0.0

        return sp.diags(
            diagonals=[main, east_west, east_west, north_south, north_south],
            offsets=[0, -1, +1, -nxi, +nxi],
            shape=(n, n),
            format="csr",
        )

    def current_density_from_strips(
        self,
        strips: Iterable[Strip],
        currents_A: Iterable[float],
    ) -> NDArray[np.float64]:
        """Return J_s(y,x) on the grid from rectangular strip currents."""
        J = np.zeros((self.ny, self.nx), dtype=float)
        X, Y = np.meshgrid(self.x, self.y, indexing="xy")

        for strip, current in zip(strips, currents_A):
            mask = (
                (np.abs(X - strip.x0) <= strip.width / 2.0)
                & (np.abs(Y - strip.y0) <= strip.thickness / 2.0)
            )
            n_cells = int(mask.sum())
            if n_cells == 0:
                raise ValueError(
                    f"Strip {strip.name!r} contains no grid cells. "
                    "Refine the mesh or enlarge strip dimensions."
                )
            # Normalize so the discrete current integral is exactly current.
            area_discrete = n_cells * self.dx * self.dy
            J[mask] += float(current) / area_discrete

        return J

    def solve_A(self, J_s: NDArray[np.float64]) -> NDArray[np.float64]:
        """Solve ∇² A_s = -mu0 J_s with A_s=0 on the boundary."""
        if J_s.shape != (self.ny, self.nx):
            raise ValueError(f"Expected J_s shape {(self.ny, self.nx)}, got {J_s.shape}.")
        rhs = -self.mu0 * J_s[1:-1, 1:-1].ravel()
        A_interior = self._solve(rhs)
        A = np.zeros_like(J_s)
        A[1:-1, 1:-1] = A_interior.reshape((self.ny - 2, self.nx - 2))
        return A

    def fields_from_A(
        self,
        A_s: NDArray[np.float64],
        *,
        convention: Literal["right_handed", "left_handed_flip"] = "right_handed",
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return Bx, By from A_s.

        right_handed:
            Bx =  dA/dy, By = -dA/dx

        left_handed_flip:
            Bx = -dA/dy, By =  dA/dx

        Choose one convention and use the same convention for the target data.
        """
        dA_dy, dA_dx = np.gradient(A_s, self.y, self.x, edge_order=2)
        if convention == "right_handed":
            return dA_dy, -dA_dx
        if convention == "left_handed_flip":
            return -dA_dy, dA_dx
        raise ValueError(f"Unknown convention {convention!r}.")

    def strip_unit_fields(
        self,
        strip: Strip,
        *,
        convention: Literal["right_handed", "left_handed_flip"] = "right_handed",
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """Return A_s, Bx, By for one ampere in one strip."""
        J = self.current_density_from_strips([strip], [1.0])
        A = self.solve_A(J)
        Bx, By = self.fields_from_A(A, convention=convention)
        return A, Bx, By

    def sample_field(
        self,
        field: NDArray[np.float64],
        points_xy: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Bilinearly interpolate a grid field at points [[x,y], ...]."""
        interp = RegularGridInterpolator(
            (self.y, self.x),
            field,
            bounds_error=False,
            fill_value=np.nan,
        )
        pts_yx = np.column_stack([points_xy[:, 1], points_xy[:, 0]])
        vals = interp(pts_yx)
        if np.isnan(vals).any():
            raise ValueError("Some sample points lie outside the Poisson domain.")
        return vals

    def response_matrix(
        self,
        strips: list[Strip],
        points_xy: NDArray[np.float64],
        *,
        component: Literal["Ix", "Iy"] = "Iy",
        convention: Literal["right_handed", "left_handed_flip"] = "right_handed",
    ) -> NDArray[np.float64]:
        """Build Q for first field integral at sample points.

        component="Ix" means length * Bx.
        component="Iy" means length * By.

        Q has shape (n_sample_points, n_strips), in T*m per ampere.
        """
        Q = np.zeros((len(points_xy), len(strips)), dtype=float)
        for j, strip in enumerate(strips):
            _, Bx, By = self.strip_unit_fields(strip, convention=convention)
            B = Bx if component == "Ix" else By
            Q[:, j] = strip.length * self.sample_field(B, points_xy)
        return Q


def solve_currents(
    Q: NDArray[np.float64],
    target: NDArray[np.float64],
    *,
    ridge: float = 0.0,
    current_limit_A: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Solve Q @ currents ~= target with optional ridge and current bounds.

    ridge is dimensionless after column normalization; values like 1e-6 to
    1e-3 are often useful to suppress large canceling currents.
    """
    target = np.asarray(target, dtype=float).ravel()
    if Q.shape[0] != target.size:
        raise ValueError("Q row count must match target length.")

    # Scale columns for numerical conditioning. Undo scaling after solve.
    col_norm = np.linalg.norm(Q, axis=0)
    col_norm[col_norm == 0.0] = 1.0
    Qs = Q / col_norm

    if ridge > 0:
        A_aug = np.vstack([Qs, np.sqrt(ridge) * np.eye(Qs.shape[1])])
        b_aug = np.concatenate([target, np.zeros(Qs.shape[1])])
    else:
        A_aug = Qs
        b_aug = target

    if current_limit_A is None:
        sol_scaled, *_ = np.linalg.lstsq(A_aug, b_aug, rcond=None)
    else:
        lim_scaled = current_limit_A * col_norm
        res = lsq_linear(A_aug, b_aug, bounds=(-lim_scaled, +lim_scaled), lsmr_tol="auto")
        if not res.success:
            print(f"Warning: bounded least-squares did not fully converge: {res.message}")
        sol_scaled = res.x

    currents = sol_scaled / col_norm
    fitted = Q @ currents
    return currents, fitted


def main() -> None:
    # --- Geometry close to the published NSLS-II current-strip examples.
    # 2012 paper example: 40 current strips, 6.5 m long, 2 mm wide,
    # 0.3 mm thick, 1 mm spacing. Here represented as 20 upper + 20 lower.
    strips = make_two_strip_rows(
        n_per_row=20,
        pitch=3.0e-3,       # width 2 mm + spacing 1 mm
        width=2.0e-3,
        thickness=0.3e-3,
        gap=10.7e-3,
        length=6.5,
    )

    solver = PoissonStripSolver(
        xlim=(-0.070, 0.070),
        ylim=(-0.035, 0.035),
        nx=241,
        ny=141,
    )

    # Sample in the horizontal mid-plane, as in the field-integral comparison
    # plots. Replace this target with your ID-derived correction integral.
    x_sample = np.linspace(-0.040, 0.040, 121)
    y_sample = np.zeros_like(x_sample)
    points = np.column_stack([x_sample, y_sample])

    # Synthetic target, in G*cm, only for demonstration.
    # Replace by the negative of the equivalent first field integral derived
    # from your nonlinear ID kick map, e.g. from RADIA / measured data.
    x_mm = 1.0e3 * x_sample
    target_gcm = (
        230.0 * np.exp(-((x_mm - 8.0) / 9.0) ** 2)
        - 190.0 * np.exp(-((x_mm + 12.0) / 11.0) ** 2)
        + 55.0 * np.sin(2.0 * np.pi * x_mm / 45.0)
    )
    target_Iy_Tm = gcm_to_t_m(target_gcm)

    # Build response matrix for vertical field integral Iy = ∫ By ds.
    # For a different local convention, set convention="left_handed_flip"
    # AND define the target using the same sign convention.
    Q = solver.response_matrix(
        strips,
        points,
        component="Iy",
        convention="right_handed",
    )

    currents_A, fitted_Iy_Tm = solve_currents(
        Q,
        target_Iy_Tm,
        ridge=1.0e-5,
        current_limit_A=5.0,
    )

    residual_Tm = fitted_Iy_Tm - target_Iy_Tm
    rms_gcm = float(np.sqrt(np.mean(t_m_to_gcm(residual_Tm) ** 2)))
    peak_gcm = float(np.max(np.abs(t_m_to_gcm(residual_Tm))))

    print(f"Fitted {len(currents_A)} strip currents.")
    print(f"RMS residual:  {rms_gcm:.3g} G*cm")
    print(f"Peak residual: {peak_gcm:.3g} G*cm")
    print(f"Peak current:  {np.max(np.abs(currents_A)):.3g} A")

    # Save current table.
    out = np.array(
        [
            [j, s.x0, s.y0, s.width, s.thickness, s.length, currents_A[j]]
            for j, s in enumerate(strips)
        ],
        dtype=object,
    )
    header = "index,x0_m,y0_m,width_m,thickness_m,length_m,current_A"
    np.savetxt(
        "strip_currents.csv",
        out,
        delimiter=",",
        header=header,
        comments="",
        fmt=["%d", "%.9e", "%.9e", "%.9e", "%.9e", "%.9e", "%.9e"],
    )

    # Plot target/fitted first integral and current distribution.
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), constrained_layout=True)

    axes[0].plot(1e3 * x_sample, t_m_to_gcm(target_Iy_Tm), label="target")
    axes[0].plot(1e3 * x_sample, t_m_to_gcm(fitted_Iy_Tm), "--", label="fitted")
    axes[0].plot(1e3 * x_sample, t_m_to_gcm(residual_Tm), ":", label="residual")
    axes[0].set_xlabel("x [mm]")
    axes[0].set_ylabel(r"$\int B_y\,ds$ [G cm]")
    axes[0].set_title("Current-strip first-integral fit")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    xs_mm = np.array([s.x0 for s in strips]) * 1e3
    ys_mm = np.array([s.y0 for s in strips]) * 1e3
    top = ys_mm > 0
    bottom = ~top
    axes[1].stem(xs_mm[top], currents_A[top], linefmt="C0-", markerfmt="C0o", basefmt="k-", label="top row")
    axes[1].stem(xs_mm[bottom], currents_A[bottom], linefmt="C1-", markerfmt="C1s", basefmt="k-", label="bottom row")
    axes[1].set_xlabel("strip x [mm]")
    axes[1].set_ylabel("current [A]")
    axes[1].set_title("Optimized strip currents")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    fig.savefig("current_strip_fit.png", dpi=180)
    print("Wrote strip_currents.csv and current_strip_fit.png")


if __name__ == "__main__":
    main()
