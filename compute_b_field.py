import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.interpolate import RegularGridInterpolator, interp1d
from idstools.cocos import IDS_COCOS, COCOS


def compute_b_field_components(eq, itime=0, i1=0, overwrite=False):
    """
    Compute B field components from psi2d on a rectilinear (R,Z) mesh.

    Parameters
    ----------
    eq : IMAS 'equilibrium' IDS handle
        Equilibrium IDS containing the magnetic field data.
    itime : int, default 0
        Time-slice index.
    i1 : int, default 0
        Index into profiles_2d[i1] for 2D fields.
    overwrite : bool, default False
        Controls magnetic field component computation behavior:
        - If True: Always compute and store magnetic field components, even if they already exist
        - If False: Only compute and store magnetic field components if they don't already exist

    Returns
    -------
    eq : IMAS 'equilibrium' IDS handle
        Modified equilibrium IDS with computed magnetic field components.
    """
    # Get data from IDS
    ts = eq.time_slice[itime]
    vf = eq.vacuum_toroidal_field
    p2 = ts.profiles_2d[i1]

    r2d = p2.r
    z2d = p2.z
    psi2d = p2.psi

    # Check if magnetic field components are already available
    has_br = p2.b_field_r.has_value
    has_bz = p2.b_field_z.has_value
    has_bphi = p2.b_field_phi.has_value

    # Compute magnetic field components if not available and overwrite is True
    if (not has_br or not has_bz or not has_bphi) or overwrite:

        r_axis_1d = r2d[:, 0]
        z_axis_1d = z2d[0, :]

        dpsi_dr = np.gradient(psi2d, r_axis_1d, axis=0, edge_order=2)
        dpsi_dz = np.gradient(psi2d, z_axis_1d, axis=1, edge_order=2)

        R0 = vf.r0
        B0 = vf.b0[0]

        cv = COCOS(
            index={
                "COCOS": IDS_COCOS,
                "ipsign": np.sign(ts.global_quantities.ip),
                "b0sign": np.sign(B0),
            }
        )
        bsgn = cv.sigma_rphi_z * cv.sigma_bp / np.power(2.0 * np.pi, cv.exp_bp)

        Br2d = +1.0 * bsgn / r2d * dpsi_dz
        Bz2d = -1.0 * bsgn / r2d * dpsi_dr
        Bphi2d = B0 * R0 / r2d

        # Store computed components in the equilibrium IDS
        if (not has_br) or overwrite:
            p2.b_field_r = Br2d
        if (not has_bz) or overwrite:
            p2.b_field_z = Bz2d
        if (not has_bphi) or overwrite:
            p2.b_field_phi = Bphi2d

    return eq


def compute_b_field_fsa(
    eq,
    itime: int = 0,
    i1: int = 0,
    npoints: int = 512,
    eps_bp: float = 1e-12,
    return_loops: bool = False,
    overwrite: bool = False,
):
    """
    Compute flux-surface average (FSA) of magnetic field magnitude |B| on psi=const contours,
    using the axisymmetric Jacobian weighting 1/Bp:
        <|B|>(psi) = (∮ |B|/Bp dl) / (∮ 1/Bp dl),  with Bp = sqrt(Br^2 + Bz^2).

    The function uses equilibrium.time_slice[itime].profiles_1d.psi as the coordinate
    and populates b_field_average, b_field_min, and b_field_max in the same coordinate.

    ASSUMPTION: (R,Z) mesh is rectilinear and strictly monotonic. We use only
    RegularGridInterpolator for field sampling.

    Magnetic axis (R_axis, Z_axis) is taken from:
        equilibrium.time_slice[itime].global_quantities.magnetic_axis.(r, z)

    Parameters
    ----------
    eq : IMAS 'equilibrium' IDS handle
    itime : int, default 0
        Time-slice index.
    i1 : int, default 0
        Index into profiles_2d[i1] for 2D fields.


    npoints : int, default 512
        Number of uniformly arc-length sampled points per closed loop.
    eps_bp : float, default 1e-12
        Lower bound for Bp to avoid division-by-zero near X-points.
    return_loops : bool, default False
        If True, include per-loop sampled arrays {"R","Z","Bp","B"} in the output.
    overwrite : bool, default False
        Controls magnetic field profiles_1d population behavior:
        - If True: Always compute and store b_field_average, b_min, and b_max, even if they already exist
        - If False: Only compute and store magnetic field profiles if they don't already exist

    Returns
    -------
    eq : IMAS 'equilibrium' IDS handle
        Modified equilibrium IDS with populated profiles_1d magnetic field data
        (if overwrite=True).
    out : dict
        {
          "psi": (M,),                      # psi values that produced valid closed loops
          "psi_norms": (M,),                # psi_norm = sqrt((psi-psi_axis_ref)/(psi_sep_ref-psi_axis_ref))
          "fsa": (M,),                      # <|B|>(psi) with 1/Bp weighting
          "fmin": (M,), "fmax": (M,),       # min/max |B| along each loop
          "R_axis": float, "Z_axis": float,
          "psi_axis_ref": float, "psi_sep_ref": float,
          "loops": [ { "R":(Ni,), "Z":(Ni,), "Bp":(Ni,), "B":(Ni,) }, ... ],
        }

    Notes
    -----
    * For each level, we select the unique closed loop whose centroid is closest
      to the magnetic axis and resample it uniformly by arc length so that Δl is constant.
    * FSA then reduces to a weighted mean with weights 1/Bp on the sampled loop.
    * This function computes flux surface averages of magnetic field magnitude |B| only.
    """

    # --- 0) Pull data from IDS ---
    # Ensure magnetic field components are available
    eq = compute_b_field_components(eq, itime=itime, i1=i1, overwrite=False)

    ts = eq.time_slice[itime]
    vf = eq.vacuum_toroidal_field
    p1 = ts.profiles_1d
    p2 = ts.profiles_2d[i1]

    # Check grid type if available
    if grid_index != 1:
        warnings.warn(
            f"Non-rectilinear grid detected (index={grid_index})."
            UserWarning
        )

    # Check if magnetic field profiles are already available
    has_b_min = p1.b_field_min.has_value
    has_b_max = p1.b_field_max.has_value
    has_b_average = p1.b_field_average.has_value
    if (has_b_min and has_b_max and has_b_average) and (not overwrite):
        return eq, None

    R2d = p2.r
    Z2d = p2.z
    psi2d = p2.psi
    Br2d = p2.b_field_r
    Bz2d = p2.b_field_z
    Bphi2d = p2.b_field_phi

    # psgn > 0 as Psi increasing, psgn < 0 as decreasing
    cv = COCOS(
        index={
            "COCOS": IDS_COCOS,
            "ipsign": np.sign(ts.global_quantities.ip),
            "b0sign": np.sign(vf.b0[0]),
        }
    )
    psgn = cv.sigma_bp * cv.sigma_ip

    # --- 1) Rectilinear, strictly monotonic mesh (required) ---
    r_axis_1d = R2d[:, 0]
    z_axis_1d = Z2d[0, :]

    # --- 2) Magnetic axis and psi refs ---
    R_axis = float(ts.global_quantities.magnetic_axis.r)
    Z_axis = float(ts.global_quantities.magnetic_axis.z)

    # Always use profiles_1d.psi as the coordinate
    if not ts.profiles_1d.psi.has_value:
        raise ValueError(
            "profiles_1d.psi is required but not available in the equilibrium IDS"
        )

    psi1d = ts.profiles_1d.psi
    psi_axis_ref = float(psi1d[0])  # assumed near axis
    psi_sep_ref = float(psi1d[-1])  # assumed near separatrix

    if psi1d.size == 0:
        raise ValueError("No valid psi levels found strictly inside the 2D psi range.")

    # Normalized psi for reporting
    psi_norm_all = (psi1d - psi_axis_ref) / (psi_sep_ref - psi_axis_ref)

    # --- 3) Interpolators for fields on (R,Z) ---
    Br_i = RegularGridInterpolator(
        (r_axis_1d, z_axis_1d), Br2d * psgn, bounds_error=False, fill_value=np.nan
    )
    Bz_i = RegularGridInterpolator(
        (r_axis_1d, z_axis_1d), Bz2d * psgn, bounds_error=False, fill_value=np.nan
    )
    Bphi_i = RegularGridInterpolator(
        (r_axis_1d, z_axis_1d), Bphi2d * psgn, bounds_error=False, fill_value=np.nan
    )

    # --- 4) Build contours (we harvest paths only) ---
    fig, ax = plt.subplots()
    try:
        # flipping sign due to the requirement of the module contour
        cs = ax.contour(R2d, Z2d, psi2d * psgn, levels=psi1d * psgn)
    finally:
        plt.close(fig)

    # --- 5) Process each level: choose loop closest to axis, resample, evaluate f, compute FSA ---
    valid_levels, psi_norm_out = [], []
    fsa_list, fmin_list, fmax_list = [], [], []
    loops = []

    for lev, collection, psi_norm in zip(cs.levels, cs.collections, psi_norm_all):
        paths = collection.get_paths()
        if not paths:
            continue

        # Pick unique closed loop whose centroid is closest to magnetic axis
        chosen_v, dmin = None, np.inf
        for path in paths:
            v = path.vertices  # (N,2): columns are (R,Z)
            if v.shape[0] < 3:
                continue
            Rc, Zc = v[:, 0].mean(), v[:, 1].mean()
            d = np.hypot(Rc - R_axis, Zc - Z_axis)
            if d < dmin:
                dmin, chosen_v = d, v
        if chosen_v is None:
            continue

        # Uniform arc-length resampling of the closed loop
        x, y = chosen_v[:, 0], chosen_v[:, 1]
        x_closed = np.r_[x, x[0]]
        y_closed = np.r_[y, y[0]]
        ds = np.sqrt(np.diff(x_closed) ** 2 + np.diff(y_closed) ** 2)
        s = np.concatenate([[0.0], np.cumsum(ds)])
        L = s[-1]
        if L <= 0:
            continue

        s_uniform = np.linspace(0, L, npoints, endpoint=False)
        fx = interp1d(s, x_closed, kind="linear")
        fy = interp1d(s, y_closed, kind="linear")
        Ru = fx(s_uniform)
        Zu = fy(s_uniform)

        # Interpolate fields along the loop
        ptsu = np.column_stack([Ru, Zu])
        Br_u = Br_i(ptsu)
        Bz_u = Bz_i(ptsu)
        Bphi_u = Bphi_i(ptsu)

        mask = np.isfinite(Br_u) & np.isfinite(Bz_u) & np.isfinite(Bphi_u)
        if mask.sum() < max(8, npoints // 4):
            continue

        Ru, Zu = Ru[mask], Zu[mask]
        Br_u, Bz_u, Bphi_u = Br_u[mask], Bz_u[mask], Bphi_u[mask]

        Bp_u = np.sqrt(Br_u**2 + Bz_u**2)
        Bp_u = np.maximum(Bp_u, eps_bp)  # avoid division by zero
        B_u = np.sqrt(Br_u**2 + Bz_u**2 + Bphi_u**2)

        # FSA: <|B|> = (Σ |B|/Bp) / (Σ 1/Bp) on uniformly arc-length sampled loop
        w = 1.0 / Bp_u
        num = np.sum(B_u * w)
        den = np.sum(w)
        if den <= 0 or not np.isfinite(num / den):
            continue

        b_fsa = float(num / den)

        valid_levels.append(float(lev) * psgn)
        psi_norm_out.append(float(psi_norm))
        fsa_list.append(b_fsa)
        fmin_list.append(float(np.min(B_u)))
        fmax_list.append(float(np.max(B_u)))

        if return_loops:
            loops.append({"R": Ru, "Z": Zu, "Bp": Bp_u, "B": B_u})

    if len(valid_levels) == 0:
        raise RuntimeError(
            "No valid closed contours found for the requested psi levels."
        )

    # Ensure arrays match the profiles_1d.psi coordinate
    # Create arrays with the same size as psi1d, filling missing values with NaN
    n_psi = len(psi1d)
    b_avg_full = np.full(n_psi, np.nan)
    b_min_full = np.full(n_psi, np.nan)
    b_max_full = np.full(n_psi, np.nan)

    # Map computed values to the correct psi indices
    for i, psi_val in enumerate(valid_levels):
        # Find the closest index in psi1d
        idx = np.argmin(np.abs(psi1d - psi_val))
        b_avg_full[idx] = fsa_list[i]
        b_min_full[idx] = fmin_list[i]
        b_max_full[idx] = fmax_list[i]

    # Helper function to interpolate/extrapolate NaN values
    def fill_nan_values(data_array, coordinate_array):
        """Fill NaN values using best available interpolation method."""
        valid_idx = ~np.isnan(data_array)
        valid_count = np.sum(valid_idx)

        # Determine best interpolation method based on available data points
        if valid_count >= 4 and np.any(~valid_idx):
            kind = "cubic"
        elif valid_count >= 3 and np.any(~valid_idx):
            kind = "quadratic"
        elif valid_count >= 2 and np.any(~valid_idx):
            kind = "linear"
        else:
            return data_array  # Not enough points for interpolation

        interp_func = interp1d(
            coordinate_array[valid_idx],
            data_array[valid_idx],
            kind=kind,
            fill_value="extrapolate",
        )
        data_array[~valid_idx] = interp_func(coordinate_array[~valid_idx])

        return data_array

    # Interpolate/extrapolate to fill NaN values for all arrays
    b_avg_full = fill_nan_values(b_avg_full, psi_norm_all)
    b_min_full = fill_nan_values(b_min_full, psi_norm_all)
    b_max_full = fill_nan_values(b_max_full, psi_norm_all)

    # --- 6) Prepare output ---
    # Populate profiles_1d if not available or overwrite is True
    if (not has_b_average) or overwrite:
        p1.b_field_average = b_avg_full
    if (not has_b_min) or overwrite:
        p1.b_field_min = b_min_full
    if (not has_b_max) or overwrite:
        p1.b_field_max = b_max_full

    out = {
        "psi": np.asarray(valid_levels),
        "psi_norm": np.asarray(psi_norm_out),
        "fsa": np.asarray(fsa_list),
        "fmin": np.asarray(fmin_list),
        "fmax": np.asarray(fmax_list),
        "R_axis": R_axis,
        "Z_axis": Z_axis,
        "psi_axis_ref": psi_axis_ref,
        "psi_sep_ref": psi_sep_ref,
    }

    # Add loops information if return_loops is True
    if return_loops:
        out["loops"] = loops

    return eq, out
