import argparse
import os
from datetime import datetime
import fipy as fp
import numpy as np
from tqdm import tqdm

from config import Config
from utils import create_gaussian_blob, get_peak_concentration


def run_sim(cfg: Config, identifier: str, save: bool = True, no_tqdm: bool = False):
    run_id = datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + identifier

    nx, ny = 50, 50
    L = 1000
    A = float(L**2)
    dx, dy = L / nx, L / ny

    mesh = fp.PeriodicGrid2D(dx, dy, nx, ny)
    x, y = mesh.cellCenters

    m_ref = get_peak_concentration(L, cfg.conc_ocean, 0.10)
    ns_ref = get_peak_concentration(L, cfg.count_small, 0.2)
    nl_ref = get_peak_concentration(L, cfg.count_large, 0.05)

    def div_zero(a, b):
        if b == 0:
            return 0
        return a / b

    m_arr = create_gaussian_blob(
        mesh, div_zero(cfg.conc_ocean, m_ref), center_ratio=(0.5, 0.5), sigma_ratio=0.10
    )
    m = fp.CellVariable(name="Ocean Plastics (Normalized)", mesh=mesh, value=m_arr, hasOld=True)

    ns_arr = create_gaussian_blob(
        mesh, div_zero(cfg.count_small, ns_ref), center_ratio=(0.3, 0.7), sigma_ratio=0.20
    )
    ns = fp.CellVariable(name="Small Animals (Normalized)", mesh=mesh, value=ns_arr, hasOld=True)

    nl_arr = create_gaussian_blob(
        mesh, div_zero(cfg.count_large, nl_ref), center_ratio=(0.8, 0.2), sigma_ratio=0.05
    )
    nl = fp.CellVariable(name="Large Animals (Normalized)", mesh=mesh, value=nl_arr, hasOld=True)

    ps = fp.CellVariable(name="Plastics in Small (Normalized)", mesh=mesh, value=div_zero(1e-15, m_ref), hasOld=True)
    pl = fp.CellVariable(name="Plastics in Large (Normalized)", mesh=mesh, value=div_zero(1e-15, m_ref), hasOld=True)

    variables = (ns, nl, m, ps, pl)

    adj_cs_ns = cfg.clear_small * ns_ref
    adj_cl_nl = cfg.clear_large * nl_ref
    adj_beta_nl = cfg.predation * nl_ref

    eq_ns = fp.TransientTerm(var=ns) == fp.DiffusionTerm(coeff=cfg.diff_small, var=ns)
    eq_nl = fp.TransientTerm(var=nl) == fp.DiffusionTerm(coeff=cfg.diff_large, var=nl)

    eq_m = (
            fp.TransientTerm(var=m) == fp.DiffusionTerm(coeff=cfg.diff_ocean, var=m)
            + fp.ImplicitSourceTerm(coeff=-(adj_cs_ns * ns + adj_cl_nl * nl), var=m)
            + fp.ImplicitSourceTerm(coeff=cfg.egest_small, var=ps)
            + fp.ImplicitSourceTerm(coeff=cfg.egest_large, var=pl)
    )

    eq_ps = (
            fp.TransientTerm(var=ps) == fp.DiffusionTerm(coeff=cfg.diff_small, var=ps)
            + fp.ImplicitSourceTerm(coeff=(adj_cs_ns * ns), var=m)
            + fp.ImplicitSourceTerm(coeff=-cfg.egest_small, var=ps)
            + fp.ImplicitSourceTerm(coeff=-(adj_beta_nl * nl), var=ps)
    )

    eq_pl = (
            fp.TransientTerm(var=pl) == fp.DiffusionTerm(coeff=cfg.diff_large, var=pl)
            + fp.ImplicitSourceTerm(coeff=(adj_cl_nl * nl), var=m)
            + fp.ImplicitSourceTerm(coeff=(adj_beta_nl * nl), var=ps)
            + fp.ImplicitSourceTerm(coeff=-cfg.egest_large, var=pl)
    )

    coupled_eq = eq_ns & eq_nl & eq_m & eq_ps & eq_pl

    initial_dt = 5 # s
    dt = initial_dt
    max_dt = 15 # s
    total_time = 7 * 24 * 60 * 60 # s

    max_sweeps = 10
    tolerance = 1e-2

    t = 0
    step_count = 0

    history = {
        "t": [0],
        "m": [m.value.copy() * m_ref],
        "ps": [ps.value.copy() * m_ref],
        "pl": [pl.value.copy() * m_ref],
        "ns": [ns.value.copy() * ns_ref],
        "nl": [nl.value.copy() * nl_ref],
        "res": [0],
    }

    solver = fp.solvers.LinearLUSolver(tolerance=1e-5)

    for var in variables:
        var.updateOld()

    if not no_tqdm:
        pbar = tqdm(total=total_time, desc="Simulating", unit="s", mininterval=10,
                    bar_format="{l_bar}{bar}| {n:.1f}/{total:.1f}s [{elapsed}<{remaining}, {postfix}]")

    while t < total_time:
        if t + dt > total_time:
            dt = total_time - t

        for var in variables:
            var.value = var.old.copy()

        residual = 1.0
        sweeps = 0

        while residual > tolerance and sweeps < max_sweeps:
            residual = coupled_eq.sweep(dt=dt, solver=solver)
            sweeps += 1

        if residual <= tolerance:
            t += dt
            for var in variables:
                var.updateOld()

            if not no_tqdm:
                pbar.update(dt)
                pbar.set_postfix(sweeps=sweeps, res=f"{residual:.5f}", dt=f"{dt:.4f}s", status="OK")

            history["t"].append(t)
            history["m"].append(m.value.copy() * m_ref)
            history["ps"].append(ps.value.copy() * m_ref)
            history["pl"].append(pl.value.copy() * m_ref)
            history["ns"].append(ns.value.copy() * ns_ref)
            history["nl"].append(nl.value.copy() * nl_ref)
            history["res"].append(residual)

            step_count += 1

            if sweeps <= 4:
                dt = min(dt * 1.5, max_dt)
            elif sweeps >= 7:
                dt = dt * 0.8
        else:
            dt = dt * 0.5
            if not no_tqdm:
                pbar.set_postfix(sweeps=sweeps, res=f"{residual:.5f}", dt=f"{dt:.1f}s", status="RETRYING")
    if not no_tqdm:
        pbar.close()

    if save:
        os.makedirs("./data/runs", exist_ok=True)
        np.savez(f"./data/runs/run_{run_id}.npz", **history)
    return history

def parse_args():
    parser = argparse.ArgumentParser(description="Microplastic Simulation")

    parser.add_argument("-c", "--config", type=str, help="Config file for the simulation")
    parser.add_argument("-i", "--ident", type=str, help="Identifier for the run", default="")

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    run_sim(Config.load(args.config), args.ident)