import numpy as np
import fipy as fp

from fipy.meshes.abstractMesh import AbstractMesh

def load_data(path: str) -> tuple:
    data = np.load(path)
    t = data["t"]
    m = data["m"]

    n = int(len(m[0]) ** 0.5)

    d = dict()
    d["t"] = data["t"]
    d["res"] = data["res"]
    d["m"] = data["m"].reshape(-1, n, n)
    d["ps"] = data["ps"].reshape(-1, n, n)
    d["pl"] = data["pl"].reshape(-1, n, n)
    d["ns"] = data["ns"].reshape(-1, n, n)
    d["nl"] = data["nl"].reshape(-1, n, n)

    return d, n

def get_peak_concentration(L: float, mean_conc: float, sigma_ratio: float):
    sigma_x = L * sigma_ratio
    sigma_y = L * sigma_ratio
    return (mean_conc * L**2) / (2 * np.pi * sigma_x * sigma_y)

def create_gaussian_blob(
        mesh: AbstractMesh,
        mean_concentration: float,
        center_ratio: tuple[float, float] = (0.5, 0.5),
        sigma_ratio:float=0.1
) -> np.ndarray:
    x, y = mesh.cellCenters
    x_max, x_min = float(x.max()), float(x.min())
    y_max, y_min = float(y.max()), float(y.min())

    L_x = x_max - x_min
    L_y = y_max - y_min
    total_area = L_x * L_y

    V = mean_concentration * total_area

    x0 = x_min + (L_x * center_ratio[0])
    y0 = y_min + (L_y * center_ratio[1])

    sigma_x = L_x * sigma_ratio
    sigma_y = L_y * sigma_ratio

    A = V / (2 * np.pi * sigma_x * sigma_y)

    gaussian_array = A * fp.numerix.exp(
        -(((x - x0) ** 2 / (2 * sigma_x ** 2)) + ((y - y0) ** 2 / (2 * sigma_y ** 2)))
    )

    return gaussian_array