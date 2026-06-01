from __future__ import annotations

import json
from pathlib import Path
from typing import NamedTuple

class Config(NamedTuple):
    conc_ocean: float
    conc_small: float
    conc_large: float

    diff_ocean: float
    diff_small: float
    diff_large: float

    egest_small: float
    egest_large: float

    clear_small: float
    clear_large: float

    count_small: float
    count_large: float

    predation: float

    @staticmethod
    def load(path: str | Path):
        path = Path(path)
        with open(path, 'r') as f:
            data = json.load(f)
        return Config(**data)

    @staticmethod
    def save(config: Config, path: str | Path):
        path = Path(path)
        with open(path, 'w') as f:
            json.dump(config._asdict(), f, indent=4)