#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Export a double-precision analytic test functional, not a Skala model."""

import argparse
import json
from typing import Dict

import torch


class GridTestFunctional(torch.nn.Module):
    __constants__ = ["geometry"]

    def __init__(self, geometry: bool = False):
        super().__init__()
        self.geometry = geometry

    @torch.jit.export
    def get_exc_density(self, fields: Dict[str, torch.Tensor]) -> torch.Tensor:
        rho = fields["density"]
        tau = fields["kin"]
        grad = fields["grad"]
        weights = fields["atomic_grid_weights"]
        sizes = fields["atomic_grid_sizes"]
        # The model ABI uses (spin, point) and (spin, Cartesian axis, point).
        energy = (0.7 * rho**2 + 0.3 * tau**2 + 0.2 * rho * tau).sum(0)
        energy = energy + 0.1 * (grad**2).sum(0).sum(0)
        begin = 0
        # A nonlocal contribution within each atom block exercises descriptor adjoints.
        for block, size in enumerate(sizes):
            end = begin + int(size)
            local_weight = weights[begin:end]
            mean = (rho[:, begin:end].sum(0) * local_weight).sum() / local_weight.sum()
            energy[begin:end] = energy[begin:end] + 0.05 * mean**2
            if self.geometry:
                delta = (
                    fields["grid_coords"][begin:end]
                    - fields["coarse_0_atomic_coords"][block]
                )
                centre = (delta * local_weight[:, None]).sum(0) / local_weight.sum()
                energy[begin:end] = (
                    energy[begin:end] + 0.04 * (delta**2).sum(1) + 0.03 * (centre**2).sum()
                )
            begin = end
        return energy


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--geometry", action="store_true")
    args = parser.parse_args()
    features = ["density", "grad", "kin", "atomic_grid_weights", "atomic_grid_sizes"]
    if args.geometry:
        features += ["grid_coords", "coarse_0_atomic_coords"]
    model = torch.jit.script(GridTestFunctional(args.geometry))
    torch.jit.save(
        model,
        args.output,
        _extra_files={
            "protocol_version": "2",
            "features": json.dumps(features),
        },
    )
