#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Export a double-precision analytic test functional, not a Skala model."""

import json
import sys
from typing import Dict

import torch


class GridTestFunctional(torch.nn.Module):
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
        for size in sizes:
            end = begin + int(size)
            local_weight = weights[begin:end]
            mean = (rho[:, begin:end].sum(0) * local_weight).sum() / local_weight.sum()
            energy[begin:end] = energy[begin:end] + 0.05 * mean**2
            begin = end
        return energy


if __name__ == "__main__":
    model = torch.jit.script(GridTestFunctional())
    torch.jit.save(
        model,
        sys.argv[1],
        _extra_files={
            "protocol_version": "2",
            "features": json.dumps(
                ["density", "grad", "kin", "atomic_grid_weights", "atomic_grid_sizes"]
            ),
        },
    )
