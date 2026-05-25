# Corrected UAV-MEC Simulation Codebase

This codebase implements the corrected two-stage UAV-assisted MEC paper:
joint trajectory optimisation, adaptive relay scheduling, and resource
allocation. It replaces the older `Claude__2_` code, which did NOT match
the paper (it had no `phi`, no `b[n]` as a variable, and used a quadratic
centroid-pull heuristic instead of SCA).

## Requirements
- MATLAB (tested target R2025b)
- CVX 2.2 on the MATLAB path  (you already have this)

## How to run
1. Put all `.m` files in one folder.
2. In MATLAB, `cd` into that folder.
3. Run:  `main_simulation`
4. Figures are saved as PNG in the same folder; results in `results.mat`.

## Files
- `parameters.m`            - all configuration; regime documented inline
- `trajectory_init.m`       - initial straight-line trajectory
- `channel_model.m`         - MD->UAV channel gains and rates
- `relay_rate.m`            - UAV->TBS relay rate per slot
- `energy_model.m`          - five-component energy (E1..E5)
- `SP1_resource_allocation.m` - Algorithm 1: inner BCD over {f},{tau,b},{l},{phi}
- `SP1_uniform_relay.m`     - SP1 variant with UNIFORM b[n] (baseline)
- `SP2_trajectory_SCA.m`    - Algorithm 2: SCA with linearised log-rates
- `bcd_sca_solve.m`         - Algorithm 3: outer BCD-SCA loop
- `baseline_scheme.m`       - local / full / static / equal / uniform baselines
- `main_simulation.m`       - top-level script; produces all figures
- `plot_results.m`          - figure generation

## What changed vs the old code (and why)
1. Added `phi` (relay ratio). The old code had no `phi`, so UAV computation
   energy E3 was structurally zero. Now the two-stage split is real.
2. `b[n]` is a true optimisation variable, solved inside SP1 (Algorithm 1,
   Block 2 = LP). The old `relay_scheduling.m` heuristic is removed.
3. SP1 is a 4-block inner BCD exactly as Algorithm 1 describes.
4. SP2 now uses first-order Taylor linearisation of the log-rates
   (paper Section V-D-3). The old quadratic centroid-pull is removed.
5. Delay model is slotted: T = N*dt is the deadline. The contradictory
   per-user T_u term and the wrong relay-delay equation are gone.
6. Power and bandwidth are fixed constants (P_u, B), matching the paper.
7. The convergence curve is honest: no line-search / revert trick.
   If SP2 is correct, the energy decreases on its own.

## Parameter regime note (for your defence)
Horizon T is short (40 s) and task load high (70 Mbit) so the UAV is
time-pressured and must position itself well to meet deadlines. The
propulsion constants are for a light rotary-wing UAV so that propulsion
energy E5 is the SAME ORDER as communication+relay energy E2+E4
(ratio ~0.8). This makes the trajectory trade-off visible and defensible.
If an examiner asks "why these numbers", the answer is: they are chosen
so that no single energy term trivially dominates, which is the regime
in which joint trajectory + resource optimisation is meaningful.

## Expected behaviour
- Convergence curve: total energy decreases over the outer BCD iterations
  and flattens (it should NOT be a flat line).
- Trajectories: UAVs bend away from the straight line toward user
  clusters / the TBS, because doing so raises rates and helps meet the
  tight deadline.
- Adaptive vs uniform relay: adaptive should show lower relay energy E4,
  because it concentrates relay bits in high-rate slots.

## Honest caveat
This code is written carefully but has NOT been executed here (no MATLAB
in the authoring environment). Expect to debug 1-3 rounds: run it, copy
any error message, and fix. The most likely first errors are CVX DCP
complaints in SP2 or SP1 Block-2/3/4 -- those are fixable by adjusting
how an expression is written, not by changing the method.

If after a correct run the adaptive-relay gain is small (say 5-10%),
that is still a valid result. Report it honestly. Do not re-introduce
any trick to inflate it.
