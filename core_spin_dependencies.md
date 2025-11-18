# Spinach core ↔ utilities dependencies

This note summarises how the main **core workflow** functions depend on the **spin utility** functions (and vice versa). It is meant as a text equivalent of the canvas graph.

Core workflow set:
- `create.m`, `basis.m`, `hamiltonian.m`, `relaxation.m`, `operator.m`, `propagator.m`, `step.m`, `evolution.m`, `steady.m`, `thermalize.m`

Spin utilities set:
- `reduce.m`, `correlation.m`, `coherence.m`, `orientation.m`, `carrier.m`, `frqoffset.m`, `rotframe.m`, `spin.m`, `spinlock.m`, `decouple.m`, `homospoil.m`, `average.m`, `utilities/intrep.m`

---

## Core workflows using spin utilities

### `evolution.m` → `reduce.m`
- **Where**: `kernel/evolution.m:128,134,140,146,680`.
- **How**: before time stepping, `evolution` calls `reduce(spin_system,L,rho)` (and variants with `L'` and `destination`/`coil`) to compute **projectors** into smaller invariant subspaces.
- **Purpose**: shrink the Liouville space to the part that actually participates in the dynamics/observables, then propagate only in those reduced subspaces.

### `evolution.m` → `propagator.m`
- **Where**: `kernel/evolution.m` multiple call sites (`P=propagator(spin_system,...)`).
- **How**: for segments where an explicit propagator is convenient, `evolution` asks `propagator` to compute `P=exp(L * timestep)` (or with local L and dt), then applies `P` repeatedly to states.
- **Purpose**: implement time evolution via explicit propagator matrices rather than always using `step`-style exp-by-action.

### `thermalize.m` → `propagator.m`
- **Where**: `kernel/thermalize.m:84` (`R=R*propagator(spin_system,HLSPS,1i*beta);`).
- **How**: in the DiBari–Levitt branch, `thermalize` uses `propagator` on the Hamiltonian superoperator with an imaginary time step (`1i*beta`) to build an equilibrium-weighting operator.
- **Purpose**: incorporate Boltzmann weighting (via imaginary-time evolution) into the thermalised relaxation superoperator.

---

## Spin utilities using core workflows

### `average.m` → `propagator.m`
- **Where**: `kernel/average.m:160` (`P=propagator(spin_system,isergen(HL,HM,HR,time_step),time_step)*P;`).
- **How**: for matrix-log / exact averaging, `average` constructs an effective generator for one period and passes it to `propagator` to get the period propagator.
- **Purpose**: derive **effective average Hamiltonians** from stroboscopic propagators over a modulation cycle.

### `carrier.m` → `operator.m`
- **Where**: `kernel/carrier.m` loop over spins.
- **How**: for each selected spin, `carrier` calls `operator(spin_system,{''Lz''},{spin_index},operator_type)` to build the single-spin `Lz` operator/superoperator.
- **Purpose**: assemble the **carrier Hamiltonian** by summing Larmor-frequency-weighted `Lz` operators over spins.

### `frqoffset.m` → `operator.m`
- **Where**: `kernel/frqoffset.m` inside both offset-application branches.
- **How**: `frqoffset` calls `operator(spin_system,''Lz'',spin_label)` for each spin/channel, and adds `2*pi*offset * Lz` to the supplied Hamiltonian/superoperator.
- **Purpose**: apply **Larmor frequency offsets** in an approximate way by shifting the Hamiltonian with `Lz` terms.

### `spinlock.m` → `step.m` and `homospoil.m`
- **Where**: `kernel/spinlock.m:40–49`.
- **How**:
  - Uses `step(spin_system,Lx_or_Ly,rho,±pi/2)` to rotate the state into and out of the locking frame.
  - Uses `homospoil(spin_system,rho,''destroy'')` between the rotations to wipe all components except those aligned with the locking axis.
- **Purpose**: provide an **analytical approximation to a spin-lock** experiment by composing core propagation (`step`) with a filtering utility (`homospoil`).

### `rotframe.m` → `spin.m` and `utilities/intrep.m`
- **Where**: `kernel/rotframe.m`.
- **How**:
  - Calls `spin(isotope)` to obtain γ and multiplicity, then computes the rotation period `T` from γ and `spin_system.inter.magnet`.
  - Calls `intrep(spin_system,H0,H,T,order)` to perform the interaction-representation / rotating-frame transform.
- **Purpose**: build a **rotating-frame Hamiltonian** using physical spin data (from `spin`) and the generic interaction-representation engine (`intrep`).

### `utilities/intrep.m` → `propagator.m`
- **Where**: `kernel/utilities/intrep.m:39` (`P=propagator(spin_system,H0,T);`).
- **How**: `intrep` uses `propagator` to construct the evolution operator over one period under `H0`, then uses it to transform `H` into the interaction/rotating frame and/or to build effective Hamiltonians.
- **Purpose**: implement interaction-picture and period-averaged dynamics in terms of **matrix exponentials** provided by `propagator`.

---

## Utilities that mainly depend on data structures

The following utilities rely primarily on the **structures built by core workflows** (especially `create.m`, `basis.m`, `hamiltonian.m`, `relaxation.m`) rather than calling those functions directly:

- `reduce.m`: uses `spin_system.bas.*`, `spin_system.sys.*`, and the Liouvillian `L` to perform symmetry factorisation, zero-track elimination, and connected-component detection.
- `correlation.m`: uses `spin_system.bas.basis` and the supplied `orders`/`spins` to mask basis states by correlation order.
- `coherence.m`: uses `spin_system.bas.basis` and `lin2lm` to compute projection quantum numbers and filter by coherence order.
- `orientation.m`: consumes `Q` (built by `hamiltonian.m`) and Euler angles to reconstruct anisotropic `H`, but does not call core functions itself.
- `decouple.m`: uses `spin_system.bas.basis` and `L`/`rho` to zero all states involving selected spins, without calling other core functions.
- `homospoil.m`: uses `spin_system.bas.formalism`, `spin_system.bas.basis`, and `spin_system.inter.basefrqs` to decide which components of `rho` to keep/kill.
- `spin.m`: is a standalone database mapping isotope/particle names to γ and multiplicity; other functions call **it**, but it does not call cores.

These still conceptually depend on the **core setup** (they require a valid `spin_system` and basis), but the dependency is structural rather than through explicit function calls.

---

## Mental model

- **Setup layer**: `create` → `basis` → `hamiltonian` / `relaxation` / `operator` construct `spin_system`, basis, and generators.
- **Propagation layer**: `propagator`, `step`, and `evolution` implement different flavours of time evolution.
- **State-conditioning utilities**: `reduce`, `correlation`, `coherence`, `homospoil`, `decouple`, `spinlock` modify which parts of the state space are considered or populated.
- **Frame / averaging utilities**: `carrier`, `frqoffset`, `orientation`, `rotframe`, `intrep`, `average` reshape Hamiltonians and frames, often relying on `propagator` under the hood.
- **Steady/thermal layer**: `steady` and `thermalize` sit on top, using propagators and relaxation to fix steady or thermal states.

This file focuses only on **how** functions call one another; for detailed algorithm descriptions, see `core_workflows_docs.md` and `spin_utilities_docs.md`.
