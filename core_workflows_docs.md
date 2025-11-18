# Core Workflows Module Documentation (Spinach Kernel View)

In this project, the conceptual **`core_workflows` module** refers to the main kernel functions that define the standard dataflow:

1. Build `spin_system`.
2. Construct basis and operators/Hamiltonians.
3. Assemble relaxation.
4. Propagate states or density matrices.
5. Compute steady/thermal states.

The functions documented here are:

- `create.m`
- `basis.m`
- `hamiltonian.m`
- `relaxation.m`
- `operator.m`
- `propagator.m`
- `step.m`
- `evolution.m`
- `steady.m`
- `thermalize.m`

For each function, boilerplate such as input validation, logging, and simple defaults is omitted; the focus is on the main algorithm and data flow. Return “types” are described semantically (MATLAB is dynamically typed).

---

## `create.m`

**Signature**

```matlab
spin_system = create(sys, inter)
```

- **Inputs**
  - `sys`: spin system and instrument specification structure.
  - `inter`: interaction specification structure.
- **Output**
  - `spin_system`: fully initialised Spinach kernel object.

**Core Logic and Purpose**

- Locates the Spinach installation root, runs global sanity checks, and applies any auto-executed overrides (`autoexec`).
- Normalises and records global options from `sys` (e.g., output mode, enabled/disabled features, scratch paths).
- Merges `sys` and `inter` into a single `spin_system` structure:
  - Stores system composition, interaction tensors, fields, assumptions, and numerical tolerances.
  - Sets up paths and scratch directories to be used by later kernel stages.
- Delegates detailed validation to helper routines but centralises the creation of the canonical `spin_system` object on which all other kernel functions operate.

---

## `basis.m`

**Signature**

```matlab
spin_system = basis(spin_system, bas)
```

- **Inputs**
  - `spin_system`: Spinach data structure from `create.m`.
  - `bas`: basis set specification structure.
- **Output**
  - `spin_system`: updated with basis information and derived structures.

**Core Logic and Purpose**

- Stores the basis specification (`bas`) inside `spin_system`, reports key basis settings, and then constructs the working basis according to the selected formalism (e.g. `sphten-liouv`, Zeeman variants).
- For `sphten-liouv`, performs connectivity analysis (e.g. IK‑DNP approximation):
  - Builds interaction graphs (electron–electron, electron–nuclear, nuclear–nuclear) from the coupling tensors, thresholding by `inter_cutoff`.
  - Uses these graphs to determine connected components and admissible correlations, which drive basis truncation and spin grouping.
- Constructs the actual basis arrays (e.g. `spin_system.bas.basis`) that encode which spins participate in which basis states and how tensor components are organised.
- The resulting basis data structure sets the dimensions and sparsity patterns for operators, Hamiltonians, and state vectors used throughout the rest of the kernel.

---

## `hamiltonian.m`

**Signature**

```matlab
[H, Q] = hamiltonian(spin_system, operator_type)
```

- **Inputs**
  - `spin_system`: fully initialised system with interactions and basis.
  - `operator_type` (Liouville only): `'left'`, `'right'`, `'comm'` (default), or `'acomm'`.
- **Outputs**
  - `H`: rotationally invariant part of the Hamiltonian (or corresponding Liouville superoperator).
  - `Q`: rotational basis (irreducible tensor components) for the anisotropic part.

**Core Logic and Purpose**

- Preallocates descriptor arrays that encode, for each spin and contribution slot:
  - Spin indices, operator labels, and isotropic coefficients.
  - Spherical tensor coefficients and irreducible components for anisotropic parts.
- Processes contributions from:
  - Zeeman interactions (including different strength models such as `'full'`, `'secular'`, etc.).
  - Quadrupolar interactions (NQI).
  - Spin–spin couplings (dipolar, J-coupling, etc.).
- For each interaction, decomposes the Hamiltonian into:
  - Isotropic part added directly into `H`.
  - Anisotropic parts encoded as spherical tensor components that populate `Q{rank}{k,m}`.
- Once descriptors are built, uses parallelised operator/superoperator generation calls (via `operator`) to assemble the actual matrices according to the requested `operator_type`.
- The key algorithmic idea is to separate rotationally invariant and anisotropic pieces so that arbitrary orientations can be reconstructed later by combining `Q` with `orientation.m`.

---

## `relaxation.m`

**Signature**

```matlab
R = relaxation(spin_system, euler_angles)
```

- **Inputs**
  - `spin_system`: system structure with relaxation settings.
  - `euler_angles` (optional): ZYZ Euler angles for orientation-dependent theories.
- **Output**
  - `R`: relaxation superoperator.

**Core Logic and Purpose**

- Starts from a zero (or preallocated) relaxation matrix and adds contributions from whatever relaxation theories are listed in `spin_system.rlx.theories`.
- For each supported theory (e.g. extended T1/T2, Redfield, cross-correlation, spectral density–based models):
  - Calls the appropriate specialised kernel, passing orientation information when required (`euler_angles` for anisotropic rates).
  - Receives one or more contribution matrices and sums them into the total `R`.
- Ensures that contributions are combined consistently (superoperators on the same state space) and that optional terms (e.g. cross-relaxation, relaxation anisotropy) are included only when requested.
- The central algorithm is an additive assembly of relaxation contributions from multiple theories into a single superoperator representing dissipative dynamics.

---

## `operator.m`

**Signature**

```matlab
A = operator(spin_system, operators, spins, operator_type, format)
```

- **Inputs**
  - `spin_system`: system and basis definition.
  - `operators`: operator label(s), either a string or cell array (e.g. `'Lz'`, `{'Lz','L+'}`).
  - `spins`: spin selection, either by isotope names (`'1H'`, `'13C'`, `'all'`) or by indices.
  - `operator_type` (Liouville): `'left'`, `'right'`, `'comm'`, or `'acomm'` (ignored in Hilbert space).
  - `format`: `'csc'` (sparse matrix) or `'xyz'` (triplet representation).
- **Output**
  - `A`: requested operator/superoperator in the specified format.

**Core Logic and Purpose**

- Parses calls into three main modes:
  1. Single operator on all spins of a given type (e.g. `'Lz'` on all `'13C'` spins).
  2. Single operator on an explicit list of spin indices.
  3. Product operators, where `operators` and `spins` are matched cell arrays describing multi-spin products.
- For each requested single-spin operator:
  - Maps labels (e.g. `'Lz'`, `'L+'`, `'Tl,m'`, central transition labels) into basis-level operator templates.
  - Uses the system basis and multiplicities to construct the sparse matrix for that operator.
  - If in Liouville space, wraps Hilbert operators into the desired superoperator type (left/right multiplication, commutator, or anticommutator).
- For multi-spin products, builds tensor products (or equivalent Kronecker structures) combined according to the basis representation rather than naive dense Kronecker products, preserving sparsity.
- Optionally uses an operator cache (`'op_cache'`) keyed by operator descriptors to reuse previously constructed operators.
- The algorithm balances combinatorial operator construction with sparse storage formats and, in Liouville space, with superoperator mapping logic.

---

## `propagator.m`

**Signature**

```matlab
P = propagator(spin_system, L, timestep)
```

- **Inputs**
  - `spin_system`: system with tolerances and GPU/cache options.
  - `L`: Hamiltonian or Liouvillian matrix.
  - `timestep`: propagation time step.
- **Output**
  - `P`: exponential propagator `exp(-1i * L * timestep)`.

**Core Logic and Purpose**

- Forms the generator `A = -1i * timestep * L`, then cleans it with `clean_up` to drop insignificant elements and adjust sparse/full storage.
- Estimates `||A||` with `cheap_norm`:
  - If the norm is very large, tightens the exponentiation tolerance or errors out.
  - Chooses the scaling factor `2^n_squarings` via `ceil(log2(mat_norm))`, then rescales `A` by this factor.
- Computes `exp(A / 2^n_squarings)` via a reordered Taylor series:
  - On GPU (when enabled and the matrix is large): runs a sparse-aware Taylor expansion in GPU memory, summing terms until the next term becomes structurally zero after clean-up.
  - On CPU: runs the same Taylor expansion, with multiplication order adapted for sparse vs full `A` to minimise cost and memory.
- After the series converges, cleans the resulting `P` again and then performs `n_squarings` squaring steps (`P = P * P`), optionally on GPU, to undo the initial scaling.
- Integrates with a hash-based disk cache (`'prop_cache'`) for repeated exponentiations with the same `L`, `timestep`, and tolerance.
- The main algorithm is a **scaled-and-squared Taylor series** for the matrix exponential with adaptive sparsity handling and optional GPU acceleration.

---

## `step.m`

**Signature**

```matlab
rho = step(spin_system, L, rho, time_step)
```

- **Inputs**
  - `spin_system`: system and basis definition.
  - `L`: Liouvillian or Hamiltonian; can be a single matrix, `{L_left, L_right}` (piecewise-linear), or `{L_left, L_mid, L_right}` (piecewise-quadratic).
  - `rho`: state vector or density matrix (or cell stack thereof).
  - `time_step`: length of time step.
- **Output**
  - `rho`: propagated state after one step.

**Core Logic and Purpose**

- Decides whether the evolution is:
  - `expm(A) * rho` (state vectors / wavefunctions in Liouville or wavefunction formalisms), or
  - `expm(A) * rho * expm(-A)` style (density matrices in Hilbert space), using formalism tests.
- Moves `L` and `rho` to GPU if requested and not already on GPU.
- For density-matrix evolution, converts the two- or three-point quadrature specification into an effective generator for the step, combining left/right/midpoint generators according to the quadrature rule.
- Uses a reordered Taylor expansion (similar in spirit to `propagator`, but applied as an action on `rho` instead of explicitly forming `P`), exploiting:
  - Sparsity-aware multiplication order.
  - Reuse of intermediate vectors/matrices to minimise memory footprint.
- For vector evolution, applies the same Taylor process directly to `rho`.
- The algorithm is effectively **`expm`-by-action** via reordered Taylor series and quadrature over time-dependent generators, avoiding explicit construction of full propagator matrices.

---

## `evolution.m`

**Signature**

```matlab
answer = evolution(spin_system, L, coil, rho, timestep, nsteps, output, destination)
```

- **Inputs (Liouville space)**
  - `L`: Liouvillian during evolution.
  - `coil`: detection state (or stack) for observable calculations.
  - `rho`: initial state vector or horizontal stack of states.
  - `timestep`: step duration.
  - `nsteps`: number of steps.
  - `output`: mode (`'final'`, `'trajectory'`, `'total'`, `'refocus'`, `'observable'`, `'multichannel'`).
  - `destination` (optional): state used for destination-based state-space screening.
- **Inputs (Hilbert space)**
  - `L`: Hamiltonian matrix.
  - `coil`, `rho`, `timestep`, `nsteps`, `output`: analogous roles for density matrices.
- **Output**
  - `answer`: depends on `output` (final state, trajectory, integrated observable, time series, etc.).

**Core Logic and Purpose**

- Detects whether the problem is in Liouville or Hilbert space, and dispatches appropriately.
- For Liouville space:
  - Optionally performs **trajectory-level state-space reduction** via `reduce` using `rho` (source screening) and/or `destination` (destination screening) to obtain projectors.
  - Projects `L`, `rho`, and `coil` into reduced subspaces and then runs time stepping:
    - Uses `propagator` and/or `step` to apply the evolution over each time step.
    - Depending on `output`, accumulates final state, full trajectory, refocused states, integrated observable, or multichannel observables.
- For Hilbert space:
  - Uses exponentiation-based propagation of density matrices and constructs trajectories or observables analogously, but in matrix form.
- The algorithm organises **high-level propagation workflows**, combining generator construction, optional reduction, time stepping, and observable evaluation into a single interface.

---

## `steady.m`

**Signature**

```matlab
rho = steady(spin_system, P, rho, tol, method)
```

- **Inputs**
  - `spin_system`: system definition with tolerances.
  - `P`: propagator (exponential of a thermalised Liouvillian), or product of such propagators for repeating blocks.
  - `rho` (optional): initial guess for the steady state.
  - `tol` (optional): convergence tolerance on the 2-norm difference between subsequent iterates.
  - `method` (optional): `'newton'` (default) or `'squaring'`.
- **Output**
  - `rho`: steady state under repeated application of `P`.

**Core Logic and Purpose**

- For the **`'squaring'` method**:
  - Iteratively squares the propagator `P` (`P ← P*P`) while monitoring the maximum absolute change between `P` and its square.
  - Once the difference falls below `tol`, interprets the resulting near-idempotent propagator as encoding the steady-state behaviour.
- For the **`'newton'` method**:
  - Uses a Newton–Raphson style fixed-point iteration for `rho = P * rho`, with the constraint that the trace component remains unity.
  - Iteratively updates `rho` until the change between iterations falls below `tol`.
- Both methods exploit the fact that the steady state is a fixed point of the propagator and differ mainly in whether they iterate on the propagator (`'squaring'`) or directly on the state (`'newton'`), trading off robustness and computational cost.

---

## `thermalize.m`

**Signature**

```matlab
R = thermalize(spin_system, R, HLSPS, T, rho_eq, method)
```

- **Inputs**
  - `spin_system`: system with basis and interaction assumptions.
  - `R`: relaxation superoperator that currently drives towards the zero state.
  - `HLSPS`: lab-frame Hamiltonian left-side product superoperator (for DiBari–Levitt).
  - `T`: absolute temperature (for DiBari–Levitt).
  - `rho_eq`: target equilibrium state (for IME).
  - `method`: `'IME'` or `'dibari'`.
- **Output**
  - `R`: modified (thermalised) relaxation superoperator.

**Core Logic and Purpose**

- For `method = 'IME'` (inhomogeneous master equation):
  - Constructs the **unit state** `U` appropriate to the formalism:
    - In `sphten-liouv`, a single basis state with unit population of the unit operator component.
    - In `zeeman-liouv`, a vectorised identity matrix.
  - Applies an inhomogeneous correction to `R` so that `rho_eq` becomes a stationary solution:

    ```matlab
    R = R - kron(U', R * rho_eq)
    ```

- For `method = 'dibari'` (DiBari–Levitt):
  - Uses the lab-frame Hamiltonian superoperator `HLSPS` and temperature `T` to construct an equilibrium-driven inhomogeneous term based on detailed balance.
  - Adjusts `R` so that the long-time limit of the dynamics corresponds to thermal equilibrium with respect to `HLSPS` at temperature `T`.
- The central idea is to transform a homogeneous relaxation superoperator that drives towards zero into one whose fixed point is either a user-specified `rho_eq` (IME) or a physically motivated thermal equilibrium state (DiBari–Levitt).
