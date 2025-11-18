# Spin Utilities Module Documentation (Spinach Kernel View)

In this project, the conceptual **`spin_utilities` module** refers to commonly used high-level utilities that manipulate states, operators, and Hamiltonians for analysis or control of spin dynamics. The functions documented here are:

- `reduce.m`
- `correlation.m`
- `coherence.m`
- `orientation.m`
- `carrier.m`
- `frqoffset.m`
- `rotframe.m`
- `spin.m`
- `spinlock.m`
- `decouple.m`
- `homospoil.m`
- `average.m`

For each function, we state its signature and summarise its core algorithm, ignoring boilerplate (validation, logging, etc.).

---

## `reduce.m`

**Signature**

```matlab
projectors = reduce(spin_system, L, rho)
```

- **Inputs**
  - `spin_system`: system with basis and symmetry information.
  - `L`: Liouvillian matrix.
  - `rho`: initial or destination state (used for screening).
- **Output**
  - `projectors`: cell array of projectors into independently evolving reduced subspaces.

**Core Logic and Purpose**

- If trajectory-level restriction is enabled, performs **multi-stage state-space reduction**:
  1. **Symmetry factorisation** (when permutation symmetry information is present and not disabled):
     - Uses symmetry labels (e.g. irreducible representation indices) to block-diagonalise the Liouvillian into symmetry sectors.
     - Builds projectors that isolate each symmetry block.
  2. **Zero-track elimination**:
     - Examines whether `rho` (or destination states) has support in each sector, and removes those sectors that are dynamically irrelevant.
  3. **Disconnected subspace identification**:
     - Treats the Liouvillian as a graph (non-zero couplings between states) and finds connected components via path tracing.
     - For each connected component, builds a projector onto its subspace.
- Returns projectors `P` that can be used as `L_reduced = P' * L * P`, `rho_reduced = P' * rho`, significantly shrinking the problem dimension without changing the dynamics of interest.

---

## `correlation.m`

**Signature**

```matlab
rho = correlation(spin_system, rho, orders, spins)
```

- **Inputs**
  - `spin_system`: system with `sphten-liouv` basis.
  - `rho`: state vector or horizontal stack of state vectors.
  - `orders`: row vector of correlation orders to keep.
  - `spins`: spin subset specification (`'all'`, isotope name, or numeric indices).
- **Output**
  - `rho`: same-shaped state vector(s) with only selected correlation orders retained.

**Core Logic and Purpose**

- Reshapes `rho` so that the first dimension corresponds to spin basis states and remaining dimensions (if any) represent additional indices.
- Parses `spins` into a list of spin indices, either by isotope name or explicit numbers.
- Computes the **order of spin correlation** for each basis state:
  - Uses `spin_system.bas.basis(:, spins)` to inspect which spins participate in each basis vector.
  - Sums the logical occupancy across the chosen spins to get the correlation order per basis state.
- Builds a boolean mask that is `true` only for basis states whose correlation order is in the `orders` list, and zeroes all other rows of `rho`.
- Restores the original shape of `rho`.
- Conceptually, this analytically enforces a correlation-order filter, often used as an alternative to long phase cycles.

---

## `coherence.m`

**Signature**

```matlab
rho = coherence(spin_system, rho, spec)
```

- **Inputs**
  - `spin_system`: system with `sphten-liouv` basis.
  - `rho`: state vector or horizontal stack thereof.
  - `spec`: cell array specifying which coherence orders to keep on which spins, e.g. `{{'13C',[1 -1]},{'1H',-1}}`.
- **Output**
  - `rho`: state vector(s) with only specified coherence orders retained.

**Core Logic and Purpose**

- Reshapes `rho` to separate spin basis dimension from extra dimensions.
- Computes projection quantum numbers for each basis state using `lin2lm(spin_system.bas.basis)`, obtaining `M` (magnetic quantum number contributions).
- For each element of `spec`:
  - Interprets the spin selection:
    - If the first element is a string, matches `'all'` or a specific isotope name to indices in `spin_system.comp.isotopes`.
    - If numeric, uses the indices directly.
  - Sums the relevant columns of `M` to obtain coherence order per basis state for the chosen spins.
  - Marks basis states whose coherence order lies in the specified set for that spec entry.
- Intersects the masks from all spec entries so that only states satisfying all coherence conditions are kept; all others are zeroed out in `rho`.
- Restores `rho` to its original dimensions.
- This implements an analytical **coherence-order filter**, matching typical NMR coherence selection logic.

---

## `orientation.m`

**Signature**

```matlab
H = orientation(Q, euler_angles)
```

- **Inputs**
  - `Q`: rotational basis (cell array of irreducible components) from `hamiltonian.m`.
  - `euler_angles`: `[alpha, beta, gamma]` (ZYZ) Euler angles in radians.
- **Output**
  - `H`: anisotropic part of the Hamiltonian at the specified orientation.

**Core Logic and Purpose**

- Initialises an empty sparse Hamiltonian `H`.
- Loops over spherical ranks `r`:
  - Computes the Wigner D‑matrix `D_r(alpha, beta, gamma)` for rank `r`.
  - For each pair of tensor indices `(k,m)`, if `Q{r}{k,m}` is non-zero, adds `D_r(k,m) * Q{r}{k,m}` to `H`.
- After summing contributions for all ranks, symmetrises the result as `(H + H') / 2` to enforce Hermiticity.
- This is a straightforward implementation of the **irreducible tensor rotation formula**, reconstructing the anisotropic Hamiltonian at arbitrary orientations using precomputed tensor components.

---

## `carrier.m`

**Signature**

```matlab
H = carrier(spin_system, spins, operator_type)
```

- **Inputs**
  - `spin_system`: system with Zeeman frequencies and basis.
  - `spins`: isotope name (e.g. `'1H'`), `'all'`, or similar.
  - `operator_type`: Liouville operator type (`'left'`, `'right'`, `'comm'`, `'acomm'`), ignored in Hilbert space.
- **Output**
  - `H`: carrier Hamiltonian (Hilbert operator or Liouville superoperator).

**Core Logic and Purpose**

- Identifies the indices of spins matching the `spins` specification.
- Initialises a sparse Hamiltonian/superoperator `H`.
- For each selected spin:
  - Takes its carrier (base) Larmor frequency from `spin_system.inter.basefrqs`.
  - Constructs the corresponding `Lz` operator (or superoperator) on that spin via `operator(spin_system, {'Lz'}, {spin_index}, operator_type)`.
  - Adds the frequency-weighted contribution to `H`.
- After summing contributions from all spins, enforces symmetry and cleans small elements via `clean_up`, returning `H` as the “carrier” part of the Zeeman Hamiltonian used in rotating frame transforms and average Hamiltonian theory.

---

## `frqoffset.m`

**Signature**

```matlab
H = frqoffset(spin_system, H, parameters)
```

- **Inputs**
  - `spin_system`: system with isotope list.
  - `H`: Hamiltonian operator or commutation superoperator.
  - `parameters`: structure with fields:
    - `spins`: cell array of spin labels (e.g. `{'1H','13C'}`).
    - `offset`: vector of frequency offsets (Hz) corresponding to `spins`.
- **Output**
  - `H`: Hamiltonian with Larmor frequency offsets applied.

**Core Logic and Purpose**

- Detects whether there are multiple channels referring to the same isotope:
  - If each entry in `parameters.spins` is unique:
    - For each spin/channel, adds `2*pi*offset * Lz(spin)` to `H`.
  - If some spins appear multiple times:
    - Uses `unique` and index mappings to ensure that all channels referring to the same spin share the same offset; errors if they do not.
    - Applies the unique offsets once per distinct spin.
- `Lz` operators are built via `operator(spin_system, 'Lz', spin_label)`.
- Conceptually, this function adds an approximate frequency shift to the Hamiltonian by applying Larmor offsets to specified spins, used in liquid-state NMR style treatments.

---

## `rotframe.m`

**Signature**

```matlab
Hr = rotframe(spin_system, H0, H, isotope, order)
```

- **Inputs**
  - `spin_system`: system with interaction and assumption information.
  - `H0`: carrier Hamiltonian (with respect to which the rotating frame is defined).
  - `H`: total lab-frame Hamiltonian `H0 + H1` to be transformed.
  - `isotope`: string specifying which spins define the rotating frame (e.g. `'1H'`).
  - `order`: perturbation theory order (can be `inf`).
- **Output**
  - `Hr`: rotating-frame Hamiltonian.

**Core Logic and Purpose**

- Determines the **rotation period** `T` associated with the selected isotope and formalism:
  - In Liouville space (`'zeeman-liouv'`, `'sphten-liouv'`), uses `T = -2*pi / (spin(isotope) * spin_system.inter.magnet)`.
  - In Hilbert space (`'zeeman-hilb'`), uses `T = -4*pi / (spin(isotope) * spin_system.inter.magnet)`.
- Calls the interaction representation helper:

  ```matlab
  Hr = intrep(spin_system, H0, H, T, order)
  ```

- `intrep` performs the formal rotating-frame/interaction-picture transformation, using an auxiliary-matrix method to efficiently compute the time-averaged effective Hamiltonian up to the specified order.
- Thus, `rotframe` computes the correct period `T` from physical parameters and invokes the generic interaction-picture engine with appropriate arguments to obtain the rotating-frame Hamiltonian.

---

## `spin.m`

**Signature**

```matlab
[gamma, multiplicity] = spin(name)
```

- **Inputs**
  - `name`: isotope or particle label (e.g. `'1H'`, `'13C'`, `'15N'`), or special symbols (`'G'`, `'E'`, `'E#'`, `'C#'`, `'V#'`, `'T#'`, `'N'`, `'M'`, etc.).
- **Outputs**
  - `gamma`: magnetogyric ratio (rad / (s·T)), or `0` for non-magnetic modes.
  - `multiplicity`: number of energy or population levels.

**Core Logic and Purpose**

- Implements a table-driven and pattern-driven **spin database**:
  - If `name` matches special patterns (`E#`, `C#`, `V#`, `T#`), interprets the suffix as multiplicity and sets `gamma` according to the particle type (e.g. electron vs cavity vs phonon vs transmon, with cavities/phonons/transmons having `gamma = 0`).
  - For standard nuclei and particles (e.g. `'1H'`, `'13C'`, `'195Pt'`, `'E'`, `'G'`, `'N'`, `'M'`), uses a big `switch` statement to assign hard-coded `gamma` and multiplicity values.
- Returns magnetogyric ratios from CODATA or literature where available and warns in comments where data may be approximate.
- All other parts of Spinach use this function to convert symbolic isotope labels into quantitative parameters.

---

## `spinlock.m`

**Signature**

```matlab
rho = spinlock(spin_system, Lx, Ly, rho, direction)
```

- **Inputs**
  - `spin_system`: system (for `step` and `homospoil` calls).
  - `Lx`: X magnetisation operator on spins to lock.
  - `Ly`: Y magnetisation operator on spins to lock.
  - `rho`: state vector or bookshelf stack thereof.
  - `direction`: `'X'` or `'Y'`, indicating locked axis.
- **Output**
  - `rho`: state after approximate spin-locking filter.

**Core Logic and Purpose**

- Depending on `direction`:
  - For `'X'`:
    1. Rotates the state into the Y‑basis using `step(spin_system, Ly, rho, pi/2)`.
    2. Applies `homospoil(spin_system, rho, 'destroy')` to obliterate all components except longitudinal ones in that frame.
    3. Rotates back with `step(spin_system, Ly, rho, -pi/2)`, leaving predominantly X‑magnetisation along the locked axis.
  - For `'Y'`:
    1. Rotates via `Lx`, homospoils, and rotates back, analogous to the X case.
- This implements an **analytical approximation** to a spin-lock experiment by combining idealised rotations and homospoil filtering instead of explicitly simulating strong RF fields.

---

## `decouple.m`

**Signature**

```matlab
[L, rho] = decouple(spin_system, L, rho, spins)
```

- **Inputs**
  - `spin_system`: system with basis and composition.
  - `L`: Liouvillian superoperator (may be empty).
  - `rho`: state vector or horizontal stack thereof (may be empty).
  - `spins`: spins to decouple; either a list of isotope names (cell array) or numeric indices.
- **Outputs**
  - `L`: Liouvillian with all terms involving target spins zeroed.
  - `rho`: state(s) with all populations involving target spins set to zero.

**Core Logic and Purpose**

- Builds a boolean mask `dec_mask` for spins to be decoupled, then:
  - Constructs a **state mask** `zero_mask` over basis states where any decoupled spin participates (using `spin_system.bas.basis`).
  - For Liouvillian:
    - Interprets the Liouville space as possibly a Fokker–Planck direct product of spin and space dimensions.
    - Kronecker‑extends `zero_mask` across space dimensions to obtain a full vector mask.
    - Zeros rows and columns of `L` corresponding to that mask, effectively removing any dynamics involving the target spins.
  - For `rho`:
    - Zeros entries corresponding to `zero_mask` (extended appropriately if there are additional dimensions).
- Mimics the effect of **perfect decoupling** of specified spins by analytically removing their interactions and populations from the state space without changing the basis or rebuilding `L` from scratch.

---

## `homospoil.m`

**Signature**

```matlab
rho = homospoil(spin_system, rho, zqc_flag)
```

- **Inputs**
  - `spin_system`: system with basis and base frequencies.
  - `rho`: state vector or density matrix (form depends on formalism).
  - `zqc_flag`: `'keep'` or `'destroy'` (ignored in Zeeman formalisms).
- **Output**
  - `rho`: filtered state, with only allowed components retained.

**Core Logic and Purpose**

- Behaviour depends on the **basis formalism**:
  - `zeeman-hilb`:
    - Keeps only the diagonal of the density matrix (longitudinal magnetisation), zeroing off-diagonal elements.
  - `zeeman-liouv`:
    - Reshapes the Liouville vector to a matrix, keeps only its diagonal, then vectorises back.
  - `sphten-liouv`:
    - Reshapes `rho` into spin × other dimensions, and uses projection quantum numbers and base frequencies:
      - Obtains `M` via `lin2lm(spin_system.bas.basis)`.
      - Forms a frequency for each basis state as a weighted sum of `spin_system.inter.basefrqs .* M`.
      - For `zqc_flag = 'keep'`: keeps only states whose carrier frequency is (approximately) zero.
      - For `zqc_flag = 'destroy'`: kills both transverse and zero-quantum components, keeping only pure longitudinal states.
- In all cases, the function implements an **idealised homospoil pulse** by analytically zeroing undesired components rather than simulating an explicit gradient or RF pulse.

---

## `average.m`

**Signature**

```matlab
H = average(spin_system, Hp, H0, Hm, omega, theory)
```

- **Inputs**
  - `spin_system`: system (used mainly for diagnostics).
  - `Hp`: part of the rotating-frame Hamiltonian with frequency `+omega`.
  - `H0`: zero-frequency part in the rotating frame.
  - `Hm`: part with frequency `-omega`.
  - `omega`: rotating-frame frequency (rad/s).
  - `theory`: averaging theory variant (`'ah_first_order'`, `'ah_second_order'`, `'ah_third_order'`, `'matrix_log'`, `'kb_first_order'`, `'kb_second_order'`, `'kb_third_order'`).
- **Output**
  - `H`: effective average Hamiltonian.

**Core Logic and Purpose**

- Depending on `theory`, constructs the average Hamiltonian with different combinations of commutators and powers of `1/omega`:
  - **Krylov–Bogolyubov (KB) orders**:
    - `kb_first_order`: `H ≈ H0 + (1/omega)*(Hp*Hm - Hm*Hp)`.
    - `kb_second_order`: adds `H2` built from nested commutators involving `Hp`, `Hm`, and `H0`.
    - Higher KB orders include further nested commutators and higher powers of `1/omega`.
  - **Average Hamiltonian (Waugh) orders**:
    - `ah_first_order`, `ah_second_order`, `ah_third_order` implement the corresponding terms from average Hamiltonian theory using combinations of commutators of `Hp`, `H0`, and `Hm`.
  - **`matrix_log`**:
    - Forms the exact Floquet/period propagator over one cycle and takes its matrix logarithm (using dense algebra) to obtain the effective Hamiltonian.
- In all cases, the function’s algorithm is to compute **effective, time-independent Hamiltonians** that approximate or exactly reproduce the stroboscopic dynamics of a periodically modulated system, using analytic commutator expansions or matrix-logarithm constructions.
