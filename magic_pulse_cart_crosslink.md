## Mapping `magic_pulse_cart.m` to Skinner–Kobzar–Luy–Bendall–Bermel–Khaneja–Glaser (2006) broadband pulse paper

This document explains how the Spinach example `examples/optimal_control/magic_pulse_cart.m` implements, in code, the optimal-control broadband excitation pulse described in

- T. E. Skinner, K. Kobzar, B. Luy, M. R. Bendall, W. Bermel, N. Khaneja, S. J. Glaser,  
  “Optimal control design of constant amplitude phase-modulated pulses: Application to calibration-free broadband excitation”,  
  J. Magn. Reson. 179 (2006) 241–249.

The key idea is that `magic_pulse_cart.m` sets up **the same control problem**—calibration-free broadband excitation with RF inhomogeneity and resonance offsets—but expressed in Spinach’s Liouville-space / GRAPE framework.

---

### 1. Physical target: broadband, RF‑robust excitation

**Code (top comments and system setup in `magic_pulse_cart.m`):**

- Comments describe a 13C 90° excitation at 28.18 T, over ~200 ppm (~60 kHz) bandwidth, with pulse duration short enough that worst‑case 13C–1H J-coupling (~200 Hz) is negligible.
- The magnet and chemical shift range:
  - `sys.magnet=28.18;`
  - `n_spins=100; sys.isotopes=cell(n_spins,1);`
  - `inter.zeeman.scalar=num2cell(linspace(-100,100,n_spins));`

**Paper crosslink:**

- The paper’s design goal: broadband excitation over a large chemical shift range with high tolerance to RF inhomogeneity/miscalibration:
  - “Design criteria were transformation of Iz → Ix over resonance offsets of ±25 kHz for constant RF amplitude anywhere in the range 10–20 kHz…”
- The code mirrors this:
  - Many offset samples (`n_spins` over −100…+100 ppm) form a discrete ensemble of resonance offsets.
  - A strong field (28.18 T) and wide bandwidth (~60 kHz) express the same “broadband excitation” requirement, adapted to this example.

---

### 2. Liouville‑space formulation and basis

**Code:**

- Basis choice and system construction:
  - `bas.formalism='sphten-liouv';`
  - `bas.approximation='IK-2';`
  - `spin_system=create(sys,inter);`
  - `spin_system=basis(spin_system,bas);`

**Paper crosslink:**

- The paper uses the Bloch / Liouville–von Neumann equation for magnetization under a time‑dependent Hamiltonian with RF controls and offsets.
- Spinach implements this in Liouville space (`sphten-liouv`), with a reduced basis (`IK-2`) that keeps single‑spin operators like Ix, Iy, Iz and discards higher‑order correlations, matching the effectively single‑spin control picture of the paper.

---

### 3. Initial and target states: encoding the excitation/UR requirements

**Code:**

- Build and normalize Cartesian single‑spin operators:
  - `Sx=state(spin_system,'Lx','13C');`
  - `Sy=state(spin_system,'Ly','13C');`
  - `Sz=state(spin_system,'Lz','13C');`
  - Each is normalized by its 2‑norm.
- Define multiple initial and target states for the optimization:
  - `control.rho_init={ Sx Sy Sz};`
  - `control.rho_targ={-Sz Sy Sx};`

**Paper crosslink:**

- The paper’s primary transformation is **Iz → Ix** (Mz → Mx) with good phase behavior across offsets and RF inhomogeneity.
- The Spinach example encodes this and adds extra constraints:
  - `Sz → Sx` corresponds directly to Iz → Ix (desired excitation).
  - `Sy → Sy` enforces preservation of transverse magnetization along y, helping control phase and effective rotation axis.
  - `Sx → −Sz` enforces behavior consistent with a well‑defined effective rotation, improving overall rotation fidelity.
- Together, these state mappings push the optimized pulse toward a **high‑fidelity excitation with good phase properties**, paralleling the paper’s requirement that Mx ≈ M0 and phase deviations be small (<2–3° over most offsets/RF).

---

### 4. Drift and control operators

**Code:**

- Cartesian control operators and drift Hamiltonian:
  - `Lx=operator(spin_system,'Lx','13C');`
  - `Ly=operator(spin_system,'Ly','13C');`
  - `H=hamiltonian(assume(spin_system,'nmr'));`
- Packaging into the control structure:
  - `control.drifts={{H}};`
  - `control.operators={Lx,Ly};`

**Paper crosslink:**

- The paper’s Hamiltonian has:
  - A drift from resonance offset: Δω · Iz.
  - Control fields: ω1x(t) Ix + ω1y(t) Iy.
- In Spinach:
  - `H` (with the Zeeman terms from `inter.zeeman.scalar`) is the drift, which includes the frequency offsets.
  - `Lx` and `Ly` are the RF control operators corresponding to Ix and Iy.
- Thus, `control.drifts` and `control.operators` directly implement the same Hamiltonian structure as in the analytical optimal control derivation.

---

### 5. Ensemble over RF amplitudes (B1 inhomogeneity / miscalibration)

**Code:**

- Time discretization:
  - `control.pulse_dt=1e-6*ones(1,40);`  → 40 steps of 1 µs (40 µs total).
- RF amplitude ensemble:
  - `control.pwr_levels=2*pi*linspace(50e3,70e3,10);`
- Ensemble optimization:
  - `control.parallel='ensemble';`

**Paper crosslink:**

- The paper targets robust performance against RF inhomogeneity/miscalibration:
  - RF amplitude anywhere in 10–20 kHz, i.e. a wide range of B1 values.
  - Performance plots (e.g. Fig. 3–4) show excitation quality over both offset and RF amplitude.
- In the code:
  - `pulse_dt` discretizes the pulse into piecewise‑constant elements, matching the time‑discretized OC formulation.
  - `pwr_levels` defines multiple RF amplitudes (50–70 kHz here) used as different ensemble members.
  - `parallel='ensemble'` ensures the optimization considers all these RF amplitudes and offsets simultaneously, designing **one pulse that works across the whole ensemble**.
- This is the Spinach implementation of the paper’s **“calibration‑free” broadband excitation** design: robustness to substantial RF variations is enforced by optimizing over an explicit RF ensemble.

---

### 6. Cost functional and regularization (penalties)

**Code:**

- Optimization method and penalties:
  - `control.method='lbfgs';`
  - `control.penalties={'NS','SNS'};`
  - `control.p_weights=[0.01 10.0];`
  - `control.max_iter=200;`

**Paper crosslink:**

- The paper uses optimal control to maximize the overlap between final and target states over the ensemble, subject to constraints that enforce constant amplitude and realistic pulse shapes. This is expressed through a cost functional with both **fidelity** and **constraint/penalty** terms.
- In Spinach:
  - `lbfgs` corresponds to a quasi‑Newton gradient‑based optimizer, in the same family as the gradient‑based schemes derived in the paper and in related GRAPE literature.
  - Penalties `'NS'` and `'SNS'` (with weights `p_weights`) represent Spinach’s internal regularization terms:
    - They penalize excessive control power and/or roughness in the waveform.
    - They push the solution toward smooth, experimentally realizable pulses consistent with the constant‑amplitude / phase‑modulated structure discussed in the paper.
- The exact functional forms are implementation details, but conceptually they match the paper’s use of penalty terms to enforce pulse constraints while maximizing excitation fidelity.

---

### 7. GRAPE optimization step (the core optimal‑control algorithm)

**Code:**

- Set up the optimal control problem:
  - `spin_system=optimcon(spin_system,control);`
- Initial guess for the waveform:
  - `guess=(1/4)*ones(2,40);`
- GRAPE optimization:
  - `xy_profile=fminnewton(spin_system,@grape_xy,guess);`

**Paper crosslink:**

- The paper derives an optimal control algorithm tailored to constant‑amplitude, phase‑modulated pulses:
  - The pulse is discretized in time.
  - Gradients of the cost functional with respect to the control parameters are computed.
  - A gradient‑based iterative scheme (essentially GRAPE/Newton‑type) updates the pulse.
- In Spinach:
  - `optimcon` prepares all the forward/backward propagators and ensemble structures required for GRAPE.
  - `grape_xy` is the objective + gradient function for Cartesian x/y controls.
  - `fminnewton` is a Newton‑type optimizer that uses these gradients to refine the control.
- The main conceptual difference is parameterization:
  - The paper parameterizes the pulse via phase φ(t) at fixed amplitude A.
  - `magic_pulse_cart.m` parameterizes the two Cartesian components, CLx(t) and CLy(t), and enforces suitable behavior (smoothness, effective constant amplitude) via penalties rather than explicit polar coordinates.
  - Despite this, the underlying optimal control structure is the same as in the paper.

---

### 8. Extracting and applying the optimized pulse

**Code:**

- Scale and split the optimized waveform:
  - `xy_profile=mean(control.pwr_levels)*xy_profile;`
  - `CLx=xy_profile(1,:);`
  - `CLy=xy_profile(2,:);`
- Apply the optimized pulse:
  - `rho_init=state(spin_system,'Lz','13C');`
  - `rho=shaped_pulse_xy(spin_system,H,{Lx,Ly},{CLx,CLy}, ... 'expv-pwc');`

**Paper crosslink:**

- The paper analyzes the **final magnetization** after applying the optimized pulse, for many offsets and RF amplitudes (e.g. Mx > 0.99 M0, phase errors <2–3°).
- In the Spinach script:
  - `CLx` and `CLy` are the time‑dependent RF amplitudes along x and y, analogous to A cosφ(t) and A sinφ(t).
  - `shaped_pulse_xy` uses these waveforms to propagate the initial state under the time‑dependent Hamiltonian using an exponential integrator (`expv-pwc` for piecewise‑constant controls).
- This reproduces the “apply optimized pulse and inspect the result” step of the paper, but within Spinach’s simulation framework.

---

### 9. Building spectra and comparing with a hard pulse

**Code:**

- After the optimized pulse:
  - Sets acquisition parameters (`parameters.spins`, `rho0=rho`, `coil`, `sweep`, `npoints`, etc.).
  - Calls `fid=liquid(spin_system,@acquire,parameters,'nmr');`
  - Applies apodization and Fourier transform.
  - Plots the resulting 1D spectrum (`plot_1d`).
- Comparison hard pulse:
  - Sets `parameters.rho0=state(spin_system,'Lz','13C');`
  - Defines a conventional 90° hard pulse (`pulse_pwr`, `pulse_dur`, etc.).
  - Calls `fid=liquid(spin_system,@sp_acquire,parameters,'nmr');`
  - Processes and plots that spectrum in a separate subplot.

**Paper crosslink:**

- The paper compares the optimized PM‑BEBOP pulse to:
  - Phase‑corrected hard pulses and their excitation profiles.
  - It demonstrates superior bandwidth and RF‑tolerance of the optimized pulse (e.g., Fig. 4).
- `magic_pulse_cart.m` implements a similar comparison:
  - The first simulation shows the response to the optimized magic/PM‑BEBOP‑like pulse.
  - The second shows the response to a conventional hard pulse with comparable nominal flip angle.
  - Both are plotted over a sweep (`sweep=70000`) that covers the wide targeted bandwidth.
- This directly parallels the paper’s theoretical/experimental comparison of optimized vs. hard pulses, now done entirely in simulation.

---

### 10. Conceptual summary

- **Paper problem:**  
  Design a constant‑amplitude, phase‑modulated broadband excitation pulse that:
  - Robustly performs Iz → Ix over a wide resonance‑offset range.
  - Is insensitive to significant RF miscalibration / inhomogeneity.
  - Maintains good phase properties across the ensemble.

- **Spinach implementation (`magic_pulse_cart.m`):**
  - Uses Spinach to build a Liouville‑space representation of the same single‑spin control system:
    - Offsets are encoded via Zeeman terms in `H`.
    - RF controls are Cartesian operators `Lx`, `Ly`.
  - Encodes Iz → Ix and phase/rotation constraints via multiple initial/target states:
    - `{Sx, Sy, Sz} → {−Sz, Sy, Sx}`.
  - Models RF inhomogeneity as a discrete amplitude ensemble (`control.pwr_levels`) and optimizes a single pulse for all ensemble members.
  - Employs a GRAPE‑style optimal control algorithm (`optimcon` + `grape_xy` + `fminnewton`) analogous to the method derived in the paper.
  - Validates performance by simulating and plotting spectra for both the optimized pulse and a conventional hard pulse, mirroring the paper’s performance comparisons.

In this way, `magic_pulse_cart.m` serves as a concrete Spinach implementation of the 2006 PM‑BEBOP broadband excitation design, with every major element of the paper (offset ensemble, RF inhomogeneity, Iz → Ix target, optimal control search, and hard‑pulse comparison) represented explicitly in the code.

