function coherence_transfer_through_J_coup()
% Coherence transfer through scalar J-coupling in a two-spin system. This
% script implements, in Spinach, the first SIMPSON optimal control example
% from:
%
%   Z. Tošner, C. Kehlet, N. Khaneja, S.J. Glaser, N.C. Nielsen, "Optimal
%   control in NMR spectroscopy: Numerical implementation in
%    SIMPSON", J. Magn. Reson. 197 (2009) 120–134.
%
%                 http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% Physical problem: A heteronuclear two-spin system (1H–13C) with an
% isotropic scalar coupling J, with both nuclei set on resonance (no
% chemical shift offsets). The goal is to transfer transverse magnetisation
% from the proton to the carbon:
%
%       I_x(1H)  →  I_x(13C)
%
% over a fixed evolution period T = 1/J.
%
% Authors:
%   Aditya Dev  (aditya.dev@weizmann.ac.il) Ilya Kuprov
%   (ilya.kuprov@weizmann.ac.il)

%==========================================================================
% 1. Spin system specification
%==========================================================================

% Scalar coupling strength, Hz
J = 140;

% Total evolution time: T = 1 / J (standard for J-transfer tests), s
T = 1 / J;

% Number of time slices in the control sequence
n_t_steps  = 300;
dt = T / n_t_steps;

% Static magnetic field, Tesla
sys.magnet = 14.1;

% Two-spin heteronuclear system: 1H–13C
sys.isotopes = {'1H','13C'};

% Isotropic chemical shifts, ppm (both on resonance)
inter.zeeman.scalar = {0.0, 0.0};

% Isotropic scalar coupling, Hz
inter.coupling.scalar = cell(2,2);
inter.coupling.scalar{1,2} = J;

% Basis set: full Liouville space in spherical tensor formalism
bas.formalism = 'sphten-liouv';
bas.approximation = 'none';

% Spinach housekeeping
spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

%==========================================================================
% 2. Initial and target states
%==========================================================================

% Initial state: Lx on spin 1 (1H)
rho_init = state(spin_system, {'Lx'}, {1});
rho_init = rho_init / norm(full(rho_init), 2);

% Target state: Lx on spin 2 (13C)
rho_targ = state(spin_system, {'Lx'}, {2});
rho_targ = rho_targ / norm(full(rho_targ), 2);

%==========================================================================
% 3. Drift and control operators
%==========================================================================

% RF control operators on 1H (spin 1)
LxH = operator(spin_system, 'Lx', 1);
LyH = operator(spin_system, 'Ly', 1);

% RF control operators on 13C (spin 2)
LxC = operator(spin_system, 'Lx', 2);
LyC = operator(spin_system, 'Ly', 2);

% Offset operators (z-terms) for both spins
LzH = operator(spin_system, 'Lz', 1);
LzC = operator(spin_system, 'Lz', 2);

% Scalar J-coupling operator: I1·I2 = Lx1Lx2 + Ly1Ly2 + Lz1Lz2
JOp = operator(spin_system, {'Lx','Lx'}, {1,2}) + ...
    operator(spin_system, {'Ly','Ly'}, {1,2}) + ...
    operator(spin_system, {'Lz','Lz'}, {1,2});

% Drift Hamiltonian (liquid-state NMR context)
H = hamiltonian(assume(spin_system, 'nmr'));

%==========================================================================
% 4. Optimal control setup
%==========================================================================

%-----------------------------
% 4.1 Drift and control fields
%-----------------------------

% Single drift Hamiltonian
control.drifts    = {{H}};

% Four Cartesian control operators: x/y on 1H and 13C
control.operators = {LxH, LyH, LxC, LyC};

control.off_ops  = {LzH, LzC, JOp};
control.offsets  = {[-10 0 10],[-10 0 10],[-10 0 10]};

%-----------------------------
% 4.3 State-to-state objective
%-----------------------------

control.rho_init = {rho_init};  % starting state
control.rho_targ = {rho_targ};  % destination state

%-----------------------------
% 4.4 Time grid and RF levels
%-----------------------------

% Piecewise-constant time grid (uniform slices)
control.pulse_dt = dt * ones(1, n_t_steps);

% RF power levels (ensemble over B1 mis-calibrations), rad/s Nominal RF
% nutation frequencies 1–1.5 kHz (example values).
control.pwr_levels = 2*pi * linspace(1000, 1500, 3);

% Constant amplitude profile for each RF channel; grape_phase will optimise
% the phase only. Here we have two RF channels (1H and 13C), hence 2 ×
% n_t_steps amplitudes.
control.amplitudes = ones(2, n_t_steps);

% Optimisation method (LBFGS-GRAPE in phase parametrisation)
control.method    = 'lbfgs';

% Termination condition: maximum number of iterations
control.max_iter  = 200;

% Parallelisation over ensemble members, if available
control.parallel  = 'ensemble';

% Plotting options during optimisation
control.plotting = {'correlation_order', ...
    'coherence_order',   ...
    'local_each_spin',   ...
    'spectrogram',       ...
    'phi_controls',      ...
    'xy_controls'};

%==========================================================================
% 5. Initial guess for the phase waveform
%==========================================================================

% For grape_phase, the optimisation variable is the phase profile. There
% are two RF channels (1H and 13C), hence two phase rows.
n_channels = numel(control.operators) / 2;      % x/y per channel

% Random initial phase profile, radians
guess = (pi/3) * rand(n_channels, n_t_steps);

%==========================================================================
% 6. Spinach optimal control housekeeping and optimisation
%==========================================================================

% Register control settings with Spinach
spin_system = optimcon(spin_system, control);

% Run LBFGS GRAPE phase optimisation
phase_profile = fminnewton(spin_system, @grape_phase, guess); 

end
