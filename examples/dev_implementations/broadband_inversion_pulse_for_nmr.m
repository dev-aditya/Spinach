function broadband_inversion_pulse_for_nmr()
% Broadband inversion pulse design for liquid-state NMR.
%
% This script reproduces, using Spinach, the second SIMPSON optimal control
% example from:
%
%   Z. Tošner, C. Kehlet, N. Khaneja, S.J. Glaser, N.C. Nielsen, "Optimal
%   control in NMR spectroscopy: Numerical implementation
%    in SIMPSON", J. Magn. Reson. 197 (2009) 120–134.
%
%                 http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% Physical problem: A single-spin (1H) system is considered in the rotating
% frame, on resonance (no isotropic chemical shift). The goal is to design
% a broadband inversion pulse that performs:
%
%       I_z  →  -I_z
%
% uniformly over a frequency offset range of ±50 kHz, under a fixed pulse
% duration T = 600 µs, discretised into 600 time steps of 1 µs.
%
% The control fields are Cartesian RF components (Lx, Ly) on 1H, and
% robustness is enforced by an offset ensemble and an RF-power penalty:
%
%   - Offset ensemble: a set of resonance offsets over the desired
%     ±50 kHz window.
%   - RF constraint: an 'NS' penalty that suppresses excessive RF
%     power while allowing enough freedom to achieve inversion.
% Authors:
%   Aditya Dev  (aditya.dev@weizmann.ac.il) Ilya Kuprov
%   (ilya.kuprov@weizmann.ac.il)

%==========================================================================
% 1. Spin system specification
%==========================================================================

% Total pulse duration (s) and time discretisation
T          = 600e-6;                % 600 µs total duration
n_t_steps  = 600;                   % 600 time steps
dt         = T / n_t_steps;         % 1 µs per step

% Static magnetic field (Tesla)
sys.magnet = 14.1;

% Single 1H spin at the origin
sys.isotopes = {'1H'};

% Isotropic chemical shift (ppm) – set to zero, i.e. on resonance
inter.zeeman.scalar = {0.0};

% Basis: spherical tensor Liouville space with IK-2 approximation IK-2
% keeps a complete basis on each spin but neglects higher multi-spin
% orders, which is sufficient for this one-spin problem.
bas.formalism     = 'sphten-liouv';
bas.approximation = 'IK-2';
bas.space_level   = 1;
bas.connectivity  = 'scalar_couplings';

% Spinach housekeeping
spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

%==========================================================================
% 2. Initial and target states
%==========================================================================

% Cartesian spin operators as density operators
Sx = state(spin_system, 'Lx', 1);
Sy = state(spin_system, 'Ly', 1);
Sz = state(spin_system, 'Lz', 1);

% Normalise the states (Hilbert–Schmidt norm)
Sx = Sx / norm(full(Sx), 2);
Sy = Sy / norm(full(Sy), 2);
Sz = Sz / norm(full(Sz), 2);

% State-to-state transfer: Sz  →  -Sz (longitudinal inversion)
rho_init = Sz;
rho_targ = -Sz;

%==========================================================================
% 3. Drift and control operators
%==========================================================================

% RF control operators (Cartesian components) on 1H
LxH = operator(spin_system, 'Lx', 1);
LyH = operator(spin_system, 'Ly', 1);

% Offset operator (z-component) for 1H, used to generate the offset
% ensemble
LzH = operator(spin_system, 'Lz', 1);

% Drift Hamiltonian (liquid-state NMR context)
H = hamiltonian(assume(spin_system, 'nmr'));

%==========================================================================
% 4. Optimal control setup
%==========================================================================

%-----------------------------
% 4.1 Drift and control fields
%-----------------------------

% Single drift Liouvillian
control.drifts    = {{H}};

% Control operators: x and y RF components on 1H
control.operators = {LxH, LyH};

%-----------------------------
% 4.2 Offset ensemble
%-----------------------------
% The list below should span the desired ±50 kHz range with appropriate
% sampling density.

control.off_ops   = {LzH};
control.offsets   = {2*pi*linspace(-50e3, 50e3, 10)};   % offsets (rads/s)

%-----------------------------
% 4.3 State-to-state objective
%-----------------------------

control.rho_init  = {rho_init};    % starting state Sz
control.rho_targ  = {rho_targ};    % destination state -Sz

%-----------------------------
% 4.4 Time grid and RF levels
%-----------------------------

% Piecewise-constant time grid (uniform slices)
control.pulse_dt  = dt * ones(1, n_t_steps);

% RF power specification: Here we choose a single RMS power level
% corresponding to 30 kHz.
control.pwr_levels = 2*pi * 30e3;  % rad/s

% Amplitude profile: grape_xy will optimise Cartesian amplitudes directly.
% There are two RF channels (LxH, LyH), hence 2 × n_t_steps amplitudes.
n_channels         = numel(control.operators);
control.amplitudes = ones(n_channels, n_t_steps);

%-----------------------------
% 4.5 Penalties and optimisation settings
%-----------------------------

% RF power penalty ('NS'): suppresses excessive RF power while still
% allowing the optimiser to reach high-fidelity inversion.
control.penalties  = {'NS'};
control.p_weights  = 0.01;

% Optimisation method
control.method     = 'lbfgs';

% Termination condition: maximum number of iterations
control.max_iter   = 200;

% Parallelisation over the offset ensemble, if supported
control.parallel   = 'ensemble';

% Plotting options during optimisation
control.plotting   = {'xy_controls', ...
    'amp_controls', ...
    'robustness'};

%==========================================================================
% 5. Initial guess for the RF waveform
%==========================================================================

% For grape_xy, the optimisation variable is the Cartesian RF amplitude on
% each control channel. There are two channels (LxH, LyH).
n_channels = numel(control.operators);

% Random initial Cartesian amplitude profile (in radians), scaled to a
% conservative initial range.
guess = (pi/4) * rand(n_channels, n_t_steps);

%==========================================================================
% 6. Spinach optimal control housekeeping and optimisation
%==========================================================================

% Register control settings with Spinach
spin_system = optimcon(spin_system, control);

% Run LBFGS GRAPE Cartesian optimisation
xy_profile = fminnewton(spin_system, @grape_xy, guess); %#ok<NASGU>

end
