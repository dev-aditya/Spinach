function broadband_refocusing_pulses()
% Broadband refocusing (π) pulse for liquid-state 1H NMR.
%
% Spinach implementation of the broadband refocusing example from:
%
%   Z. Tošner, C. Kehlet, N. Khaneja, S.J. Glaser, N.C. Nielsen, J. Magn.
%   Reson. 197 (2009) 120–134. http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% Task: design a 200 µs broadband x-phase π pulse
%       {Sx -> Sx, Sy -> -Sy, Sz -> -Sz}
% over an offset range of ±12.5 kHz with RF penalty.
%
% Authors:
%   Aditya Dev  (aditya.dev@weizmann.ac.il) Ilya Kuprov
%   (ilya.kuprov@weizmann.ac.il)

%---------------------------------------------------------------------------%
% Spin system
%---------------------------------------------------------------------------%

T          = 200e-6;                % pulse length, s
n_t_steps  = 600;                   % time slices
dt         = T / n_t_steps;

sys.magnet   = 14.1;               % Tesla
sys.isotopes = {'1H'};
inter.zeeman.scalar = {0.0};       % on-resonance

bas.formalism     = 'sphten-liouv';
bas.approximation = 'IK-2';
bas.space_level   = 1;
bas.connectivity  = 'scalar_couplings';

spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

%---------------------------------------------------------------------------%
% States and operators
%---------------------------------------------------------------------------%

Sx = state(spin_system, 'Lx', 1);
Sy = state(spin_system, 'Ly', 1);
Sz = state(spin_system, 'Lz', 1);

Sx = Sx / norm(full(Sx), 2);
Sy = Sy / norm(full(Sy), 2);
Sz = Sz / norm(full(Sz), 2);

% Multi-target refocusing: Rx(pi)
rho_init = {Sx,  Sy,  Sz};
rho_targ = {Sx, -Sy, -Sz};

Lx = operator(spin_system, 'Lx', 1);
Ly = operator(spin_system, 'Ly', 1);
Lz = operator(spin_system, 'Lz', 1);

H  = hamiltonian(assume(spin_system, 'nmr'));

%---------------------------------------------------------------------------%
% Optimal control settings
%---------------------------------------------------------------------------%

control.drifts    = {{H}};
control.operators = {Lx, Ly};

% Offset ensemble over ±12.5 kHz
control.off_ops   = {Lz};
control.offsets   = {linspace(-12.5e3, 12.5e3, 10)};  % Hz

control.rho_init  = rho_init;
control.rho_targ  = rho_targ;

control.pulse_dt  = dt * ones(1, n_t_steps);

% RF scale (max ~15 kHz)
control.pwr_levels = 2*pi*15e3;

% RF power penalty
control.penalties  = {'NS'};
control.p_weights  = 0.01;

control.method     = 'lbfgs';
control.max_iter   = 200;
control.parallel   = 'ensemble';

control.plotting   = {'phi_controls', ...
    'xy_controls', ...
    'amp_controls', ...
    'level_populations'};

%---------------------------------------------------------------------------%
% Initial guess and optimisation
%---------------------------------------------------------------------------%

n_channels = numel(control.operators);
guess      = (pi/3) * rand(n_channels, n_t_steps);

% Set a flat initial amplitude profile for both channels
control.amplitudes = ones(n_channels, n_t_steps);

spin_system = optimcon(spin_system, control);
xy_profile  = fminnewton(spin_system, @grape_xy, guess); %#ok<NASGU>

end
