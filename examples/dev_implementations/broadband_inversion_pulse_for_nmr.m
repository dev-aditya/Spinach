function broadband_inversion_pulse_for_nmr()
% Trying to reproduce the SIMPSON example 2

T = 600e-6; % mu s
n_t_steps = 600;
dt = T/n_t_steps;
sys.magnet=14.1; % Does it get affected by the offets values used below? 

sys.isotopes={'1H'};
inter.zeeman.scalar={0.0}; % No zeeman shifts  i.e on resonance

bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up spin states
Sx=state(spin_system,'Lx',1);
Sy=state(spin_system,'Ly',1);
Sz=state(spin_system,'Lz',1);
Sx=Sx/norm(full(Sx),2);
Sy=Sy/norm(full(Sy),2);
Sz=Sz/norm(full(Sz),2);
% Get the control operators
LxH=operator(spin_system,'Lx',1);
LyH=operator(spin_system,'Ly',1);


% Get offset operators
LzH=operator(spin_system,'Lz',1);

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                           % Drift
control.operators={LxH,LyH};              % Controls
control.off_ops={LzH};
control.offsets={linspace(-50e3, 50e3, 10)}; % offsets in Hz
control.rho_init={Sz};                      % Starting state
control.rho_targ={-Sz};                      % Destination state
control.pulse_dt = dt * ones(1, n_t_steps);      % or column vector (nsteps×1)
control.penalties={'NS'};
control.p_weights = 0.01;
control.pwr_levels=2*pi*30e3;       % rms power level
control.amplitudes=ones(1, n_t_steps);
control.method='lbfgs';                         % Optimisation method
control.max_iter=200;                           % Termination condition
control.parallel='ensemble';                    % Parallelisation

% Plotting options
%control.plotting={'xy_controls','amp_controls','robustness',};

% Initial guess
guess=rand(numel(control.operators),n_t_steps)*pi/4;

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Run LBFGS GRAPE pulse optimisation
xy_profile=fminnewton(spin_system,@grape_xy,guess);

end