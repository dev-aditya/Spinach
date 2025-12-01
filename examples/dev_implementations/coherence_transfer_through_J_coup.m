function coherence_transfer_through_J_coup()
% Trying to reproduce the SIMPSON example 1

J = 140; % Hz
T =1/J; %s
n_t_steps = 150;
dt = T/n_t_steps;
sys.magnet=14.1;

sys.isotopes={'1H' '13C'};
inter.zeeman.scalar={0.0 0.0}; % No zeeman shifts  i.e on resonance
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=J;

bas.formalism='sphten-liouv';
bas.approximation='none'; 

spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up spin states for target and init
rho_init=state(spin_system, {'Lx'}, {1});
rho_init=rho_init/norm(full(rho_init),2);

rho_targ=state(spin_system,{'Lx'},{2});
rho_targ=rho_targ/norm(full(rho_targ),2);

% Get the control operators
LxH=operator(spin_system,'Lx',1);
LyH=operator(spin_system,'Ly',1);
LxC=operator(spin_system,'Lx',2);
LyC=operator(spin_system,'Ly',2);

% Get offset operators
LzH=operator(spin_system,'Lz',1);
LzC=operator(spin_system,'Lz',2);

% Get coupling operator
JOp=operator(spin_system,{'Lx','Lx'},{1 2})+...
  operator(spin_system,{'Ly','Ly'},{1 2})+...
  operator(spin_system,{'Lz','Lz'},{1 2});

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                           % Drift
control.operators={LxH,LyH,LxC,LyC};              % Controls
control.off_ops={LzH,LzC,JOp};
control.offsets={[-10 0 10],[-10 0 10],[-10 0 10]};
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pulse_dt = dt * ones(1, n_t_steps);      % or column vector (nsteps×1)
control.pwr_levels=2*pi*[1000 1200 1500];       % Power levels from 1000 - 1500 Hz
control.method='lbfgs';                         % Optimisation method
control.max_iter=200;                           % Termination condition
control.parallel='ensemble';                    % Parallelisation

% Plotting options
control.plotting={'correlation_order','coherence_order',...
                  'xy_controls','local_each_spin',...
                  'spectrogram'};

% Initial guess
guess=rand(numel(control.operators),n_t_steps)/4;

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Run LBFGS GRAPE pulse optimisation
xy_profile=fminnewton(spin_system,@grape_xy,guess);

end