function spin_traj_analysys_fig1()

sys.magnet=28.18;

sys.isotopes={'1H'};
inter.zeeman.scalar={0};


bas.formalism='sphten-liouv';
bas.approximation='none'; 

spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up spin states
Sx=state(spin_system,'Lx','1H');
Sy=state(spin_system,'Ly','1H');
Sz=state(spin_system,'Lz','1H');
Sx=Sx/norm(full(Sx),2);
Sy=Sy/norm(full(Sy),2);
Sz=Sz/norm(full(Sz),2);

% Get the control operators
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');

% Get the offset operator
Lz=operator(spin_system,'Lz','1H');

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));
% Define control parameters
control.drifts={{H}};                           % Drift
control.operators={Lx,Ly};                      % Controls
control.rho_init={Sz};                   % Starting states
control.rho_targ={Sx};                   % Target states
control.pulse_dt=1e-3*ones(1,625)/625;              % Pulse interval grid
control.pwr_levels=2*pi*15e3*linspace(0.7, 1.3, 10); % Power levels 15kHz iwth 30 % inhomogenity
control.off_ops={Lz};
control.offsets={linspace(-25e3,25e3,100)};
control.amplitudes=ones(1,625);                  % Amplitude profile
control.method='lbfgs';                         % Optimisation method
control.max_iter=200;                           % Termination condition
control.parallel='ensemble';                    % Parallelisation

% Plotting options
control.plotting={'phi_controls','spectrogram'};

% Initial guess
guess=(pi/9)*randn(1,625);

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Run LBFGS GRAPE pulse optimisation
phi_profile=fminnewton(spin_system,@grape_phase,guess);

end