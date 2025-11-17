% A template file for the "magic pulse" optimisations. The term refers to
% a family of broadband NMR pulses that are tolerant to resonance offsets
% and power calibration errors:
%
%               http://dx.doi.org/10.1016/j.jmr.2005.12.010
%
% Consider a 13C 90-degree excitation pulse in a 28.18 Tesla magnet. The
% pulse must uniformly excite a bandwidth of around 200 ppm (60 kHz) and
% must be short enough for the worst-case 13C-1H J-coupling (ca. 200 Hz)
% to be negligible. The latter requirement caps the duration at 1/100*J
% = 50 us. The pulse must accomplish the following transfers: {Lz -> Lx,
% Ly -> Ly, Lx -> -Lz}. A realistically achievable nutation frequency is
% between 50 kHz and 70 kHz across the RF coil.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il
% david.goodwin@inano.au.dk

function magic_pulse_phase()

% Set the magnetic field
sys.magnet=28.18;

% Single 13C spin at the origin 
sys.isotopes={'13C'};
inter.zeeman.scalar={0};

% Select a basis set - IK-2 keeps complete basis on each 
% spin in this case, but ignores multi-spin orders
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up spin states
Sx=state(spin_system,'Lx','13C');
Sy=state(spin_system,'Ly','13C');
Sz=state(spin_system,'Lz','13C');
Sx=Sx/norm(full(Sx),2);
Sy=Sy/norm(full(Sy),2);
Sz=Sz/norm(full(Sz),2);

% Get the control operators
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');

% Get the offset operator
Lz=operator(spin_system,'Lz','13C');

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                           % Drift
control.operators={Lx,Ly};                      % Controls
control.rho_init={ Sx Sy Sz};                   % Starting states
control.rho_targ={-Sz Sy Sx};                   % Target states
control.pulse_dt=2.5e-5*ones(1,100);              % Pulse interval grid
control.pwr_levels=2*pi*linspace(0.9e3,1.1e3,10); % Power levels
control.off_ops={Lz};
control.offsets={linspace(-30e3,30e3,100)};
control.amplitudes=ones(1,100);                  % Amplitude profile
control.method='lbfgs';                         % Optimisation method
control.max_iter=200;                           % Termination condition
control.parallel='ensemble';                    % Parallelisation

% Plotting options
control.plotting={'phi_controls','xy_controls',...
                  'robustness','spectrogram'};

% Initial guess
guess=(pi/5)*randn(1,100);

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Composite fidelity functional
function [traj_data,fidelity,gradient]=grape_bulu(phi_profile,spin_system)

    % Set the fidelity
    spin_system.control.fidelity='square';

    % First call: Sz -> Sx
    spin_system.control.rho_init={Sz}; 
    spin_system.control.rho_targ={Sx};
    [~,fidelity_a,gradient_a]=grape_phase(phi_profile,spin_system);

    % Second call: Sz -> Sy
    spin_system.control.rho_init={Sz}; 
    spin_system.control.rho_targ={Sy};
    [~,fidelity_b,gradient_b]=grape_phase(phi_profile,spin_system);

    % Composite fidelity and gradient
    fidelity=fidelity_a+fidelity_b;
    gradient=gradient_a+gradient_b;

    % Empty trajectory data
    traj_data=[];

end

% Run LBFGS GRAPE pulse optimisation
phi_profile=fminnewton(spin_system,@grape_bulu,guess);

% % Get Cartesian components of the pulse
% amp_profile=mean(control.pwr_levels)*control.amplitudes;
% [CLx,CLy]=polar2cartesian(amp_profile,phi_profile);
% 
% % Simulate the optimised pulse
% rho_init=state(spin_system,'Lz','13C');
% rho=shaped_pulse_xy(spin_system,H,{Lx,Ly},{CLx,CLy},...
%                     control.pulse_dt,rho_init,'expv-pwc');
% 
% % Set acquisition parameters
% parameters.spins={'13C'};
% parameters.rho0=rho;
% parameters.coil=state(spin_system,'L+','13C');
% parameters.decouple={};
% parameters.offset=0;
% parameters.sweep=70000;
% parameters.npoints=2048;
% parameters.zerofill=16384;
% parameters.axis_units='ppm';
% parameters.invert_axis=1;
% 
% % Simulate the free induction decay
% fid=liquid(spin_system,@acquire,parameters,'nmr');
% 
% % Apodisation
% fid=apodisation(spin_system,fid,{{'gauss',10}});
% 
% % Fourier transform
% spectrum=fftshift(fft(fid,parameters.zerofill));
% 
% % Plotting
% figure(2); subplot(2,1,2);
% plot_1d(spin_system,real(spectrum),parameters);
% kylabel('intensity, a.u.');
% 
% % Simulate the conventional hard pulse
% parameters.rho0=state(spin_system,'Lz','13C');
% parameters.pulse_frq=0;
% parameters.pulse_phi=pi/2;
% parameters.pulse_pwr=2*pi*60e3;
% parameters.pulse_dur=4.2e-6;
% parameters.pulse_rnk=3;
% parameters.method='expv';
% fid=liquid(spin_system,@sp_acquire,parameters,'nmr');
% fid=apodisation(spin_system,fid,{{'gauss',10}});
% spectrum=fftshift(fft(fid,parameters.zerofill));
% subplot(2,1,1); plot_1d(spin_system,real(spectrum),parameters);
% kylabel('intensity, a.u.');

end

