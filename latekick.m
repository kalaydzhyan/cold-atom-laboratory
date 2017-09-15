%% Releasing condensate from a trap
%% and applying delta-kick cooling within a long time interval

global glTime;                               
global gldeltaT;
global Timeline;
global now_stamp;
global color_limit;

% upper limit on the value of the color bar (drawing)
color_limit = 1.3E-5;

% creating folder for the output
now_stamp = datestr(datetime('now'),'yyyy-MM-dd-HH-mm-ss');
mkdir(now_stamp);

% Add path to the subroutines 
addpath(genpath(pwd));

f0 = 100;                                    % Units of frequency in Hz
omega0 = 2*pi*f0;
t0 = 1/omega0;                               % Units of time in s

% Initial trap parameters:
Omega_x = 1;                                 % Frequencies as measured in omega_0
Omega_y = 1;                                 
Alpha_x = 0;                                 % Cubic anharmonicity of the potential
Alpha_y = 0;

Natoms = 1000;                               % Number of atoms in the condensate
m0 = 87;                                     % Mass of atoms in atomic units (Rb = 87)
hbar = 1.05*10^(-34);                        % hbar in J*s
ma = 1.66*10^(-27);                          % Unit of atomic mass in kg
a0 = (10^6)*sqrt(hbar/(omega0*m0*ma));       % Unit of length in um
aBohr = 5.3*10^(-5);                         % Bohr radius in um
as = 104.5*aBohr;                            % Scattering length for the given species

Beta = 4*pi*as*Natoms/a0;                    % Parameter of nonlinearity in GPE
alpha0 = m0*ma*omega0^2/a0;                  % Unit of cubic nonlinearity in kg/(um*s^2)

%% Setting the geometry
Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 1e-2;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-4};
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
Lsize = 400;
xmin = -Lsize;
xmax = Lsize;
ymin = -Lsize;
ymax = Lsize;
Nx = 2^10+1;
Ny = 2^10+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);

%% Setting the physical problem
Delta = 0.5;
%Beta = 26; % N = 1000
mass = 1;

Physics2D = Physics2D_Var2d(Method, Delta, Beta);
Physics2D = Dispersion_Var2d(Method, Physics2D);
% Physics2D = Potential_Var2d(Method, Physics2D,...
%    @(X,Y) quadratic_plus_quartic_potential2d(2,1,1,0,X,Y));
Physics2D = Potential_Var2d(Method, Physics2D,...
    @(X,Y) quadratic_plus_cubic_potential2d(mass, Omega_x, Omega_y, Alpha_x, Alpha_y, X, Y));
Physics2D = Nonlinearity_Var2d(Method, Physics2D);


%% Setting the initial data
InitialData_Choice = 2;
Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_Choice);

%% Setting informations and outputs
%Outputs = OutputsINI_Var2d(Method);
%Printing = 0;
%Evo = 100;
%Draw = 1;
%Print = Print_Var2d(Printing,Evo,Draw);

%Physics2D = Potential_Var2d(Method, Physics2D,...
%    @(X,Y) quadratic_plus_cubic_potential2d(mass, 0, 0, Alpha_x, Alpha_y, X, Y));
%Physics2D = Nonlinearity_Var2d(Method, Physics2D);


%[Phi, Outputs] = GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);

%-----------------------------------------------------------
% Dynamic part
%-----------------------------------------------------------

Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
%Type = 'Relaxation';
Totalt = 0.110;                                    % Total simulation time, in s
dt = 1E-4;                                       % Time step of the simulation in s
dkkt = 0.1;                                      % Delta-kick cooling time, in s
dkkd = 0.0015;                                   % Delta-kick half-duration, in s
% the interval of time [TL, TR] corresponds to the free expansion and can
% be accelerated computationally by manking the time step larger
TL = floor(0.01/t0);
TR = floor((dkkt-4*dkkd)/t0);                    % within 4 STD from the DK-moment
DF = 10;                                          % factor, by which the time step for the simulation
                                                 % is increased at the free
                                                 % expansion evolution

Deltat = dt/t0;                                  % Time step of the simulation
Stop_time = TL;                     
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
Method.Splitting = 'Fourth';

Physics2D = Physics2D_Var2d(Method, Delta, Beta);
Physics2D = Dispersion_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D);

Outputs_iterations = 100;
Outputs_save = 1;                                 % Number of iterations prior to output
Outputs = OutputsINI_Var2d(Method, Outputs_iterations, Outputs_save);
Printing = 0;
Evo = 1;
Draw = 1; % KEEP ALWAYS 1
Print = Print_Var2d(Printing, Evo, Draw);

glTime = Deltat*Evo*t0;              % Global time counter
gldeltaT = Deltat*Evo*t0;            % time interval between outputs, in s
Timeline = [];
[Phi_1, Outputs] = GPELab2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, [], Print);

if Timeline(end) ~= TL*t0 
    Timeline = [Timeline TL*t0];
end

Deltat = DF*dt/t0;                                  % Time step of the simulation
Stop_time = TR-TL;                     
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
Method.Splitting = 'Fourth';

Physics2D = Physics2D_Var2d(Method, Delta, Beta);
Physics2D = Dispersion_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D);

glTime = TL*t0 + Deltat*Evo*t0;
gldeltaT = Deltat*Evo*t0;            % time interval between outputs, in s
[Phi_2, Outputs] = GPELab2d(Phi_1, Method, Geometry2D, Physics2D, Outputs, [], Print);

if Timeline(end) ~= TR*t0 
    Timeline = [Timeline TR*t0];
end

Deltat = dt/t0;                                  % Time step of the simulation
Stop_time = Totalt/t0 - TR;                     
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
Method.Splitting = 'Fourth';
% Delta-kick cooling trap parameters:
Omega_x = 0.098;                                   % set this and next to 0 for a potential turned off
Omega_y = 0.098;
Alpha_x = 0;      
Alpha_y = -2E-4;
Physics2D = Physics2D_Var2d(Method, Delta, Beta);
Physics2D = Dispersion_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D);
Physics2D = TimePotential_Var2d(Method, Physics2D,...
    @(T,X,Y) exp(-(T+TR-dkkt/t0)^2/(dkkd/t0)^2)*quadratic_plus_cubic_potential2d(mass, Omega_x, Omega_y, Alpha_x, Alpha_y, X, Y));

glTime = TR*t0 + Deltat*Evo*t0;
gldeltaT = Deltat*Evo*t0;            % time interval between outputs, in s
[Phi_3, Outputs] = GPELab2d(Phi_2, Method, Geometry2D, Physics2D, Outputs, [], Print);

if Timeline(end) ~= Totalt
    Timeline = [Timeline Totalt];
end

plot(Timeline, Outputs.x_rms{1}.*a0,'b-', Timeline, Outputs.y_rms{1}.*a0,'r-')
xlabel('Time, s')
ylabel('RMS radius, um')
title('Size of the condensate')
legend('x', 'y', 'Location','southeast')
print(strcat(now_stamp,'/expansion_plot.png'), '-dpng');

disp('Scaling factor: ')
disp(Outputs.x_rms{1}(length(Outputs.x_rms{1}))/Outputs.x_rms{1}(1))