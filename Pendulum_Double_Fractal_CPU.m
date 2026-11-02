clear
close all

% Simulation time
T = 8.9; % [s] 
tspan = linspace(0,T,4);  % SUPER  IMPORTANT FOR THE FINAL IMAGE IS HOW MANY EVALUATIONS OR INTGRATION STEPS ARE SUMMED UP FOR THE FINAL CHAOS METRIC. USE ONLY THE LAST VALUES TO HAVE MORE SHARP FEATURES!!! 
eps = 1e-5; % perturbation strength

% Pendulum design
m1 = 1;  % [kg]
m2 = 1;  % [kg]
L1 = 4;  % [m]
L2 = 0.8*4;  % [m]
g = 10;  % [m/s^2]
param = [m1;L1;m2;L2;g];

% Define grid resolution
n = 1001; % must be an ODD number so there will be a center line that is exactly at angle 0 in both dimensions

hyperparam_str = ['_T_',num2str(T),'_eps_',num2str(eps, '%.0e'),'_m1m2_',num2str(m1),'_L1L2_',num2str(L1),'_n_',num2str(n),'_CPU_ode45()_RelAbsTol_-1e6']


% Chose render mode
render_mode = 'free_zoom'; % 'symetry' or 'free_zoom'
switch render_mode
    case 'symetry'
        n_half = ceil(n/2); % simulate only on half of the first angle range since the other half is the same
        alpha1_lb =  0 * pi/180; 
        alpha1_ub =  180 * pi/180;
        alpha2_lb = -180 * pi/180;
        alpha2_ub =  180 * pi/180;
        alpha1_span = linspace(alpha1_lb, alpha1_ub, n_half);
        alpha2_span = linspace(alpha2_lb, alpha2_ub, n);
        Deviation_Mat = zeros(n, n_half);   % Pre-allocate output vector
    case 'free_zoom'
        n_half = n; % dirty hotfix
        alpha1_lb = -32  * pi/180;
        alpha1_ub =  85  * pi/180;
        alpha2_lb = -125 * pi/180;
        alpha2_ub =  35  * pi/180;

        alpha1_lb = -0  * pi/180;
        alpha1_ub =  180  * pi/180;
        alpha2_lb = -180 * pi/180;
        alpha2_ub =  0  * pi/180;
        alpha1_span = linspace(alpha1_lb, alpha1_ub, n);
        alpha2_span = linspace(alpha2_lb, alpha2_ub, n);
        Deviation_Mat = zeros(n,n);   % Pre-allocate output vector
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve ODE in ANGLE coordinates
odefun = @(t,x) pendulum_double_angle_minTrig_rhs(x,param);
% odevarfun = @(t,x_aug) pendulum_double_angle_aug_variational(x_aug, param);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 1);


tic
for i = 1:n
    % use the symetry of the system: Simulating the angles of the first link from 0 to pi is sufficient. We simply rotate the solution 
    alpha1_span = linspace(alpha1_lb, alpha1_ub, n_half); % This is here to avoid "broadcasting and communication overhead"
    % alpha1_span = linspace(0, ub, n_half);
    for j = 1:n_half
        % Perturb the angles
        alpha1_pos = alpha1_span(j) + eps;
        alpha2_pos = alpha2_span(i) + eps;

        alpha1_neg = alpha1_span(j) - eps;
        alpha2_neg = alpha2_span(i) - eps;

        % %%% Variational Approach
        % x0 = [alpha1_span(j); alpha2_span(i); 0; 0; 2*eps; 2*eps; 0; 0]; % Initial state for variational stuff
        % [~, y_pos] = ode45(odevarfun, tspan, x0, options);
        % delta_y = y_pos(:,5:8);
        % %%%

        % Form the vector of initial conditions
        x0_pos = [alpha1_pos; alpha2_pos; 0; 0];
        x0_neg = [alpha1_neg; alpha2_neg; 0; 0];

        % Solve the ODE's with the given angles
        [~, y_pos] = ode45(odefun, tspan, x0_pos, options);
        [~, y_neg] = ode45(odefun, tspan, x0_neg, options);

        % Compute difference between the two perturbed phase space trajectories
        % delta_y = y_pos(end,:) - y_neg(end,:);  % SHARP DETAILS IF ONLY THE FINAL STATE IS USED
        delta_y = y_pos - y_neg;

        % Deviation_Mat(i,j) = sum(sqrt(sum(delta_y.^2,1)),2);  % Sum of the column norms. Sum of the magnitude of the state difference vector
        Deviation_Mat(i,j) = sum(sqrt(sum(delta_y.^2,2)),1);  % Sum of phase space, row , vector norm. Faster than sum(vecnorm(delta_y,2,2)). Sum of the magnitude of the difference vector, [dx1_k,dx2_k,...,dxn_k], at each discrete time step over all difference vectors

    end
end
time_double_pendulum_ANGLE = toc
ms_per_ODE45_call = time_double_pendulum_ANGLE/(n * n * 2) * 1000

% Reshape to form the correct image
switch render_mode
    case 'symetry'
        % Construct the full image from one half
        Deviation_Mat = [flip(flip(Deviation_Mat(:,2:end),1),2),Deviation_Mat];
        % Measure elapsed time
        CPU_time_symetry_s = toc
        Gridpoint_time_ms = CPU_time_symetry_s/(n*n_half) * 1000
        alpha1_span = linspace(-pi/1, alpha1_ub, n);
    case 'free_zoom'
        % Measure elapsed time
        n
        CPU_time_free_zoom_s = toc
        Gridpoint_time_ms = CPU_time_free_zoom_s/(n^2) * 1000
end

%%% Prepare the plots
fig = figure('Position', [200, 100, 1200, 1200], 'Color','k', 'Name','Double Pendulum Fractal'); % For the dynamcis
ax = axes(fig, 'Position', [0.08, 0.08, 0.88, 0.88]);


% Convert radians to degrees
alpha1_span = alpha1_span .* 180./pi;
alpha2_span = alpha2_span .* 180./pi;
[A1, A2] = meshgrid(alpha1_span, alpha2_span);

% Plot
surf(A1,A2, log(Deviation_Mat), 'EdgeAlpha', 0)


xlabel('$\alpha_0 \; [^{\circ}]$', 'FontSize', 20, 'Color', 'w')
ylabel('$\beta_0 \;  [^{\circ}]$', 'FontSize', 20, 'Color', 'w')
zlabel(['Chaos level after ',num2str(T),'s'], 'FontSize', 20, 'Color', 'w')

xlim([alpha1_span(1), alpha1_span(end)])
ylim([alpha2_span(1), alpha2_span(end)])


% colormap(inferno)
% colormap(magma)
% colormap(viridis)
% colormap(twilight_shifted)
% colormap(hot)
% colormap(flip(hot))

colormap((slanCM('seasons')))

% colormap(flip(slanCM('fusion')))
% 

% colormap((slanCM('tokyo')))
% colormap(flip(slanCM('gist_heat')))

% colormap(slanCM('pride'))
% colormap((slanCM('copper2')))
% 
% 
%  colormap((slanCM('emergency')))
% 
% colormap(flip(slanCM('gist_ncar')))
% colormap(flip(slanCM('gnuplot')))

% ax.Position = [0, 0, 1, 1];
% ax.XTicks = []; 
% ax.YTicks = [];

ax.View = [0,90];


grid off
ax.Color = 'k';
ax.XAxis.Color = 'w';
ax.YAxis.Color = 'w';
ax.ZAxis.Color = 'w';

ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.ZAxis.FontSize = 14;

ax.View = [0,90];

set(findall(fig,'-property','Box'),'Box','off') % optional
set(findall(fig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,y] = pendulum_double_angle_RK4(x0, param, Nt, dt)

T = dt*Nt;
t = linspace(0, T, Nt+1);
% Preallocate solution matrix
y = zeros(length(x0), Nt+1);
y(:,1) = x0;

for k = 1:Nt
    Xpos = y(:,k);
    k1 = pendulum_double_angle_minTrig_rhs(Xpos, param);
    k2 = pendulum_double_angle_minTrig_rhs(Xpos + 0.5*dt*k1, param);
    k3 = pendulum_double_angle_minTrig_rhs(Xpos + 0.5*dt*k2, param);
    k4 = pendulum_double_angle_minTrig_rhs(Xpos + dt*k3, param);
    y(:,k+1) = Xpos + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

% Transpose to be compliant with an ode45() output
y = y'; % time incremetns pre row and states are in columns
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


