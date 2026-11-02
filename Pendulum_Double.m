clear
close all

% Simulation time
T = 20;
frame_rate = 1/60; % images per second
tspan = linspace(0,T,T/frame_rate);

% Pendulum design
g = 10;  % [m/s^2]
m1 = 1;  % [kg]
m2 = 1;  % [kg]
L1 = 1; % [m]
L2 = 1; % [m]
param = [m1;L1;m2;L2;g];


% % Complex wenn L1/L2 = 2 und m1=m2
% alpha1 =  142.02 * pi/180;
% alpha2 = -137.52* pi/180;

% % Die Amphore wenn L1/L2 = 1/2 und m1=m2
% alpha1 = 51.84 * pi/180;
% alpha2 = -173.88 * pi/180;

%%% Coole Winkel für m1=m2=1, L1=L2=1, g=10
% % Der Akrobat
% alpha1 = 0.001 + 2.331;
% alpha2 = 0.001 - 3.00917;

% % Die Vase
% alpha1 = 0.0 + 3.056;  
% alpha2 = 0.0 + 0.178;  

% % Der Schmetterling
% alpha1 = 0.0 + 2.9012;
% alpha2 = 0.0 + 1.91265;

% % Der andere Schmetterling
% alpha1 = 0.0 + 1.435;
% alpha2 = 0.0 + 2.633;

% % Der dicke Schmetterling
% alpha1 = 0.0 + 2.224;  
% alpha2 = 0.0 - 2.113;

% % Der Kreis
% alpha1 = 0.00 + 2.51952;
% alpha2 = 0.00 + 2.3137;

% % Die Halfpipe
% alpha1 = 0.0 + 1.761;
% alpha2 = 0.0 + 2.009;

% % Pseudostabil
% alpha1 = 0.0 + 2.85374;  % spezial island
% alpha2 = 0.0 + 2.13326;

% % Das Herz
% alpha1 = 0.00 + 2.9388;
% alpha2 = 0.00 + 2.34198;

% % Die Schlaufe
% alpha1 = 0.00 - 0.547556;
% alpha2 = 0.00 - 2.84001;

% % Der Weirdo
% alpha1 = 0.00 + 2.23895;
% alpha2 = 0.00 - 0.470031;

% % Der Weirdo II
% alpha1 = 0.00 + 2.06378;
% alpha2 = 0.00 - 0.479697;

% Die Schere - Must see
alpha1 = 0.0 + 2.30545;
alpha2 = 0.0 - 1.44392;

% % Der labile Nachbar der Schere
% alpha1 = 0.0 + 2.32567;
% alpha2 = 0.0 - 1.41002;

% % Die Biene
% alpha1 = 0.0 + 2.72594;
% alpha2 = 0.0 - 1.84508;

% % Der Hummel
% alpha1 = 0.0 + 2.22691;
% alpha2 = 0.0 - 2.11212;

% % Der Hummel
% alpha1 =90.48* pi/180;
% alpha2 =-20.88* pi/180;


% % Der Hummel
% alpha1 = 0.0 + 2.30593;
% alpha2 = 0.0 - 1.40365;




Rot_Mat1 = [cos(alpha1),-sin(alpha1); sin(alpha1),cos(alpha1)]; % counter-clockwise rotation
Rot_Mat2 = [cos(alpha2),-sin(alpha2); sin(alpha2),cos(alpha2)]; % counter-clockwise rotation
pos_bob1_init = Rot_Mat1*[0;-L1];
pos_bob2_init = Rot_Mat2*[0;-L2] + pos_bob1_init;

trace_length = 400;
etrace_length = 300;

%%% Prepare the plots
fig = figure('Position', [200, 200, 1000, 1000], 'Color','k', 'ToolBar','none', 'Name', 'Double Pendulum Dynamics'); % For the dynamcis
ax = axes(fig, 'Position', [0, 0, 1, 1]);
hold on
grid on

[lnk1_plt, lnk2_plt, trace_plt, bob1_plt, bob2_plt] = double_pendulum_plot_handles(ax, pos_bob1_init, pos_bob2_init);
[lnk1_plt2, lnk2_plt2, trace_plt2, bob1_plt2, bob2_plt2] = double_pendulum_plot_handles(ax, pos_bob1_init, pos_bob2_init);

bob1_plt2.MarkerFaceColor = [1,0,0];
bob2_plt2.MarkerFaceColor = [1,0.5,0];
trace_plt2.Color = [1,0.7,0];

max_range = L1 + L2;
xlim([-1.1 * max_range, 1.1* max_range]); ax.Toolbar.Visible = 'off';
ylim([-1.1 * max_range, 1.1* max_range]); ax.DataAspectRatio = [1,1,1];
ax.Color = 'k';
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
ax.GridColor = 'y';
ax.MinorGridColor = 'y';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';


% Figure for tracking the total, kinetic and potential energy  
fig = figure('Position', [1250, 500, 800, 400], 'Name', 'Double Pendulum Energy'); % For the dynamcis
ay = axes(fig);
hold on
grid on

% Plot Energy
E_pot1 = m1*g*(pos_bob1_init(2) + L1);  % We want energy to be zero at rest   
E_pot2 = m2*g*(pos_bob2_init(2) + L1 + L2);  % We want energy to be zero at rest
E_pot = E_pot1 + E_pot2;
E_tot_fix =  E_pot;
plot(ay, [tspan(1),tspan(end)],[E_tot_fix,E_tot_fix], 'LineWidth', 1, 'Color', [1,0,0]); % True Total Energy

[Etot_plt, Epot1_plt, Ekin1_plt, Epot2_plt, Ekin2_plt] = double_pendulum_Energy_plot_handles(ay, '-');
[Etot_plt2, Epot1_plt2, Ekin1_plt2, Epot2_plt2, Ekin2_plt2] = double_pendulum_Energy_plot_handles(ay, '--');

xlabel('Simulation time  [s]')
ylabel('Energy  [J]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve ODE in cartesian coordinates
% odefun = @(t,x) pendulum_double_cartesian_rhs(x,param);
% x0 = [pos_bob1_init; pos_bob2_init; 0; 0; 0; 0];
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); 
% tic 
% [~, y1] = ode45(odefun,tspan,x0, options);
% time_double_pendulum_CARTESIAN = toc
% 
% % Compute energy in the system
% E_kin1 = 0.5*m1*(y1(:,5).^2 + y1(:,6).^2);
% E_pot1 = m1*g*(y1(:,2) + L1);  % We want energy to be zero at rest   
% E_kin2 = 0.5*m2*(y1(:,7).^2 + y1(:,8).^2);
% E_pot2 = m2*g*(y1(:,4) + L1 + L2);  % We want energy to be zero at rest
% E_kin = E_kin1 + E_kin2;
% E_pot = E_pot1 + E_pot2;
% E_tot = E_kin + E_pot;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve ODE in ANGLE coordinates
odefun = @(t,x) pendulum_double_angle_minTrig_rhs(x,param);
x0 = [alpha1; alpha2; 0; 0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); 
tic 
[~, y1] = ode45(odefun,tspan,x0, options);
time_double_pendulum_ANGLE = toc

% Compute the position of the bob's
bob1_y = -L1*cos(y1(:, 1));
bob2_y = -L2*cos(y1(:, 2)) + bob1_y;

bob1_x_d =  L1*cos(y1(:, 1)).*y1(:, 3);
bob1_y_d =  L1*sin(y1(:, 1)).*y1(:, 3);
bob2_x_d =  L2*cos(y1(:, 2)).*y1(:, 4) + bob1_x_d;
bob2_y_d =  L2*sin(y1(:, 2)).*y1(:, 4) + bob1_y_d;

% Compute energy in the system
E_pot1 = m1*g*(bob1_y + L1);  % We want energy to be zero at rest   
E_pot2 = m2*g*(bob2_y + L1 + L2);  % We want energy to be zero at rest
E_kin1 = 0.5*m1*(bob1_x_d.^2 + bob1_y_d.^2);
E_kin2 = 0.5*m2*(bob2_x_d.^2 + bob2_y_d.^2);
E_kin = E_kin1 + E_kin2;
E_pot = E_pot1 + E_pot2;
E_tot = E_kin + E_pot;


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve ODE in ANGLE coordinates
J = @(t,x) pendulum_double_angle_J(t, x, param);  % Jacbian ∂f/∂x
x0 = [alpha1; alpha2; 0; 0];

odefun = @(t,x) pendulum_double_angle_minTrig_rhs(x,param);
% odefun = @(t,x) pendulum_double_angle_inertiamat_rhs(x,param);

% odefun = @(t,x_aug) pendulum_double_angle_aug_variational(x_aug, param);
% x0 = [alpha1; alpha2; 0; 0; 1e-5; 1e-5; 0; 0]; % Initial state for variational stuff

% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); 
% tic 
% [~, y2] = ode45(odefun,tspan,x0, options);
% y2 = y2(:,1:4);
% time_double_pendulum_ANGLE_2 = toc
% 
% time_pc = (time_double_pendulum_ANGLE-time_double_pendulum_ANGLE_2)/time_double_pendulum_ANGLE * 100
% time_ratio_pc = time_double_pendulum_ANGLE/time_double_pendulum_ANGLE_2 
% 

%%% Runge-Kutta 4
dt = 0.01;
Nt = ceil(T/dt);
tic
[t, y2] = pendulum_double_angle_RK4(x0, param, Nt, dt);
time_double_pendulum_ANGLE_RK4 = toc
y2 = interp1(t,y2,tspan, 'linear');
time_ratio_pc = time_double_pendulum_ANGLE/time_double_pendulum_ANGLE_RK4 

% %%% Velocity Verlet
% dt = 0.001;
% Nt = T/dt;
% tic
% [t, y2] = velocity_verlet_angle_integ(x0, param, Nt, dt);
% y2 = interp1(t,y2,tspan, 'linear');
% time_double_pendulum_ANGLE_verlet = toc

% Compute the position of the bob's
bob1_x =  L1*sin(y2(:, 1));
bob1_y = -L1*cos(y2(:, 1));
bob2_x =  L2*sin(y2(:, 2)) + bob1_x;
bob2_y = -L2*cos(y2(:, 2)) + bob1_y;

bob1_x_d =  L1*cos(y2(:, 1)).*y2(:, 3);
bob1_y_d =  L1*sin(y2(:, 1)).*y2(:, 3);
bob2_x_d =  L2*cos(y2(:, 2)).*y2(:, 4) + bob1_x_d;
bob2_y_d =  L2*sin(y2(:, 2)).*y2(:, 4) + bob1_y_d;

% Compute energy in the system
E_pot1_a = m1*g*(bob1_y + L1);  % We want energy to be zero at rest   
E_pot2_a = m2*g*(bob2_y + L1 + L2);  % We want energy to be zero at rest
E_kin1_a = 0.5*m1*(bob1_x_d.^2 + bob1_y_d.^2);
E_kin2_a = 0.5*m2*(bob2_x_d.^2 + bob2_y_d.^2);
E_kin_a = E_kin1_a + E_kin2_a;
E_pot_a = E_pot1_a + E_pot2_a;
E_tot_a = E_kin_a + E_pot_a;

% odefun = @(t,x) pendulum_double_cartesian_rhs(x,param);
% x0 = [pos_bob1_init; pos_bob2_init; 0; 0; 0; 0];
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); 
% tic 
% [~, y2] = ode15s(odefun,tspan,x0, options);
% time_double_pendulum_CARTESIAN_ode15s = toc
% 
% % Compute energy in the system
% E_kin1_a = 0.5*m1*(y2(:,5).^2 + y2(:,6).^2);
% E_pot1_a = m1*g*(y2(:,2) + L1);  % We want energy to be zero at rest   
% E_kin2_a = 0.5*m2*(y2(:,7).^2 + y2(:,8).^2);
% E_pot2_a = m2*g*(y2(:,4) + L1 + L2);  % We want energy to be zero at rest
% E_kin_a = E_kin1_a + E_kin2_a;
% E_pot_a = E_pot1_a + E_pot2_a;
% E_tot_a = E_kin_a + E_pot_a;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause(1.5)  % wait until the plots are up before starting the animation
%%%%
t0 = tic;
timerr = tic;
for i = 1:length(tspan)
    % Start measuring time, used to keep constant FPS
    frame_start = toc(t0);

    % Plot positions
    % update_angle_plots(y2, i, L1,L2, trace_length, lnk1_plt,lnk2_plt,trace_plt,bob1_plt,bob2_plt)
    % update_carte_plots(y1, i, trace_length, lnk1_plt2,lnk2_plt2,trace_plt2,bob1_plt2,bob2_plt2)

    update_angle_plots(y1, i, L1,L2, trace_length, lnk1_plt,lnk2_plt,trace_plt,bob1_plt,bob2_plt)
    update_angle_plots(y2, i, L1,L2, trace_length, lnk1_plt2,lnk2_plt2,trace_plt2,bob1_plt2,bob2_plt2)

    % update_carte_plots(y1, i, trace_length, lnk1_plt,lnk2_plt,trace_plt,bob1_plt,bob2_plt)
    % update_carte_plots(y2, i, trace_length, lnk1_plt2,lnk2_plt2,trace_plt2,bob1_plt2,bob2_plt2)

    % Plot Energy
    update_energy_plots(tspan, i, E_pot1, E_pot2, E_kin1, E_kin2, E_tot, etrace_length, Epot1_plt, Epot2_plt, Ekin1_plt, Ekin2_plt, Etot_plt)
    update_energy_plots(tspan, i, E_pot1_a, E_pot2_a, E_kin1_a, E_kin2_a, E_tot_a, etrace_length, Epot1_plt2, Epot2_plt2, Ekin1_plt2, Ekin2_plt2, Etot_plt2)

    if i > 1
        ay.XLim = [min(Etot_plt.XData), max(Etot_plt.XData)];
    end

    %%% Try to achieve stable FPS
    drawnow limitrate
    next_frame = i * frame_rate;
    sleep_time = next_frame - frame_start;
    if sleep_time > 0
        pause(sleep_time)
    end
end
toc(timerr)







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lnk1_plt, lnk2_plt, trace_plt, bob1_plt, bob2_plt] = double_pendulum_plot_handles(ax, pos_bob1, pos_bob2)
% The pendulum links
lnk1_plt = plot(ax, [0,pos_bob1(1)],[0,pos_bob1(2)], 'Color', 'w', 'LineWidth', 2);
lnk2_plt = plot(ax, [pos_bob1(1),pos_bob2(1)],[pos_bob1(2),pos_bob2(2)], 'Color', 'w', 'LineWidth', 2);
% The anker point pivot
scatter(ax, 0,0, 100, 'magenta', 'filled')
% Outermost bob Trace
trace_plt = plot(ax, 0,0,'Color', [0,0.8,0.8], 'LineStyle', ':', 'LineWidth', 2);
% The pendulum bob
bob1_plt = scatter(pos_bob1(1),pos_bob1(2), 600, 'green', 'filled', 'MarkerEdgeColor','w');
bob2_plt = scatter(pos_bob2(1),pos_bob2(2), 600, [0,1,1], 'filled', 'MarkerEdgeColor','w');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Etot_plt, Epot1_plt, Ekin1_plt, Epot2_plt, Ekin2_plt] = double_pendulum_Energy_plot_handles(ax, linestyle)
Etot_plt  = plot(ax, [0,0], [0,0], 'LineWidth', 2, 'Color', [0.4,0.4,0.4],'LineStyle', linestyle);
Epot1_plt = plot(ax, [0,0], [0,0], 'LineWidth', 2, 'Color', [0,0.4,1],    'LineStyle', linestyle);
Ekin1_plt = plot(ax, [0,0], [0,0], 'LineWidth', 2, 'Color', [1,0.4,0],    'LineStyle', linestyle);
Epot2_plt = plot(ax, [0,0], [0,0], 'LineWidth', 2, 'Color', [0,0.6,0.7],  'LineStyle', linestyle);
Ekin2_plt = plot(ax, [0,0], [0,0], 'LineWidth', 2, 'Color', [0.8,0.5,0],  'LineStyle', linestyle);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_angle_plots(y_data, i, L1,L2, trace_length, lnk1_plt,lnk2_plt,trace_plt,bob1_plt,bob2_plt)
    % Extract the angles
    alpha1 = y_data(i, 1);
    alpha2 = y_data(i, 2);
    % Compute the position of the bob's
    bob1_x =  L1*sin(alpha1);
    bob1_y = -L1*cos(alpha1);
    bob2_x =  L2*sin(alpha2) + bob1_x;
    bob2_y = -L2*cos(alpha2) + bob1_y;
    % Update Plots
    lnk1_plt.XData = [0,bob1_x];
    lnk1_plt.YData = [0,bob1_y];
    lnk2_plt.XData = [bob1_x,bob2_x];
    lnk2_plt.YData = [bob1_y,bob2_y];
    if i > trace_length
        trace_plt.XData =  L2*sin(y_data(i-trace_length:i, 2)) +  L1*sin(y_data(i-trace_length:i, 1));
        trace_plt.YData = -L2*cos(y_data(i-trace_length:i, 2)) + -L1*cos(y_data(i-trace_length:i, 1));
    else
        trace_plt.XData =  L2*sin(y_data(1:i, 2)) +  L1*sin(y_data(1:i, 1));
        trace_plt.YData = -L2*cos(y_data(1:i, 2)) + -L1*cos(y_data(1:i, 1));
    end
    bob1_plt.XData = bob1_x;
    bob1_plt.YData = bob1_y;
    bob2_plt.XData = bob2_x;
    bob2_plt.YData = bob2_y;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_carte_plots(y_data, i, trace_length, lnk1_plt,lnk2_plt,trace_plt,bob1_plt,bob2_plt)
    % Extract Positions
    pos_bobs = y_data(i, 1:4);
    % Update Plots
    lnk1_plt.XData = [0,pos_bobs(1)];
    lnk1_plt.YData = [0,pos_bobs(2)];
    lnk2_plt.XData = [pos_bobs(1),pos_bobs(3)];
    lnk2_plt.YData = [pos_bobs(2),pos_bobs(4)];

    if i > trace_length
        trace_plt.XData = y_data(i-trace_length:i, 3);
        trace_plt.YData = y_data(i-trace_length:i, 4);
    else
        trace_plt.XData = y_data(1:i, 3);
        trace_plt.YData = y_data(1:i, 4); 
    end

    bob1_plt.XData = pos_bobs(1);
    bob1_plt.YData = pos_bobs(2);
    bob2_plt.XData = pos_bobs(3);
    bob2_plt.YData = pos_bobs(4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_energy_plots(tspan, i, E_pot1, E_pot2, E_kin1, E_kin2, E_tot, etrace_length, Epot1_plt, Epot2_plt, Ekin1_plt, Ekin2_plt, Etot_plt)
% Plot Energy
if i > etrace_length
    % Epot1_plt.XData = tspan(i-etrace_length:i);
    % Ekin1_plt.XData = tspan(i-etrace_length:i);
    % Epot2_plt.XData = tspan(i-etrace_length:i);
    % Ekin2_plt.XData = tspan(i-etrace_length:i);
    % % YData
    % Epot1_plt.YData = E_pot1(i-etrace_length:i);
    % Ekin1_plt.YData = E_kin1(i-etrace_length:i);
    % Epot2_plt.YData = E_pot2(i-etrace_length:i);
    % Ekin2_plt.YData = E_kin2(i-etrace_length:i);

    % Total Energy
    Etot_plt.XData  = tspan(i-etrace_length:i);
    Etot_plt.YData  = E_tot(i-etrace_length:i);
else
    % Epot1_plt.XData = tspan(1:i);
    % Ekin1_plt.XData = tspan(1:i);
    % Epot2_plt.XData = tspan(1:i);
    % Ekin2_plt.XData = tspan(1:i);
    % % YData
    % Epot1_plt.YData = E_pot1(1:i);
    % Ekin1_plt.YData = E_kin1(1:i);
    % Epot2_plt.YData = E_pot2(1:i);
    % Ekin2_plt.YData = E_kin2(1:i);

    % Total Energy
    Etot_plt.XData  = tspan(1:i);
    Etot_plt.YData  = E_tot(1:i);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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



% Verlet is not for systems where acceleration depends on velocity, only position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, y] = velocity_verlet_angle_integ(x0, param, Nt, dt)

    T = dt*Nt;
    t = linspace(0, T, Nt+1);

    % Preallocate solution matrix
    y = zeros(Nt+1, length(x0));
    y(1,:) = x0;

    % Get initial acceleration
    a_curr = pendulum_double_angle_rhs(x0, param);
    a_curr1 = a_curr(3);
    a_curr2 = a_curr(4);

    % dt_sqrt_half = 0.5*dt^2;
    % dt_half = 0.5*dt;
    
    for k = 1:Nt
        % Current state
        x1 = y(k, 1);   v1 = y(k, 3); 
        x2 = y(k, 2);   v2 = y(k, 4);

        % Update positions
        x1_new = x1 + v1*dt + 0.5*a_curr1*dt^2;
        x2_new = x2 + v2*dt + 0.5*a_curr2*dt^2;

        % % First velocity estimate (half step)
        % v1_half = v1 + 0.5*a_curr1*dt;
        % v2_half = v2 + 0.5*a_curr2*dt;
        % 
        % % Compute acceleration at new position with half-step velocity
        % x_half = [x1_new; x2_new; v1_half; v2_half];
        % a_next_full = pendulum_double_angle_rhs(x_half, param);
        % a_next1 = a_next_full(3);
        % a_next2 = a_next_full(4);
        % 
        % % Complete velocity update
        % v1_new = v1_half + 0.5*a_next1*dt;
        % v2_new = v2_half + 0.5*a_next2*dt;

        % Compute acceleration at new position with half-step velocity
        x_half = [x1_new; x2_new; v1; v2];
        a_next_full = pendulum_double_angle_rhs(x_half, param);
        a_next1 = a_next_full(3);
        a_next2 = a_next_full(4);

        % Complete velocity update
        v1_new = v1 + 0.5*(a_next1 + a_curr1)*dt;
        v2_new = v2 + 0.5*(a_next2 + a_curr2)*dt;

        % Update acceleration for next step
        a_curr1 = a_next1;
        a_curr2 = a_next2;

        % Store results
        y(k+1, :) = [x1_new, x2_new, v1_new, v2_new];
         
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%