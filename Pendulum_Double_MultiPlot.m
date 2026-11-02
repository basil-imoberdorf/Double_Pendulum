clear
close all

% Simulation time
T = 20;
frame_rate = 1/24; % images per second
tspan = linspace(0,T,T/frame_rate);

% Pendulum design
g = 10;  % [m/s^2]
m1 = 1;  % [kg]
m2 = 1;  % [kg]
L1 = 5; % [m]
L2 = 5; % [m]
theta = [m1;L1;m2;L2;g];

alpha1 = 150 * pi/180;
alpha2 = 150 * pi/180;

Rot_Mat1 = [cos(alpha1),-sin(alpha1); sin(alpha1),cos(alpha1)]; % counter-clockwise rotation
Rot_Mat2 = [cos(alpha2),-sin(alpha2); sin(alpha2),cos(alpha2)]; % counter-clockwise rotation
pos_bob1_init = Rot_Mat1*[0;-L1];
pos_bob2_init = Rot_Mat2*[0;-L2] + pos_bob1_init;

trace_length = 40;
etrace_length = 300;

%%% Prepare the plots
fig = figure('Position', [200, 200, 1000, 1000], 'Color','k', 'ToolBar','none', 'Name', 'Double Pendulum Dynamics'); % For the dynamcis
ax = axes(fig, 'Position', [0, 0, 1, 1]);
hold on
grid on

n = 15;

alpha1_span = linspace(pi/2 - 0.055, pi/2-0.052, n);
alpha2_span = linspace(-pi/2 - 0.055, -pi/2-0.052, n);


plot_handles = cell(n,1);
for i = 1:n
    Rot_Mat1 = [cos(alpha1_span(i)),-sin(alpha1_span(i)); sin(alpha1_span(i)),cos(alpha1_span(i))]; % counter-clockwise rotation
    Rot_Mat2 = [cos(alpha2_span(i)),-sin(alpha2_span(i)); sin(alpha2_span(i)),cos(alpha2_span(i))]; % counter-clockwise rotation
    pos_bob1_init = Rot_Mat1*[0;-L1];
    pos_bob2_init = Rot_Mat2*[0;-L2] + pos_bob1_init;

    [lnk1_plt, lnk2_plt, trace_plt, bob1_plt, bob2_plt] = double_pendulum_plot_handles(ax, pos_bob1_init);
    plot_handles{i} = [lnk1_plt, lnk2_plt, trace_plt, bob1_plt, bob2_plt];

    max_range = L1 + L2;
    xlim([-1.1 * max_range, 1.1* max_range])
    ylim([-1.1 * max_range, 1.1* max_range])
end

ax.DataAspectRatio = [1,1,1];
ax.Color = 'k';
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
ax.GridColor = 'y';
ax.MinorGridColor = 'y';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.Toolbar.Visible = 'off';

% % Figure for tracking the total, kinetic and potential energy  
% fig = figure('Position', [1250, 500, 800, 400], 'Name', 'Double Pendulum Energy'); % For the dynamcis
% ay = axes(fig);
% hold on
% grid on
% 
% % Plot Energy
% E_pot1 = m1*g*(pos_bob1_init(2) + L1);  % We want energy to be zero at rest   
% E_pot2 = m2*g*(pos_bob2_init(2) + L1 + L2);  % We want energy to be zero at rest
% E_pot = E_pot1 + E_pot2;
% E_tot_fix =  E_pot;
% plot(ay, [tspan(1),tspan(end)],[E_tot_fix,E_tot_fix], 'LineWidth', 1, 'Color', [1,0,0]); % True Total Energy
% 
% plot_handles_energy = cell(n,1);
% for i = 1:n
%     [Etot_plt, Epot1_plt, Ekin1_plt, Epot2_plt, Ekin2_plt] = double_pendulum_Energy_plot_handles(ay, '-');
%     plot_handles_energy{i} = [Etot_plt, Epot1_plt, Ekin1_plt, Epot2_plt, Ekin2_plt];
% end
% 
% xlabel('Simulation time  [s]')
% ylabel('Energy  [J]')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve ODE in ANGLE coordinates
odefun = @(t,x) pendulum_double_angle_rhs(x,theta);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'MaxStep', 0.01); 

tic 
y = cell(n,1);
for i = 1:n
    x0 = [alpha1_span(i); alpha2_span(i); 0; 0];
    [~, y{i}] = ode45(odefun, tspan, x0, options);

    % alpha1 = y{i}(:, 1);
    % alpha2 = y{i}(:, 2);
    % alpha1_d = y{i}(:, 3);
    % alpha2_d = y{i}(:, 4);
    % 
    % % Compute the position of the bob's
    % bob1_x =  L1*sin(alpha1);
    % bob1_y = -L1*cos(alpha1);
    % bob2_x =  L2*sin(alpha2) + bob1_x;
    % bob2_y = -L2*cos(alpha2) + bob1_y;
    % 
    % bob1_x_d =  L1*cos(alpha1).*alpha1_d;
    % bob1_y_d =  L1*sin(alpha1).*alpha1_d;
    % bob2_x_d =  L2*cos(alpha2).*alpha2_d + bob1_x_d;
    % bob2_y_d =  L2*sin(alpha2).*alpha2_d + bob1_y_d;
    % 
    % % Compute energy in the system
    % E_pot1 = m1*g*(bob1_y + L1);  % We want energy to be zero at rest
    % E_pot2 = m2*g*(bob2_y + L1 + L2);  % We want energy to be zero at rest
    % E_kin1 = 0.5*m1*(bob1_x_d.^2 + bob1_y_d.^2);
    % E_kin2 = 0.5*m2*(bob2_x_d.^2 + bob2_y_d.^2);
    % E_kin = E_kin1 + E_kin2;
    % E_pot = E_pot1 + E_pot2;
    % E_tot = E_kin + E_pot;
    % E = [E_tot, E_pot1, E_pot2, E_kin1, E_kin2];
end
time_double_pendulum_ANGLE = toc
time_double_pendulum_ANGLE/n * 10000

% pause(1)
t0 = tic;
timerr = tic;
for i = 1:length(tspan)
    % Start measuring time, used to keep constant FPS
    frame_start = toc(t0);

    for k = 1:n
        update_angle_plots(y{k}, i, L1,L2, trace_length, plot_handles{k})
        % update_energy_plots(E, tspan, i, etrace_length, plot_handles_energy{k})
    end

    % if i > 1
    %     ay.XLim = [min(Etot_plt.XData), max(Etot_plt.XData)];
    % end

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
function [lnk1_plt, lnk2_plt, trace_plt, bob1_plt, bob2_plt] = double_pendulum_plot_handles(ax, pos_bob)
% The pendulum links
lnk1_plt = plot(ax, [0,pos_bob(1)],[0,pos_bob(2)], 'Color', 'w', 'LineWidth', 2);
lnk2_plt = plot(ax, [pos_bob(1),pos_bob(1)],[pos_bob(2),pos_bob(2)], 'Color', 'w', 'LineWidth', 1.4);
% The anker point pivot
scatter(ax, 0,0, 100, 'magenta', 'filled')
% Outermost bob Trace
trace_plt = plot(ax, 0,0,'Color', [0,0.8,0.8], 'LineStyle', ':', 'LineWidth', 2);
% The pendulum bob
bob1_plt = scatter(pos_bob(1),pos_bob(2), 500, [0,0.8,0], 'filled', 'MarkerEdgeColor','w');
bob2_plt = scatter(pos_bob(1),pos_bob(2), 500, [0,0.9,0.9], 'filled', 'MarkerEdgeColor','w');
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
function update_angle_plots(y_data, i, L1,L2, trace_length, plot_handles)
    % Unpack plot handles
    lnk1_plt  = plot_handles(1);
    lnk2_plt  = plot_handles(2);
    trace_plt = plot_handles(3);
    bob1_plt  = plot_handles(4);
    bob2_plt  = plot_handles(5);
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
    % if i > trace_length
    %     trace_plt.XData =  L2*sin(y_data(i-trace_length:i, 2)) +  L1*sin(y_data(i-trace_length:i, 1));
    %     trace_plt.YData = -L2*cos(y_data(i-trace_length:i, 2)) + -L1*cos(y_data(i-trace_length:i, 1));
    % else
    %     trace_plt.XData =  L2*sin(y_data(1:i, 2)) +  L1*sin(y_data(1:i, 1));
    %     trace_plt.YData = -L2*cos(y_data(1:i, 2)) + -L1*cos(y_data(1:i, 1));
    % end
    bob1_plt.XData = bob1_x;
    bob1_plt.YData = bob1_y;
    bob2_plt.XData = bob2_x;
    bob2_plt.YData = bob2_y;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_energy_plots(E, tspan, i, etrace_length, plot_handles_energy)
% Unpack plot handles
Etot_plt = plot_handles_energy(1);
Epot1_plt = plot_handles_energy(2);
Ekin1_plt = plot_handles_energy(3);
Epot2_plt = plot_handles_energy(4);
Ekin2_plt = plot_handles_energy(5);

E_tot  = E(:,1);
E_pot1 = E(:,2);
E_pot2 = E(:,3);
E_kin1 = E(:,4);
E_kin2 = E(:,5);

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



































