clear
close all



% Pendulum design
string_length = 12; % [m]
g = 10; % [m/s^2]
m = 0.1;  % [kg]
alpha = 95 * pi/180;
Rot_Mat = [cos(alpha),-sin(alpha); sin(alpha),cos(alpha)]; % counter-clockwise rotation
pos_bob_init = Rot_Mat*[0;-string_length];

F_g_mag = m*g;
F_g = [0; -F_g_mag]; % gravity acts only on vertical dimension
F_string_mag = F_g_mag * cos(0);
F_string = [F_string_mag * sin(0); F_string_mag * cos(0)];
F_net = F_string + F_g;

% Plot
fig = figure('Position', [500, 200, 800, 800]);
ax = axes(fig);
hold on
grid on

% Plot first the possible orrbit path
x_circ = linspace(-string_length, string_length, 200);
y_circ1 = sqrt(string_length^2 - x_circ.^2);
y_circ2 = -sqrt(string_length^2 - x_circ.^2);
plot(x_circ, y_circ1, 'Color', [1,0.3,0],'LineStyle', '--')
plot(x_circ, y_circ2, 'Color', [1,0.3,0],'LineStyle', '--')

% The anker point pivot
scatter(0,0, 100, 'k', 'filled')
% The pendulums string
string_plt = plot([0,pos_bob_init(1)],[0,pos_bob_init(2)], 'Color', 'k', 'LineWidth', 2);
% The pendulum bob
bob_plt = scatter(pos_bob_init(1),pos_bob_init(2), 400, 'red', 'filled');
% The force of gravity acting on the bob
% Fg_plt = quiver(pos_bob_init(1),pos_bob_init(2),0,-F_g_mag);
% The force on the bob due to the string constraining its motion
% F_string_plt = quiver(pos_bob_init(1),pos_bob_init(2),F_string(1),F_string(2));
% The net force on the pendulum
% F_net_plt = quiver(pos_bob_init(1),pos_bob_init(2),F_net(1),F_net(2));

xlim([-1.2 * string_length, 1.2* string_length])
ylim([-1.4 * string_length, 1.2* string_length])
ax.DataAspectRatio = [1,1,1];

% for i = 45:720
% % Define which angel should be plotted
% alpha = i * pi/180;
% Rot_Mat = [cos(alpha),-sin(alpha); sin(alpha),cos(alpha)]; % counter-clockwise rotation
% 
% 
% % Compute tension force in the cord
% F_string_mag = (F_g_mag * cos(alpha));
% F_string = -[F_string_mag * sin(alpha); F_string_mag * -cos(alpha)];
% 
% F_net = F_string + F_g;
% 
% % [-m*g*cos(alpha)*sin(alpha);m*g*(cos(alpha)^2 - 1)]
% 
% % Rotate the bob around the angle
% pos_bob = Rot_Mat * pos_bob_init;
% 
% pos_x = pos_bob(1);
% pos_y = pos_bob(2);
% 
% 
% % Update Plots
% string_plt.XData = [0,pos_bob(1)];
% string_plt.YData = [0,pos_bob(2)];
% bob_plt.XData = pos_bob(1);
% bob_plt.YData = pos_bob(2);
% 
% Fg_plt.XData = pos_bob(1);
% Fg_plt.YData = pos_bob(2);
% 
% F_string_plt.XData = pos_bob(1);
% F_string_plt.YData = pos_bob(2);
% F_string_plt.UData = F_string(1);
% F_string_plt.VData = F_string(2);
% 
% F_net_plt.XData = pos_bob(1);
% F_net_plt.YData = pos_bob(2);
% F_net_plt.UData = F_net(1);
% F_net_plt.VData = F_net(2);
% 
% pause(1/60)
% 
% end




%%% Solve ODE in cartesian coordinates
odefun = @(t,x) pendulum_cartesian_ODE(t,x,[string_length,g,m]);
% odefun = @(t,x) pendulum_cartesian_rhs(x,[m;string_length;g]);
x0 = [pos_bob_init;0;0];
T = 10;
tspan = linspace(0,T,30*T);
tic 
[~, y] = ode15s(odefun,tspan,x0);
toc

trace_plt = plot(0,0,'Color','k', 'LineStyle', ':', 'LineWidth', 2);

frame_rate = 1/60;
t0 = tic;
for i = 1:length(tspan)

    frame_start = toc(t0);
    pos_bob = y(i,[1,2]);

    % Update Plots
    string_plt.XData = [0,pos_bob(1)];
    string_plt.YData = [0,pos_bob(2)];
    bob_plt.XData = pos_bob(1);
    bob_plt.YData = pos_bob(2);

    % Fg_plt.XData = pos_bob(1);
    % Fg_plt.YData = pos_bob(2);

    if i > 50
        trace_plt.XData = y(i-50:i, 1);
        trace_plt.YData = y(i-50:i, 2);
    else
        trace_plt.XData = y(1:i, 1);
        trace_plt.YData = y(1:i, 2);        
    end

    % Try to achieve stable FPS
    drawnow limitrate
    next_frame = i * frame_rate;
    sleep_time = next_frame - frame_start;
    if sleep_time > 0
        pause(sleep_time)
    end
end



% %%% Solve the ODE in angle space
% odefun = @(t,x) pendulum_angle_ODE(t,x,[string_length,g,m]);
% x0 = [alpha;0];
% T = 10;
% tspan = linspace(0,T,30*T);
% tic 
% [~, y] = ode15s(odefun,tspan,x0);
% toc
% trace_plt = plot(0,0,'Color','k', 'LineStyle', ':', 'LineWidth', 2);
% frame_rate = 1/60;
% t0 = tic;
% for i = 1:length(tspan)
%     frame_start = toc(t0);
%     angle = y(i,1);
%     pos_bob(1) = string_length*sin(angle);
%     pos_bob(2) = -string_length*cos(angle);
% 
%     % Update Plots
%     string_plt.XData = [0,pos_bob(1)];
%     string_plt.YData = [0,pos_bob(2)];
%     bob_plt.XData = pos_bob(1);
%     bob_plt.YData = pos_bob(2);
% 
%     if i > 50
%         trace_plt.XData = string_length*sin(y(i-50:i, 1));
%         trace_plt.YData = -string_length*cos(y(i-50:i, 1));
%     else
%         trace_plt.XData = string_length*sin(y(1:i, 1));
%         trace_plt.YData = -string_length*cos(y(1:i, 1));        
%     end
% 
%     % Try to achieve stable FPS
%     drawnow limitrate
%     next_frame = i * frame_rate;
%     sleep_time = next_frame - frame_start;
%     if sleep_time > 0
%         pause(sleep_time)
%     end
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxdt = pendulum_cartesian_ODE(t,x,theta)
% Unpack state variables
pos_x = x(1);
pos_y = x(2);
vel_x = x(3);
vel_y = x(4);
% Model parameters
L = theta(1);
g = theta(2);
m = theta(3);

% Helper variables
vel_sq = vel_x^2 + vel_y^2;
lambda = (- m * vel_sq + m * g * pos_y) / L^2;

% ODE's
dxdt = zeros(4,1);
dxdt(1) = vel_x;
dxdt(2) = vel_y;
dxdt(3) = lambda/m*pos_x;
dxdt(4) = -g + lambda/m*pos_y;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxdt = pendulum_angle_ODE(t,x,theta)
% Unpack state variables
angle = x(1);
angle_vel = x(2);

% Model parameters
L = theta(1);
g = theta(2);

% ODE's
dxdt = zeros(2,1);
dxdt(1) = angle_vel;
dxdt(2) = -g/L*sin(angle);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%