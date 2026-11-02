clear
close all
reset(gpuDevice)


% Simulation set-up
T = 90; % [s]
dt = 0.005; % [s] RK4 step size
Nt = ceil(T/dt);
eps = 1e-5; % perturbation strength

% Pendulum design
g = 10;  % [m/s^2]
m1 = 1;  % [kg]
m2 = 1;  % [kg]
L1 = 1;  % [m]
L2 = 1;  % [m]
param = [m1;L1;m2;L2;g];

%%% Define grid resolution
n = 101; % must be an ODD number when used with render mode 'symetry' so there will be a uniqe line of pixels along the vertical axis

% Switch render mode
render_mode = 'symetry'; % 'symetry' or 'free_zoom'
switch render_mode
    % Use rotational symetry to cut down cost by half. Instead of n^2 points we compute only 1/2*n^2
    case 'symetry'
        n_half = ceil(n/2); % we will simulate only on half of the first angle range since the other half is the same
        alpha1_lb =  0 * pi/180; 
        alpha1_ub =  180 * pi/180;
        alpha2_lb = -180 * pi/180;
        alpha2_ub =  180 * pi/180;
        alpha1_span = linspace(alpha1_lb, alpha1_ub, n_half);
        alpha2_span = linspace(alpha2_lb, alpha2_ub, n);
        Deviation_flat = zeros(1, n*n_half);   % Pre-allocate output vector

    case 'free_zoom'
        alpha1_lb = -180 * pi/180;
        alpha1_ub = -10 * pi/180;
        alpha2_lb = -120 * pi/180;
        alpha2_ub =  120 * pi/180;
        alpha1_span = linspace(alpha1_lb, alpha1_ub, n);
        alpha2_span = linspace(alpha2_lb, alpha2_ub, n);
        Deviation_flat = zeros(1, n^2);   % Pre-allocate output vector
end

[A1, A2] = meshgrid(alpha1_span, alpha2_span);
% Flatten and transpose the grid
alpha1 = A1(:).';
alpha2 = A2(:).';
Np = numel(alpha1); % Number double pendulums simulated in parallel

% Create matrix of intial condition. Each column represents the state of a double pendulum
X_aug0 = [A1(:).'; A2(:).'; zeros(2,Np); eps.*ones(2,Np); zeros(2,Np)];

%%% Move the whole stack to the GPU
tic
X_aug = gpuArray(X_aug0);
Deviation_flat = gpuArray(Deviation_flat);

% Perform Runge-Kutta 4 integration on the GPU
for k = 1:Nt
    % RK4 
    K1 = pendulum_double_angle_aug_variational_vectorized(X_aug,               param);
    K2 = pendulum_double_angle_aug_variational_vectorized(X_aug + 0.5.*dt.*K1, param);
    K3 = pendulum_double_angle_aug_variational_vectorized(X_aug + 0.5.*dt.*K2, param);
    K4 = pendulum_double_angle_aug_variational_vectorized(X_aug + dt.*K3,      param);
    X_aug = X_aug + dt./6 .* (K1 + 2.*K2 + 2.*K3 + K4);
   
    % Extract current variation to sum up the deviations         
    delta_X_k = X_aug(5:8,:);
    Deviation_flat = Deviation_flat + sqrt(sum(delta_X_k.^2));

end
% Fetch data from GPU
Deviation_flat = gather(Deviation_flat);

% Reshape to form the corret image
switch render_mode
    case 'symetry'
        Deviation_HalfMat = reshape(Deviation_flat, n, n_half);
        Deviation_Mat = [flip(flip(Deviation_HalfMat(:,2:end),1),2),Deviation_HalfMat];
        % Measure elapsed time
        time_on_GPU_symetry = toc
        time_per_gridpoint_ms = time_on_GPU_symetry/(n*n_half) * 1000
        alpha1_span = linspace(-pi/1, alpha1_ub, n);
    case 'free_zoom'
        Deviation_Mat = reshape(Deviation_flat, n, n);
        % Measure elapsed time
        n
        time_on_GPU_free_zoom = toc
        time_per_gridpoint_ms = time_on_GPU_free_zoom/(n^2) * 1000
end


%%% Prepare the plots
fig = figure('Position', [200, 100, 1200, 1200], 'Color', 'k','Name', 'Double Pendulum Fractal'); % For the dynamcis
ax = axes(fig, 'Position', [0.08, 0.08, 0.88, 0.88]);

% Convert radians to degrees
alpha1_span = alpha1_span .* 180./pi;
alpha2_span = alpha2_span .* 180./pi;
[A1, A2] = meshgrid(alpha1_span, alpha2_span);

% Plot
surf(A1, A2, log(Deviation_Mat), 'EdgeAlpha', 0)

xlabel('$\phi_1 \; [^{\circ}]$', 'FontSize', 20, 'Color', 'w')
ylabel('$\phi_2 \; [^{\circ}]$', 'FontSize', 20, 'Color', 'w')

xlim([alpha1_span(1), alpha1_span(end)])
ylim([alpha2_span(1), alpha2_span(end)])

% Apply the magic
colormap(inferno)
colormap(magma)
colormap(viridis)
colormap(twilight_shifted)
colormap(hot)
colormap(flip(hot))
colormap((slanCM('seasons')))

% colormap(flip(slanCM('fusion')))
% colormap((slanCM('tokyo')))
% colormap(flip(slanCM('gist_heat')))
% colormap(slanCM('pride'))
% colormap((slanCM('copper2')))
% colormap((slanCM('emergency')))
% colormap(flip(slanCM('gist_ncar')))
% colormap(flip(slanCM('gnuplot')))

% ax.Position = [0, 0, 1, 1];
% ax.XTicks = []; 
% ax.YTicks = [];


ax.View = [0,90];


% hold on
% plot3([-pi,pi],[0,0],[10,10],'r-')
% plot3([0,0],[-pi,pi],[10,10],'r-')

grid off
ax.Color = 'k';
ax.XAxis.Color = 'w';
ax.YAxis.Color = 'w';

ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

set(findall(fig,'-property','Box'),'Box','off') % optional
set(findall(fig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
% print(hfig,fname,'-dpng','-painters')fig = figure('Position', [200, 100, 1200, 1200], 'Name', 'Double Pendulum Fractal'); % For the dynamcis





% % Plot as image
% fig = figure('Position', [1000, 100, 1200, 1200], 'Name', 'Double Pendulum Fractal'); % For the dynamcis
% ax = axes(fig, 'Position', [0.08, 0.08, 0.88, 0.88]);
% 
% % Flip and normalize
% Deviation_Mat_fliped_log = flip(log(Deviation_Mat));
% Deviation_Mat_normalized = (Deviation_Mat_fliped_log - min(Deviation_Mat_fliped_log,[],'all'))./(max(Deviation_Mat_fliped_log,[],'all') - min(Deviation_Mat_fliped_log,[],'all'));

% imagesc(Deviation_Mat_normalized)
% colormap((slanCM('seasons')))







% %%% SVD Plots
% fig = figure('Position', [600, 100, 1200, 1200], 'Color', 'k','Name', 'Double Pendulum Fractal'); % For the dynamcis
% ay = axes(fig, 'Position', [0.08, 0.08, 0.88, 0.88]);
% 
% [U,S,V] = svd(log(Deviation_Mat));
% 
% suf = surf(ay, A1, A2, U(:,1:1)*S(1:1,1:1)*V(:,1:1)', 'EdgeAlpha', 0);
% colormap(ay, (slanCM('seasons')))
% ay.View = [0,90];
% 
% pause(1)
% for i = 1:50
% 
%     lowREs = U(:,1:i)*S(1:i,1:i)*V(:,1:i)';
% 
%     suf.ZData = lowREs;
%     pause(1/10)
% end