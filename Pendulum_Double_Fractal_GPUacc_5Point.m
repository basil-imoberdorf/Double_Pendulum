clear
close all
reset(gpuDevice)

% SUPER IMPORTANT FOR THE FINAL IMAGE IS HOW MANY EVALUATIONS OF THE CHAOS METRIC ALONG THE TRAJECTORY ARE SUMMED UP FOR THE FINAL METRIC. USE ONLY THE FINAL STATE TO HAVE MORE SHARP FEATURES!!! 
% IF YOU SUM UP DIFFERENCES SAMPLED FROM THE WHOLE TRAJECTORY YOU WILL SMEAR OUT SOME EFFECTS SINCE IN THE BEGINNING ALOT OF POSITONS ARE MORE STABLE THAN LATER ON

% Simulation set-up
T = 9.9; % [s] 
dt = 0.01; % [s] RK4 step size
Nt = ceil(T/dt);
eps = 1e-5; % perturbation strength

% Pendulum design
m1 = 1;  % [kg]
m2 = 1;  % [kg]
L1 = 1;  % [m]
L2 = 1;  % [m]
g = 10;  % [m/s^2]
param = [m1;L1;m2;L2;g];

%%% Define grid resolution
n = 3201; % must be an ODD number when used with render mode 'symmetry' so there will be a center line.

hyperparam_str = ['_T_',num2str(T),'_dt_',num2str(dt),'_eps_',num2str(eps, '%.0e'),'_m1m2_',num2str(m1),'_L1L2_',num2str(L1),'_n_',num2str(n),'_GPU_RK4']

% Chose render mode
render_mode = 'symmetry'; % 'symmetry' or 'free_zoom'
switch render_mode
    case 'symmetry'
        n_half = ceil(n/2); % we will simulate only on half of the first angle range since the other half is the same
        alpha_lb = 0 * pi/180; 
        alpha_ub = 180 * pi/180;
        beta_lb = -180 * pi/180;
        beta_ub =  180 * pi/180; 
        alpha_span = linspace(alpha_lb, alpha_ub, n_half);
        beta_span  = linspace(beta_lb, beta_ub, n);
        Deviation_flat = zeros(1, n*n_half);   % Pre-allocate output vector
    case 'free_zoom'
        alpha_lb = 20  * pi/180;
        alpha_ub = 180  * pi/180;
        beta_lb = -180 * pi/180;
        beta_ub =  50  * pi/180;
        alpha_span = linspace(alpha_lb, alpha_ub, n);
        beta_span = linspace(beta_lb, beta_ub, n);
        Deviation_flat = zeros(1, n^2);   % Pre-allocate output vector
end

% Create grid of initial conditions
[A, B] = meshgrid(alpha_span, beta_span);
% Flatten and transpose the grid before adding the perturbation
alpha_pos = A(:).' + eps;
alpha_neg = A(:).' - eps;
beta_pos  = B(:).' + eps;
beta_neg  = B(:).' - eps;
alpha_nominal = A(:).';
beta_nominal  = B(:).';

Np = numel(alpha_pos); % Number of double pendula simulated in parallel

% Create matrices of intial condition. Each column represents the state of a double pendulum
zeros_padding = zeros(2,Np);
X0_beta_pos = [alpha_nominal; beta_pos; zeros_padding];
X0_beta_neg = [alpha_nominal; beta_neg; zeros_padding];
X0_alpha_neg = [alpha_neg; beta_nominal; zeros_padding];
X0_alpha_pos = [alpha_pos; beta_nominal; zeros_padding];
X0_nominal = [alpha_nominal; beta_nominal; zeros_padding];

%%% Move the necessary variables to the GPU
tic
% param = gpuArray(param); % oddly the whole thing runs almost twice as fast if we DO NOT pass the parameters to the GPU
X_beta_pos  = gpuArray(X0_beta_pos);
X_beta_neg  = gpuArray(X0_beta_neg);
X_alpha_neg = gpuArray(X0_alpha_neg);
X_alpha_pos = gpuArray(X0_alpha_pos);
X_nominal = gpuArray(X0_nominal);
    
% Perform Runge-Kutta 4 integration on the GPU
for k = 1:Nt
    % Integrate positive perturbation on the first links angle
    K1 = pendulum_double_angle_minTrig_rhs(X_beta_pos, param);
    K2 = pendulum_double_angle_minTrig_rhs(X_beta_pos + 0.5.*dt.*K1, param);
    K3 = pendulum_double_angle_minTrig_rhs(X_beta_pos + 0.5.*dt.*K2, param);
    K4 = pendulum_double_angle_minTrig_rhs(X_beta_pos + dt.*K3, param);
    X_beta_pos = X_beta_pos + dt./6 .* (K1 + 2.*K2 + 2.*K3 + K4);
    % Integrate negative perturbation on the first links angle
    K1 = pendulum_double_angle_minTrig_rhs(X_beta_neg, param);
    K2 = pendulum_double_angle_minTrig_rhs(X_beta_neg + 0.5.*dt.*K1, param);
    K3 = pendulum_double_angle_minTrig_rhs(X_beta_neg + 0.5.*dt.*K2, param);
    K4 = pendulum_double_angle_minTrig_rhs(X_beta_neg + dt.*K3, param);
    X_beta_neg = X_beta_neg + dt./6 .* (K1 + 2.*K2 + 2.*K3 + K4);
    % Integrate negative perturbation on the second links angle
    K1 = pendulum_double_angle_minTrig_rhs(X_alpha_neg, param);
    K2 = pendulum_double_angle_minTrig_rhs(X_alpha_neg + 0.5.*dt.*K1, param);
    K3 = pendulum_double_angle_minTrig_rhs(X_alpha_neg + 0.5.*dt.*K2, param);
    K4 = pendulum_double_angle_minTrig_rhs(X_alpha_neg + dt.*K3, param);
    X_alpha_neg = X_alpha_neg + dt./6 .* (K1 + 2.*K2 + 2.*K3 + K4);
    % Integrate positive perturbation on the second links angle
    K1 = pendulum_double_angle_minTrig_rhs(X_alpha_pos, param);
    K2 = pendulum_double_angle_minTrig_rhs(X_alpha_pos + 0.5.*dt.*K1, param);
    K3 = pendulum_double_angle_minTrig_rhs(X_alpha_pos + 0.5.*dt.*K2, param);
    K4 = pendulum_double_angle_minTrig_rhs(X_alpha_pos + dt.*K3, param);
    X_alpha_pos = X_alpha_pos + dt./6 .* (K1 + 2.*K2 + 2.*K3 + K4);
    % Integrate nominal trajecotry 
    K1 = pendulum_double_angle_minTrig_rhs(X_nominal, param);
    K2 = pendulum_double_angle_minTrig_rhs(X_nominal + 0.5.*dt.*K1, param);
    K3 = pendulum_double_angle_minTrig_rhs(X_nominal + 0.5.*dt.*K2, param);
    K4 = pendulum_double_angle_minTrig_rhs(X_nominal + dt.*K3, param);
    X_nominal = X_nominal + dt./6 .* (K1 + 2.*K2 + 2.*K3 + K4);
end

% Fetch data from GPU
X_beta_pos  = gather(X_beta_pos);
X_beta_neg  = gather(X_beta_neg);
X_alpha_neg = gather(X_alpha_neg);
X_alpha_pos = gather(X_alpha_pos);
X_nominal   = gather(X_nominal);

% Save the raw GPU output
save('X_beta_pos__T9.9_3201n.mat', 'X_beta_pos');
save('X_beta_neg__T9.9_3201n.mat', 'X_beta_neg');
save('X_alpha_neg_T9.9_3201n.mat', 'X_alpha_neg');
save('X_alpha_pos_T9.9_3201n.mat', 'X_alpha_pos');
save('X_nominal___T9.9_3201n.mat', 'X_nominal');

%%% Compute the chaos metrics
% First only using the position of the two angles
delta_beta_pos__POSITION = sqrt(sum((X_nominal(1:2,:) - X_beta_pos(1:2,:) ).^2, 1));
delta_beta_neg__POSITION = sqrt(sum((X_nominal(1:2,:) - X_beta_neg(1:2,:) ).^2, 1));
delta_alpha_pos_POSITION = sqrt(sum((X_nominal(1:2,:) - X_alpha_neg(1:2,:)).^2, 1));
delta_alpha_neg_POSITION = sqrt(sum((X_nominal(1:2,:) - X_alpha_pos(1:2,:)).^2, 1));
Position_Deviation_flat = mean([delta_beta_pos__POSITION; delta_beta_neg__POSITION; delta_alpha_pos_POSITION; delta_alpha_neg_POSITION]);
% Then only using the angular velocity of the two angles
delta_beta_pos__VELOCITY = sqrt(sum((X_nominal(3:4,:) - X_beta_pos(3:4,:) ).^2, 1));
delta_beta_neg__VELOCITY = sqrt(sum((X_nominal(3:4,:) - X_beta_neg(3:4,:) ).^2, 1));
delta_alpha_pos_VELOCITY = sqrt(sum((X_nominal(3:4,:) - X_alpha_neg(3:4,:)).^2, 1));
delta_alpha_neg_VELOCITY = sqrt(sum((X_nominal(3:4,:) - X_alpha_pos(3:4,:)).^2, 1));
Velocity_Deviation_flat = mean([delta_beta_pos__VELOCITY; delta_beta_neg__VELOCITY; delta_alpha_pos_VELOCITY; delta_alpha_neg_VELOCITY]);
% Finally using the full state
delta_beta_pos__FULLSTATE = sqrt(sum((X_nominal - X_beta_pos ).^2, 1));     %delta_beta_pos_FULLSTATE =  sum(abs(X_nominal - X_beta_pos ), 1);
delta_beta_neg__FULLSTATE = sqrt(sum((X_nominal - X_beta_neg ).^2, 1));     %delta_beta_neg_FULLSTATE =  sum(abs(X_nominal - X_beta_neg ), 1);
delta_alpha_pos_FULLSTATE = sqrt(sum((X_nominal - X_alpha_neg).^2, 1));     %delta_alpha_pos_FULLSTATE = sum(abs(X_nominal - X_alpha_neg), 1);
delta_alpha_neg_FULLSTATE = sqrt(sum((X_nominal - X_alpha_pos).^2, 1));     %delta_alpha_neg_FULLSTATE = sum(abs(X_nominal - X_alpha_pos), 1);
FullState_Deviation_flat = mean([delta_beta_pos__FULLSTATE; delta_beta_neg__FULLSTATE; delta_alpha_pos_FULLSTATE; delta_alpha_neg_FULLSTATE]);

% Compute Lyapunov Exponent
Lyapunov_beta_pos_FULLSTATE  = 1./T .* log( delta_beta_pos__FULLSTATE ./ sqrt(sum((X0_nominal - X0_beta_pos ).^2, 1)) );
Lyapunov_beta_neg_FULLSTATE  = 1./T .* log( delta_beta_neg__FULLSTATE ./ sqrt(sum((X0_nominal - X0_beta_neg ).^2, 1)) );
Lyapunov_alpha_pos_FULLSTATE = 1./T .* log( delta_alpha_pos_FULLSTATE ./ sqrt(sum((X0_nominal - X0_alpha_neg).^2, 1)) );
Lyapunov_alpha_neg_FULLSTATE = 1./T .* log( delta_alpha_neg_FULLSTATE ./ sqrt(sum((X0_nominal - X0_alpha_pos).^2, 1)) );
FullState_Lyapunov_flat_AVG = mean([Lyapunov_beta_pos_FULLSTATE; Lyapunov_beta_neg_FULLSTATE; Lyapunov_alpha_pos_FULLSTATE; Lyapunov_alpha_neg_FULLSTATE], 1);
FullState_Lyapunov_flat_MAX = max([Lyapunov_beta_pos_FULLSTATE; Lyapunov_beta_neg_FULLSTATE; Lyapunov_alpha_pos_FULLSTATE; Lyapunov_alpha_neg_FULLSTATE], [],1);


% Reshape to form the correct image
switch render_mode
    case 'symmetry'
        Deviation_HalfMat_POSITION  = reshape(Position_Deviation_flat,  n, n_half);
        Deviation_HalfMat_VELOCITY  = reshape(Velocity_Deviation_flat,  n, n_half);
        Deviation_HalfMat_FULLSTATE = reshape(FullState_Deviation_flat, n, n_half);
        Deviation_Mat_POSITION  = [flip(flip(Deviation_HalfMat_POSITION(:,2:end),1),2), Deviation_HalfMat_POSITION ];
        Deviation_Mat_VELOCITY  = [flip(flip(Deviation_HalfMat_VELOCITY(:,2:end),1),2), Deviation_HalfMat_VELOCITY ];
        Deviation_Mat_FULLSTATE = [flip(flip(Deviation_HalfMat_FULLSTATE(:,2:end),1),2),Deviation_HalfMat_FULLSTATE];
        Lyapunov_AVG_HalfMat_FULLSTATE = reshape(FullState_Lyapunov_flat_AVG, n, n_half);
        Lyapunov_MAX_HalfMat_FULLSTATE = reshape(FullState_Lyapunov_flat_MAX, n, n_half);
        Lyapunov_AVG_Mat_FULLSTATE = [flip(flip(Lyapunov_AVG_HalfMat_FULLSTATE(:,2:end),1),2),Lyapunov_AVG_HalfMat_FULLSTATE];
        Lyapunov_MAX_Mat_FULLSTATE = [flip(flip(Lyapunov_MAX_HalfMat_FULLSTATE(:,2:end),1),2),Lyapunov_MAX_HalfMat_FULLSTATE];
        % Measure elapsed time
        GPU_total_time_symetry = toc
        GPU_Gridpoint_time_ms = GPU_total_time_symetry/(n*n_half) * 1000
        alpha_span = linspace(-pi/1, alpha_ub, n);
    case 'free_zoom'
        Deviation_Mat_POSITION = reshape(Position_Deviation_flat, n, n);
        Deviation_Mat_VELOCITY = reshape(Velocity_Deviation_flat, n, n);
        Deviation_Mat_FULLSTATE = reshape(FullState_Deviation_flat, n, n);
        Lyapunov_AVG_Mat_FULLSTATE = reshape(FullState_Lyapunov_flat_AVG, n, n);
        Lyapunov_MAX_Mat_FULLSTATE = reshape(FullState_Lyapunov_flat_MAX, n, n);
        % Measure elapsed time
        n
        GPU_total_time_free_zoom = toc
        GPU_Gridpoint_time_ms = GPU_total_time_free_zoom/(n^2) * 1000
end

% Take the logarithm so the colormap spans more evenly on the highly nonlinear data
Deviation_Mat_FULLSTATE    = log(Deviation_Mat_FULLSTATE   );
Deviation_Mat_POSITION     = log(Deviation_Mat_POSITION    );
Deviation_Mat_VELOCITY     = log(Deviation_Mat_VELOCITY    );
Lyapunov_MAX_Mat_FULLSTATE = ((Lyapunov_MAX_Mat_FULLSTATE));
Lyapunov_AVG_Mat_FULLSTATE = ((Lyapunov_AVG_Mat_FULLSTATE));

% Min-Max scaling to facilitate
Deviation_Mat_FULLSTATE    = (Deviation_Mat_FULLSTATE    - min(Deviation_Mat_FULLSTATE(:)))    ./ (max(Deviation_Mat_FULLSTATE(:))    - min(Deviation_Mat_FULLSTATE(:))); 
Deviation_Mat_POSITION     = (Deviation_Mat_POSITION     - min(Deviation_Mat_POSITION(:)))     ./ (max(Deviation_Mat_POSITION(:))     - min(Deviation_Mat_POSITION(:)));
Deviation_Mat_VELOCITY     = (Deviation_Mat_VELOCITY     - min(Deviation_Mat_VELOCITY(:)))     ./ (max(Deviation_Mat_VELOCITY(:))     - min(Deviation_Mat_VELOCITY(:)));
Lyapunov_MAX_Mat_FULLSTATE = (Lyapunov_MAX_Mat_FULLSTATE - min(Lyapunov_MAX_Mat_FULLSTATE(:))) ./ (max(Lyapunov_MAX_Mat_FULLSTATE(:)) - min(Lyapunov_MAX_Mat_FULLSTATE(:)));    
Lyapunov_AVG_Mat_FULLSTATE = (Lyapunov_AVG_Mat_FULLSTATE - min(Lyapunov_AVG_Mat_FULLSTATE(:))) ./ (max(Lyapunov_AVG_Mat_FULLSTATE(:)) - min(Lyapunov_AVG_Mat_FULLSTATE(:)));    

%%% Plot results
[srf1, ax1, fig1] = plot_fractal(Deviation_Mat_FULLSTATE, alpha_span, beta_span, [20,   50, 1200, 1200], 'FULLSTATE');
[srf2, ax2, fig2] = plot_fractal(Deviation_Mat_POSITION,  alpha_span, beta_span, [300, 130, 1200, 1200], 'POSITION'); 
[srf3, ax3, fig3] = plot_fractal(Deviation_Mat_VELOCITY,  alpha_span, beta_span, [600, 200, 1200, 1200], 'VELOCITY'); 
[srf4, ax4, fig4] = plot_fractal(Lyapunov_MAX_Mat_FULLSTATE,  alpha_span, beta_span, [1050, 50,  1200, 1200], 'LYAPUNOV - MAX');
[srf5, ax5, fig5] = plot_fractal(Lyapunov_AVG_Mat_FULLSTATE,  alpha_span, beta_span, [1350, 200, 1200, 1200], 'LYAPUNOV - AVG');  
all_axes = [ax1, ax2, ax3, ax4, ax5];


% colormap(flip(slanCM('fusion')))
% colormap((slanCM('tokyo')))
% colormap(slanCM('pride'))
% colormap((slanCM('copper2')))
% colormap((slanCM('emergency')))
% colormap(flip(slanCM('gist_ncar')))
% colormap(flip(slanCM('gnuplot')))



% arrayfun(@(ax) colormap(ax, slanCM('tokyo')), all_axes);
% arrayfun(@(ax) colormap(ax, flip(hot)), all_axes);
% arrayfun(@(ax) colormap(ax, slanCM('seasons')), all_axes);
% arrayfun(@(ax) colormap(ax, slanCM('seasons')), all_axes);
% arrayfun(@(ax) colormap(ax, slanCM('emergency')), all_axes);
% arrayfun(@(ax) colormap(ax, flip(slanCM('chroma'))), all_axes);
% arrayfun(@(ax) colormap(ax, flip(slanCM('gnuplot2'))), all_axes);


% cmap = slanCM('emergency');   % N x 3 RGB, assumed in [0,1]
% cmap = flip(slanCM('heat'));   % N x 3 RGB, assumed in [0,1]
% cmap = flip(slanCM('gnuplot2'));   % N x 3 RGB, assumed in [0,1]
% 
% cmap = flip(slanCM('afmhot'));
% cmap = flip(bone);
% cmap = flip(hot);
% cmap = flip(slanCM('amber'));
cmap = slanCM('seasons');
% cmap = [hot;flip(hot)];
% cmap = [bone;flip(bone)];
% cmap = [hot;flip(bone)];
% cmap = [hot;flip(slanCM('afmhot'))];
% cmap = [hot;flip(slanCM('ember'))];
% cmap = [bone;flip(slanCM('ember'))];
% cmap = [slanCM('heat');flip(bone)];
% cmap = [slanCM('sunburst');flip(bone)];
% cmap = [slanCM('torch');flip(bone)];
% cmap = [slanCM('afmhot');(slanCM('deep'))];
% cmap = [slanCM('afmhot');flip(magma)];
% cmap = [slanCM('afmhot');flip(bone)];
% cmap = [slanCM('ember');flip(bone)];
% cmap = [slanCM('afmhot');flip(slanCM('voltage'))];
% cmap = [hot;flip(slanCM('voltage'))];
% cmap = [slanCM('ember');flip(bone)];
% cmap = [hot;flip(bone)]; % FINAL [hot;flip(bone)]

Sat = 1.4; % FINAL 6
HSV = rgb2hsv(cmap);
HSV(:,2) = min(HSV(:,2) * Sat, 1);  % saturation factor
HSV(:,3) = min(HSV(:,3) * 1, 1);  % saturation factor
cmap_sat = hsv2rgb(HSV);
arrayfun(@(ax) colormap(ax, cmap_sat), all_axes);

% arrayfun(@(ax) caxis(ax, [0,1]), all_axes);

% arrayfun(@(ax) colormap(ax, flip(hot)), all_axes);
% arrayfun(@(ax) caxis(ax, [0,0.45]), all_axes);

% arrayfun(@(ax) caxis(ax, [0,1]), all_axes);
% arrayfun(@(ax) caxis(ax, [0.05,0.7]), all_axes);
% arrayfun(@(ax) caxis(ax, [0.12,0.62]), all_axes);
% arrayfun(@(ax) caxis(ax, [0.15,0.65]), all_axes);

 % arrayfun(@(ax) caxis(ax, [0.02,0.72]), all_axes); % FINAL [0.02,0.72]

% % Plot as image
% fig = figure('Position', [1000, 100, 1200, 1200], 'Name', 'Double Pendulum Fractal'); % For the dynamcis
% ax = axes(fig, 'Position', [0.08, 0.08, 0.88, 0.88]);
% 
% % Flip and normalize
% Deviation_Mat_fliped_log = flip(log(Deviation_Mat));
% Deviation_Mat_normalized = (Deviation_Mat_fliped_log - min(Deviation_Mat_fliped_log,[],'all'))./(max(Deviation_Mat_fliped_log,[],'all') - min(Deviation_Mat_fliped_log,[],'all'));

% imagesc(Deviation_Mat_normalized)
% colormap((slanCM('seasons')))


ticks  = [-pi, -pi/2, 0, pi/2, pi] * 180/pi;
labels = { ...
    '$-\pi$', ...
    '$-\frac{\pi}{2}$', ...
    '$0$', ...
    '$\frac{\pi}{2}$', ...
    '$\pi$' ...
};
arrayfun(@(ax) set(ax, ...
    'XTick', ticks, ...
    'XTickLabel', labels, ...
    'YTick', ticks, ...
    'YTickLabel', labels, ...
    'TickLabelInterpreter','latex', ...
    'FontSize', 20, ...
    'FontWeight','bold', ...
    'LineWidth',2, ...
    'Box','on', ...
    'XMinorTick','on', ...
    'YMinorTick','on', ...
    'TickDir','out', ...
    'Projection','orthographic'), ...
    all_axes);
 
arrayfun(@(ax) set(ax.Title, ...
    'FontSize', 30, ...
    'FontWeight', 'bold', ...
    'Position', [0 185 0.5]), ...
    all_axes);

arrayfun(@(ax) set(ax, ...
    'Position', [0.09, 0.09, 0.86, 0.86]), ...
    all_axes);

arrayfun(@(ax) set(ax.XLabel, ...
    'Position', [0, -195, 0]), ...
    all_axes);
arrayfun(@(ax) set(ax.YLabel, ...
    'Position', [-195, 0, 0]), ...
    all_axes);



% 
% arrayfun(@(ax) set(ax.Title, ...
%     'Units','normalized', ...
%     'Position', ax.Title.Position - [0 0.001 0]), ...
%     all_axes);
% 
% 
% [90 185 0]
% 
% set(ax1.Title, 'Position', [0.5 1.1 0])



arrayfun(@(ax) set(ax.XLabel, ...
    'FontSize', 35), ...
    all_axes);
arrayfun(@(ax) set(ax.YLabel, ...
    'FontSize', 35), ...
    all_axes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [srf, ax, fig] = plot_fractal(Deviation_Mat, alpha_span, beta_span, fig_pos, title_str)
%%% Prepare the plots
fig = figure('Position', fig_pos, 'Color', 'k','Name', ['Double Pendulum Fractal_',title_str(2:7)]); % For the dynamcis
ax = axes(fig, 'Position', [0.07, 0.06, 0.88, 0.88]);

% Convert radians to degrees for easier interpretation of angles
alpha_span_deg = alpha_span .* 180./pi;
beta_span_deg = beta_span .* 180./pi;
[A_deg, B_deg] = meshgrid(alpha_span_deg, beta_span_deg);

% Plot
srf = surf(A_deg, B_deg, Deviation_Mat, 'EdgeAlpha', 0);

xlabel('$\alpha_0$', 'FontSize', 22, 'Color', 'w')  % [^{\circ}]
ylabel('$\beta_0 $', 'FontSize', 22, 'Color', 'w')  % [^{\circ}]

xlim([alpha_span_deg(1), alpha_span_deg(end)])
ylim([beta_span_deg(1), beta_span_deg(end)])

title(title_str, 'FontSize', 20, 'Color', 'w')

%%% Apply the magic
colormap((slanCM('seasons')))

grid off
ax.Color = 'k';
ax.XAxis.Color = 'w';
ax.YAxis.Color = 'w';
ax.ZAxis.Color = 'w';

ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.ZAxis.FontSize = 14;

ax.View = [0,90];
% ax.Position = [0, 0, 1, 1];
% ax.XTicks = [];
% ax.YTicks = [];


set(findall(fig,'-property','Box'),'Box','off') % optional
set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
% print(hfig,fname,'-dpng','-painters')

% clim([min(log(Deviation_Mat),[],"all"), 0.01* max(log(Deviation_Mat),[],"all")])


end







% Deviation_Mat = Lyapunov_MAX_Mat_FULLSTATE;
% fig = figure('Position', [200, 100, 1200, 1200], 'Color', 'k','Name', 'Double Pendulum Fractal'); % For the dynamcis
% ax = axes(fig, 'Position', [0.08, 0.08, 0.88, 0.88]);
% 
% % Convert radians to degrees for easier interpretation of angles
% alpha_span_deg = alpha_span .* 180./pi;
% beta_span_deg = beta_span .* 180./pi;
% [A_deg, B_deg] = meshgrid(alpha_span_deg, beta_span_deg);
% 
% % Plot
% surf(A_deg, B_deg, Deviation_Mat, 'EdgeAlpha', 0)
% 
% xlabel('$\alpha_0 \; [^{\circ}]$', 'FontSize', 20, 'Color', 'w')
% ylabel('$\beta_0 \;  [^{\circ}]$', 'FontSize', 20, 'Color', 'w')
% 
% xlim([alpha_span_deg(1), alpha_span_deg(end)])
% ylim([beta_span_deg(1), beta_span_deg(end)])
% 
% %%% Apply the magic
% % colormap(inferno)
% % colormap(magma)
% % colormap(viridis)
% % colormap(twilight_shifted)
% % colormap(hot)
% colormap(flip(hot))
% 
% % colormap((slanCM('seasons')))
% % colormap(flip(slanCM('gist_heat')))
% 
% % colormap(flip(slanCM('fusion')))
% % colormap((slanCM('tokyo')))
% % colormap(slanCM('pride'))
% % colormap((slanCM('copper2')))
% % colormap((slanCM('emergency')))
% % colormap(flip(slanCM('gist_ncar')))
% % colormap(flip(slanCM('gnuplot')))
% 
% % hold on
% % plot3([-pi,pi],[0,0],[10,10],'r-')
% % plot3([0,0],[-pi,pi],[10,10],'r-')
% 
% grid off
% ax.Color = 'k';
% ax.XAxis.Color = 'w';
% ax.YAxis.Color = 'w';
% 
% ax.XAxis.FontSize = 14;
% ax.YAxis.FontSize = 14;
% 
% ax.View = [0,90];
% % ax.Position = [0, 0, 1, 1];
% % ax.XTicks = [];
% % ax.YTicks = [];
% 
% set(findall(fig,'-property','Box'),'Box','off') % optional
% set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
% set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% % set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% % pos = get(fig,'Position');
% % set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% %print(hfig,fname,'-dpdf','-painters','-fillpage')
% % print(hfig,fname,'-dpng','-painters')



%%%
fig_final = figure( ...
    'Visible','off', ...
    'Units','pixels', ...
    'Position',[0 0 4200 4200], ...
    'PaperUnits','inches', ...
    'PaperPosition',[0 0 42 42], ...
    'PaperPositionMode','manual', ...
    'Color', 'k', ...
    'Name', ['Double Pendulum Fractal_',hyperparam_str(2:7)]);
ax = axes(fig_final, 'Position', [0.09, 0.09, 0.86, 0.86]);


% Convert radians to degrees for easier interpretation of angles
alpha_span_deg = alpha_span .* 180./pi;
beta_span_deg = beta_span .* 180./pi;
[A_deg, B_deg] = meshgrid(alpha_span_deg, beta_span_deg);

% Plot
srf = surf(ax, A_deg, B_deg, Deviation_Mat_FULLSTATE, 'EdgeAlpha', 0);

xlabel(ax, '$\alpha_0$', 'FontSize', 35, 'Color', 'w', 'Position', [0, -195, 0])  % [^{\circ}]
ylabel(ax, '$\beta_0 $', 'FontSize', 35, 'Color', 'w', 'Position', [-195, 0, 0])  % [^{\circ}]

xlim(ax, [alpha_span_deg(1), alpha_span_deg(end)])
ylim(ax, [beta_span_deg(1), beta_span_deg(end)])

title(ax, 'LYAPUNOV STABILITY', 'FontSize', 20, 'Color', 'w', 'Position', [0 185 0.5])

%%% Apply the magic
cmap = [hot;flip(bone)]; % FINAL [hot;flip(bone)]

Sat = 6; % FINAL 6
HSV = rgb2hsv(cmap);
HSV(:,2) = min(HSV(:,2) * Sat, 1);  % saturation factor
HSV(:,3) = min(HSV(:,3) * 1, 1);  % saturation factor
cmap_sat = hsv2rgb(HSV);

colormap(ax, cmap_sat);
caxis(ax, [0.02,0.72]); % FINAL [0.02,0.72]

grid off
ax.Color = 'k';
ax.XAxis.Color = 'w';
ax.YAxis.Color = 'w';
ax.ZAxis.Color = 'w';

ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.ZAxis.FontSize = 14;

ax.DataAspectRatio = [1,1,1];
ax.View = [0,90];


ticks  = [-pi, -pi/2, 0, pi/2, pi] * 180/pi;
labels = { ...
    '$-\pi$', ...
    '$-\frac{\pi}{2}$', ...
    '$0$', ...
    '$\frac{\pi}{2}$', ...
    '$\pi$'};

set(ax, ...
    'XTick', ticks, ...
    'XTickLabel', labels, ...
    'YTick', ticks, ...
    'YTickLabel', labels, ...
    'TickLabelInterpreter','latex', ...
    'FontSize', 20, ...
    'FontWeight','bold', ...
    'LineWidth',2, ...
    'Box','on', ...
    'XMinorTick','on', ...
    'YMinorTick','on', ...
    'TickDir','out', ...
    'Projection','orthographic')





% set(findall(fig_final,'-property','Box'),'Box','off') % optional
set(findall(fig_final,'-property','Interpreter'),'Interpreter','latex')
set(findall(fig_final,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

pause(1)
exportgraphics(fig_final,'double_pendulum_ggfractal_axes.png');
close(fig_final)


