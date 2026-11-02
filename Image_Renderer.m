clear
close all

% Dimension of target image
n = 3201; % must be an ODD number when used with render mode 'symmetry' so there will be a center line.
n_half = ceil(n/2); 
T = 8.9;

% Load the matrices
% X_beta_pos  = load('X_beta_pos__T8.9_4001n.mat');
% X_beta_neg  = load('X_beta_neg__T8.9_4001n.mat');
% X_alpha_neg = load('X_alpha_neg_T8.9_4001n.mat');
% X_alpha_pos = load('X_alpha_pos_T8.9_4001n.mat');
% X_nominal   = load('X_nominal___T8.9_4001n.mat');
% 
% X_beta_pos  = load('X_beta_pos__T8.9_3201n.mat');
% X_beta_neg  = load('X_beta_neg__T8.9_3201n.mat');
% X_alpha_neg = load('X_alpha_neg_T8.9_3201n.mat');
% X_alpha_pos = load('X_alpha_pos_T8.9_3201n.mat');
% X_nominal   = load('X_nominal___T8.9_3201n.mat');

X_beta_pos  = load('X_beta_pos__T9.9_3201n.mat');
X_beta_neg  = load('X_beta_neg__T9.9_3201n.mat');
X_alpha_neg = load('X_alpha_neg_T9.9_3201n.mat');
X_alpha_pos = load('X_alpha_pos_T9.9_3201n.mat');
X_nominal   = load('X_nominal___T9.9_3201n.mat');

X_beta_pos  = X_beta_pos.X_beta_pos;
X_beta_neg  = X_beta_neg.X_beta_neg;
X_alpha_neg = X_alpha_neg.X_alpha_neg;
X_alpha_pos = X_alpha_pos.X_alpha_pos;
X_nominal   = X_nominal.X_nominal;

% X_beta_pos  = readmatrix('X_sym_alpha_neg_T_10_dt_0.001_eps_1e-05_m1m2_1_L1L2_1_n_3001_GPU_RK4.csv');
% X_beta_neg  = readmatrix('X_sym_alpha_pos_T_10_dt_0.001_eps_1e-05_m1m2_1_L1L2_1_n_3001_GPU_RK4.csv');
% X_alpha_neg = readmatrix('X_sym_beta_neg__T_10_dt_0.001_eps_1e-05_m1m2_1_L1L2_1_n_3001_GPU_RK4.csv');
% X_alpha_pos = readmatrix('X_sym_beta_pos__T_10_dt_0.001_eps_1e-05_m1m2_1_L1L2_1_n_3001_GPU_RK4.csv');
% X_nominal   = readmatrix('X_sym_nominal___T_10_dt_0.001_eps_1e-05_m1m2_1_L1L2_1_n_3001_GPU_RK4.csv');


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

% Reshape to form the correct image
Deviation_HalfMat_POSITION  = reshape(Position_Deviation_flat,  n, n_half);
Deviation_HalfMat_VELOCITY  = reshape(Velocity_Deviation_flat,  n, n_half);
Deviation_HalfMat_FULLSTATE = reshape(FullState_Deviation_flat, n, n_half);
Deviation_Mat_POSITION  = [flip(flip(Deviation_HalfMat_POSITION(:,2:end),1),2), Deviation_HalfMat_POSITION ];
Deviation_Mat_VELOCITY  = [flip(flip(Deviation_HalfMat_VELOCITY(:,2:end),1),2), Deviation_HalfMat_VELOCITY ];
Deviation_Mat_FULLSTATE = [flip(flip(Deviation_HalfMat_FULLSTATE(:,2:end),1),2),Deviation_HalfMat_FULLSTATE];

% Take the logarithm so the colormap spans more evenly on the highly nonlinear data
Deviation_Mat_FULLSTATE    = log(Deviation_Mat_FULLSTATE   );
Deviation_Mat_POSITION     = log(Deviation_Mat_POSITION    );
Deviation_Mat_VELOCITY     = log(Deviation_Mat_VELOCITY    );

% Min-Max scaling to facilitate
Deviation_Mat_FULLSTATE    = (Deviation_Mat_FULLSTATE    - min(Deviation_Mat_FULLSTATE(:)))    ./ (max(Deviation_Mat_FULLSTATE(:))    - min(Deviation_Mat_FULLSTATE(:))); 
Deviation_Mat_POSITION     = (Deviation_Mat_POSITION     - min(Deviation_Mat_POSITION(:)))     ./ (max(Deviation_Mat_POSITION(:))     - min(Deviation_Mat_POSITION(:)));
Deviation_Mat_VELOCITY     = (Deviation_Mat_VELOCITY     - min(Deviation_Mat_VELOCITY(:)))     ./ (max(Deviation_Mat_VELOCITY(:))     - min(Deviation_Mat_VELOCITY(:)));


Raw_Image = flip(Deviation_Mat_FULLSTATE,1);

% 
% %%% MAP THE PIXEL VALUES TO RGB COLORS
% fig = figure('Position', [300, 130, 1200, 1200],'Color', [0.2,0.2,0.2],'Name', 'Double Pendulum Fractal'); % For the dynamcis
% ax = axes(fig);
% 
% 
% %%% Convert 2D to RGB image by applying a colormap: This expects normalized values
% cmap =  [hot;flip(bone)];
% Sat = 6; % Saturation
% HSV = rgb2hsv(cmap);
% HSV(:,2) = min(HSV(:,2) * Sat, 1);  % s
% cmap = hsv2rgb(HSV);
% cLim = [0.02,0.72];
% Raw_Image = (Raw_Image - cLim(1)) ./ (cLim(2) - cLim(1));
% Raw_Image = min(max(Raw_Image, 0), 1);
% 
% 
% nColors = size(cmap, 1);  % This are two 8-bit colormaps concatenated
% idx = round(Raw_Image .* (nColors-1)) + 1;  % assigns each pixel the closest color index. + 1 to start from 1
% 
% % Initialize RGB output
% [n, m] = size(Raw_Image);
% Image_RGB = zeros(n, m, 3);
% 
% % Assign RGB values
% for ch = 1:3
%     channel = cmap(:,ch);
%     Image_RGB(:,:,ch) = channel(idx);
% end
% 
% 
% %%% Padding
% padding = 60; % Even number pls! number rows/columns added to the matrix on the left and right and rounded corners
% padding_color = [1, 1, 1];
% 
% % % Padding with rounded corners
% Padded_Image = round_padding(Image_RGB, n, padding, padding_color);
% 
% % % Alternative Padding
% % Major_TickScale = 4;
% % Minor_TickScale = 2;
% % Padded_Image = tick_padding(Image_RGB, n, padding, padding_color, Major_TickScale, Minor_TickScale);
% 
% 
% % Show image
% % imshow(Padded_Image)
% % imwrite(Padded_Image, 'double_pendulum_gg_final_3201x3201_TickPad_stahl.png');
% % imwrite(Image_RGB, 'double_pendulum_gg_final_3201x3201_NoPad.png','-dpng','-r300');
% 
% dpi = 200;
% pixelsPerMeter = dpi / 0.0254;
% imwrite(Image_RGB, 'double_pendulum_gg_final_3201x3201_NoPad_200DPI.png', ...
%         'XResolution', pixelsPerMeter, ...
%         'YResolution', pixelsPerMeter, ...
%         'ResolutionUnit', 'meter');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert radians to degrees for easier interpretation of angles
alpha_span = linspace(-pi, pi, n);
beta_span = linspace(-pi, pi, n);
alpha_span_deg = alpha_span .* 180./pi;
beta_span_deg = beta_span .* 180./pi;
[A_deg, B_deg] = meshgrid(alpha_span_deg, beta_span_deg);



%%% Make figure with decorations
fig_final = figure( ...
    'Name', 'HurenFraktal',...
    'Visible','off', ...
    'Units','pixels', ...
    'Position',[0 0 3401 3401], ...
    'PaperUnits','inches', ...
    'PaperPosition',[0 0 13.2 13.2], ...
    'PaperPositionMode','manual', ...
    'Color', 'w');
ax = axes(fig_final, 'Position', [0.095, 0.095, 0.83, 0.83]);
% ax = axes(fig_final, 'Position', [0.08, 0.08, 0.84, 0.84]);
% ax = axes(fig_final, 'Position', [0.07, 0.07, 0.85, 0.85]);
% ax = axes(fig_final, 'Position', [0.06, 0.075, 0.89, 0.905]);

% Plot
% srf = surf(ax, A_deg, B_deg, Deviation_Mat_FULLSTATE, 'EdgeAlpha', 0);
srf = imagesc(ax, Raw_Image);

ticks_rad  = [-pi, -pi/2, 0, pi/2, pi];
labels = { ...
    '$-\pi$', ...
    '$-\frac{\pi}{2}$', ...
    '$0$', ...
    '$\frac{\pi}{2}$', ...
    '$\pi$'};

ticks = (ticks_rad + pi) .* n/(2*pi);

set(ax, ...
    'XTick', ticks, ...
    'XTickLabel', labels, ...
    'YTick', ticks, ...
    'YTickLabel',flip(labels), ...
    'TickLabelInterpreter','latex', ...
    'FontWeight','bold', ...
    'LineWidth',3, ...
    'Box','on', ...
    'XMinorTick','on', ...
    'YMinorTick','on', ...
    'TickDir','out', ...
    'Projection','orthographic')


title(ax, 'LYAPUNOV', 'FontSize', 72, 'Color', 'k', 'Position', ax.Title.Position + [0, 15, 0]) % 'Position', [0 185 0.5]
% title(ax, 'LYAPUNOV', 'FontSize', 100, 'Color', 'w', 'Position', ax.Title.Position + [0, 10, 0]) % 'Position', [0 185 0.5]


%%% Apply the magic
cmap = [hot;flip(bone)]; % FINAL [hot;flip(bone)]
cmap = slanCM('seasons');

Sat = 1.4; % FINAL 6
HSV = rgb2hsv(cmap);
HSV(:,2) = min(HSV(:,2) * Sat, 1);  % saturation factor
HSV(:,3) = min(HSV(:,3) * 1, 1);  % saturation factor
cmap_sat = hsv2rgb(HSV);

colormap(ax, cmap_sat);
% caxis(ax, [0.02,0.72]); % FINAL [0.02,0.72]

grid off
ax.Color = 'k';
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
ax.ZAxis.Color = 'k';

ax.XAxis.FontSize = 25;
ax.YAxis.FontSize = 25;
ax.ZAxis.FontSize = 25;


% ax.DataAspectRatio = [1,1,1];
ax.View = [0,90];

xlabel(ax, '$\alpha_0$', 'FontSize', 43, 'Color', 'k', 'Position', ax.XLabel.Position + [0, 30, 0])  % [^{\circ}] 'Position', [0, -195, 0]
ylabel(ax, '$\beta_0 $', 'FontSize', 43, 'Color', 'k', 'Position', ax.YLabel.Position + [-55, 0, 0])  % [^{\circ}] 'Position', [-195, 0, 0]


xlim(ax, [0-0.01, n+0.01])
ylim(ax, [0-0.01, n+0.01])


% set(findall(fig_final,'-property','Box'),'Box','off') % optional
set(findall(fig_final,'-property','Interpreter'),'Interpreter','latex')
set(findall(fig_final,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')


% exportgraphics(fig_final,'double_pendulum_ggfractal_axes6.png', 'BackgroundColor','none');

% print(fig_final,'double_pendulum_gg_TitleLableTicksR180_3201p_13inch.png','-dpng','-r180')
print(fig_final,'double_pendulum_gg_TitleLableTicksR200_3201p_13inch_Final.png','-dpng','-r200')
% print(fig_final,'double_pendulum_gg_TitleLableTicksR300_3201p_13inch.png','-dpng','-r300')



% % export_fig('double_pendulum_ggfractal_axes4.png', [fig_final]);
% 
% % % drawnow;                   % ensure the figure is fully rendered
% % F = getframe(fig_final);     % capture entire figure window
% % imwrite(F.cdata, 'double_pendulum_ggfractal_axes3.png');
% % size(F.cdata)
% 
% % set(fig_final,'Units','inches','Position',[3 3 21 21])
% % pos = get(fig_final,'Position');
% % set(fig_final,'PaperPositionMode','Auto','PaperUnits','inches','PaperSize',[pos(3), pos(4)])
% % % print(fig_final,'double_pendulum_ggfractal_axes3.png','-dpdf','-painters','-fillpage')
% % print(fig_final,'double_pendulum_ggfractal_axes3.png','-dpng','-image')
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function padded_image = tick_padding(Image_RGB, n, pad_size, color, Major_TickScale, Minor_TickScale)
% 
% % Presets
% Major_spike_breite = pad_size; 
% Minor_spike_breite = pad_size; 
% num_MajorTicks = 5;
% num_MinorTicks = 5*num_MajorTicks;
% 
% % Create padding matrix for one side. will be flipped for the others
% side_padd = zeros(n,pad_size*Major_TickScale);
% % Apply the color to the padding backbone
% side_padd(:, 1:pad_size, 1) = color(1);
% side_padd(:, 1:pad_size, 2) = color(2);
% side_padd(:, 1:pad_size, 3) = color(3);
% 
% % Get the indices where the tick marks will be drawn. Give indices closest to the fracs
% Majorfracs = (0:(num_MajorTicks-1)) ./ (num_MajorTicks-1);
% MajorTick_idx = round((n-1)*Majorfracs + 1);
% Minorfracs = (0:(num_MinorTicks-1)) ./ (num_MinorTicks-1);
% MinorTick_idx = round((n-1)*Minorfracs + 1);
% 
% % Major Ticks
% for i = 1:num_MajorTicks
%     range = floor(Major_spike_breite./2);
%     if i==1
%         spike_mat_idx = (MajorTick_idx(i)):(MajorTick_idx(i) + range);
%     elseif i==num_MajorTicks
%         spike_mat_idx = (MajorTick_idx(i) - range):(MajorTick_idx(i));
%     else
%         spike_mat_idx = (MajorTick_idx(i) - range):(MajorTick_idx(i)+ range);
%     end
%     side_padd(spike_mat_idx,:, 1) = color(1);  % Red channel
%     side_padd(spike_mat_idx,:, 2) = color(2);  % Green channel
%     side_padd(spike_mat_idx,:, 3) = color(3);  % Blue channel
% end
% 
% % Minor Ticks
% for i = 1:num_MinorTicks
%     range = floor(Minor_spike_breite./2);
%     if i~=1 && i~=num_MinorTicks
%         spike_mat_idx = (MinorTick_idx(i) - range):(MinorTick_idx(i)+ range);
%     end
%     side_padd(spike_mat_idx,1:floor(Minor_TickScale*pad_size), 1) = color(1);  % Red channel
%     side_padd(spike_mat_idx,1:floor(Minor_TickScale*pad_size), 2) = color(2);  % Green channel
%     side_padd(spike_mat_idx,1:floor(Minor_TickScale*pad_size), 3) = color(3);  % Blue channel
% end
% 
% % Create the edge matrix for the right bottom corner. Then flip for the others
% edge_square_rb = zeros(pad_size*Major_TickScale); 
% tickhalfwidth = floor(Major_spike_breite./2);
% % Red channel
% edge_square_rb(1:pad_size,1:pad_size, 1) = color(1); % Backbone
% edge_square_rb(1:tickhalfwidth,: ,    1) = color(1); % upper tick half
% edge_square_rb(:, 1:tickhalfwidth,    1) = color(1); % left tick half
% % Green channel
% edge_square_rb(1:pad_size,1:pad_size, 2) = color(2);
% edge_square_rb(1:tickhalfwidth,: ,    2) = color(2); 
% edge_square_rb(:, 1:tickhalfwidth,    2) = color(2);
% % Blue channel
% edge_square_rb(1:pad_size,1:pad_size, 3) = color(3);
% edge_square_rb(1:tickhalfwidth,: ,    3) = color(3); 
% edge_square_rb(:, 1:tickhalfwidth,    3) = color(3);
% 
% % Construct the padded image
% padded_image = [ flip(flip(edge_square_rb, 2),1), flip(pagetranspose(side_padd),1), flip(edge_square_rb, 1) ; flip(side_padd,2), Image_RGB, side_padd; flip(edge_square_rb, 2) , pagetranspose(side_padd), edge_square_rb];
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function padded_image = round_padding(Image_RGB, n, padding, color)

% R = 1;
% [x,y] = meshgrid(linspace(-1,1,2*padding));
% r = sqrt(x.^2 + y.^2);
% aa_width = 1 / padding;   % transition thickness (tune this)
% alpha = 1 - min(max((r - R) / aa_width, 0), 1);
% alpha = max(alpha, 0);
% 
% corner_padding(:,:,1) = color(1);
% corner_padding(:,:,2) = color(2);
% corner_padding(:,:,3) = color(3);
% corner_padding = corner_padding .* alpha; 


R = padding;   % radius in pixels

[x,y] = meshgrid(0:2*padding-1, 0:2*padding-1);

% shift so corner center is at (padding, padding)
x = x - padding + 0.5;
y = y - padding + 0.5;

r = sqrt(x.^2 + y.^2);

% distance from content boundary
d = max(R - r, 0);

% normalized gradient (same definition as edges)
g = 1 - min(d / padding, 1);

corner_padding = zeros(2*padding, 2*padding, 3);
corner_padding(:,:,1) = color(1) .* g;
corner_padding(:,:,2) = color(2) .* g;
corner_padding(:,:,3) = color(3) .* g;
corner_padding = flip(flip(corner_padding,1),2);

% % Colored Paddings
% leftright_padding = ones(n,padding, 3);
% leftright_padding(:,:,1) = color(1);
% leftright_padding(:,:,2) = color(2);
% leftright_padding(:,:,3) = color(3);
% topbottom_padding = ones(padding,n, 3);
% topbottom_padding(:,:,1) = color(1);
% topbottom_padding(:,:,2) = color(2);
% topbottom_padding(:,:,3) = color(3);

% Gradient Paddings
gradient = flip((1:padding)./padding);
leftright_padding(:,:,1) = repmat(gradient, n, 1);
leftright_padding(:,:,2) = repmat(gradient, n, 1);
leftright_padding(:,:,3) = repmat(gradient, n, 1);

topbottom_padding(:,:,1) = repmat(gradient', 1, n);
topbottom_padding(:,:,2) = repmat(gradient', 1, n);
topbottom_padding(:,:,3) = repmat(gradient', 1, n);


% Apply the padding
padded_image = [corner_padding(1:padding,1:padding,:),topbottom_padding,corner_padding(1:padding,padding+1:end,:);...
                leftright_padding, Image_RGB, flip(leftright_padding,2);...
                corner_padding(padding+1:end,1:padding,:),flip(topbottom_padding,1),corner_padding(padding+1:end,padding+1:end,:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%