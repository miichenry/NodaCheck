close all
clear

%% Load grid and kernel data 
% (produced by "create_grid_and_kernel_from_stations.m" script)
fid = fopen('base_path.txt', 'r');  % Open the file for reading
folder_name = fscanf(fid, '%s');    % Read the content of the file into folder_name
fclose(fid);                        % Close the file
disp(['Using folder: ', folder_name]);

load([folder_name, '/kernel.mat'])
load([folder_name, '/stat_grid.mat'])
load([folder_name, '/map_matrix.mat'])
disp('.mat files loaded successfully');

%% Define some parameters

% Synthetic model parameters
v_ref = 2.5; % Reference velocity in km/s
v_perturb = 30/100; % Velocity perturbation of the anomalies for checkerboard model use 30/100
nb_el = 1;  % number of elements defining the size of the checkerboard pattern. 
% if nb_el = 2 and grid spacing is 500m, then anomalies are squares of 1 km x 1 km

% Period at which to measure2T = 1; % period in seconds
T = str2double(getenv('T'));
disp(['The period T is: ', num2str(T)]);
%T = 4; % period in seconds
% Station pairs for which their inter-station distance is smaller than
% mul_lambda_thresh multiples of the wavelength apart are excluded (nb
% wavelength is phase velocity times period)
mul_lambda_thresh = 2; %to be more conservative increase this to around 3

if T == 1
    vs_ave = 2.0;
elseif T == 2
    vs_ave = 2.0;
elseif T == 3
    vs_ave = 2.0; %2.5
elseif T == 4
    vs_ave = 2.0; %2.8
elseif T == 5
    vs_ave = 2.0; %3.0
else
    disp('Invalid period value');
end
disp(['Using V = ', num2str(vs_ave)]);

% Load nb_ray from the text file
fid = fopen('nb_ray.txt', 'r');  % Open file for reading
nb_ray = fscanf(fid, '%d');      % Read the value of nb_ray
fclose(fid);                     % Close the file

% Display the loaded value (optional)
disp(['Number of rays (nb_ray) loaded: ', num2str(nb_ray)]);

% Inversion parameters
rel_err = 5/100; % relative error on data. keep whatever, 5%, since it's not real data...
LC = 0.1; % correlation distance. higher values cause more spatial smoothing
sigma = 5; % trade-off parameter between data fit and model complexity. Higher values fit data more.

%% Cut G to station pairs that satisfy distance threshold for given period (multiples of wavelengths)

%vs_ave = 3.0; % Average Vs for the area under study in the upper crust
lambda = vs_ave * T; % wavelength = velocity * period
ind_keep = find(interstation_distance > mul_lambda_thresh*lambda);
G = G_mat(ind_keep, :);
disp(['Number of ray paths that are respect ' num2str(mul_lambda_thresh) '*wavelength distance threshold: ' num2str(size(G,1)) ' out of ' num2str(nb_ray)])

%% Make synthetic checkerboard model

V_2D = zeros(size(X_GRID)); % initialize synthetic velocity model
v_up = (1 + v_perturb) * v_ref; % velocity of positive anomalies
v_down = (1 - v_perturb) * v_ref; % Velocity of negative anomalies
ind_x = 1 : length(x_grid); 
ind_y = 1 : length(y_grid);
bool_x = mod(floor((ind_x - 1)/ nb_el),2) * 2 - 1;
bool_y = mod(floor((ind_y - 1)/ nb_el),2) * (-2) + 1;
bool_check = logical(reshape((bool_x' * bool_y + 1) /2, [1 numel(X_GRID)])); % index of positive anomalies
V_2D(bool_check) = v_up; % Assign positive anomalies
V_2D(~bool_check) = v_down; % Assign positive anomalies
S_2D = 1./V_2D; % slowness synthetic model
S_lin = reshape(S_2D, [numel(X_GRID), 1]); % slowness synthetic model as a vector

%% Calculate synthetic group travel times
TAU = G * S_lin; % TAU contains the group traveltimes between each station pair

%% Perform inversion

d = TAU;  % d = travel time data
N_d = length(d); % number of data points
N_m = size(G,2); % number of model cells

% Define data prior covariance matrix
Cd_vec = (rel_err * d).^2;
CD = diag(Cd_vec);
CD_inv = diag(1 ./ Cd_vec);

% Define model prior
s_prior = 1 / v_ref; % prior homogeneous slowness
m_prior = s_prior * ones(N_m,1);

% Define model prior covariance
L0 = sqrt(dx_grid^2 + dy_grid^2); % size of model cells
x_cell = reshape(X_GRID,[N_m, 1]); y_cell=reshape(Y_GRID,[N_m, 1]);
X_CELL = repmat(x_cell,[1 N_m]); Y_CELL=repmat(y_cell,[1 N_m]);
DIST_CELL = sqrt((X_CELL-X_CELL').^2+(Y_CELL-Y_CELL').^2);
CM = (sigma * L0 / LC)^2 * exp(-1 / LC * DIST_CELL);
CM_inv = inv(CM);

% Perform Tarantola-Valette inversion
m_est = m_prior + (G' * CD_inv * G + CM_inv) \ G' * CD_inv * (d - G * m_prior);
S_map = reshape(m_est,[length(x_grid), length(y_grid)]); % recovered slowness map
V_map = 1./S_map; % recovered velocity map

%% Make a mask for plotting, to hide cells where no ray is sampling
min_density = 1; % require at least 3 ray paths to cross the cell to keep it
thres_dist = 0.01; % distance threshold in km to count a cell as being crossed by ray
G3D = reshape(G_mat',[length(x_grid) length(y_grid) size(G_mat',2)]);
G_count = zeros(size(G3D));
ind_G_ray = G3D(:) > thres_dist; % count ray if >100m in cell
G_count(ind_G_ray) = 1;
G_sum = sum(G_count,3);
mask = zeros(size(G_sum));
mask(G_sum > min_density) = 0.4; % define level of transparency

%% Plot

% in kernel, grid nodes defined at bottom left of cell (should dblcheck); 
% in imagesc, node is at center of cell; these new effective axes compensate for this
x_grid_eff = x_grid + dx_grid/2; 
y_grid_eff = y_grid + dy_grid/2; 


figure(1);
clf;
set(gcf, 'color', 'w','Position', [100, 100, 1600, 600]); % Set figure size (increased width and height)


colormap('jet')

% Adjust the tiled layout with 'loose' padding to create more space at the top
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');  % 'loose' padding for more space

% Plot station map
nexttile;
title(' ', 'FontSize', 1); % Blank title for extra spacing
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'layer', 'top', 'CLim', [v_down, v_up])
hold on; box on
imagesc(x_map, y_map, map); % Plot background
axis equal
set(gca, 'xlim', [min(x_grid), max(x_grid)], 'ylim', [min(y_grid), max(y_grid)]);
plot(x_stat, y_stat, 'vw', 'linewidth', 0.5, 'markersize', 5, 'MarkerFaceColor', 'k')

cb = colorbar;
ylabel(cb, 'Velocity [km/s]', 'FontSize', 14)

title('Station Layout', 'FontSize', 14)

% Plot synthetic model
nexttile;
title(' ', 'FontSize', 1); % Blank title for extra spacing
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'layer', 'top', 'CLim', [v_down, v_up])
hold on; box on
imagesc(x_map, y_map, map); % Plot background
axis equal
set(gca, 'xlim', [min(x_grid), max(x_grid)], 'ylim', [min(y_grid), max(y_grid)]);
im = pcolor(x_grid_eff, y_grid_eff, V_2D'); % plot data
set(im, 'facealpha', 'flat', 'alphadata', mask'); % apply mask
shading('flat');
plot(x_stat, y_stat, 'vw', 'linewidth', 0.5, 'markersize', 5, 'MarkerFaceColor', 'k')

cb = colorbar;
ylabel(cb, 'Velocity [km/s]', 'FontSize', 14)
title('Synthetic model', 'FontSize', 14)

% Plot recovered model
nexttile;
title(' ', 'FontSize', 1); % Blank title for extra spacing
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'layer', 'top', 'CLim', [v_down, v_up])
hold on; box on
imagesc(x_map, y_map, map); % Plot background
axis equal
set(gca, 'xlim', [min(x_grid), max(x_grid)], 'ylim', [min(y_grid), max(y_grid)]);
im = pcolor(x_grid_eff, y_grid_eff, V_map'); % plot data
set(im, 'facealpha', 'flat', 'alphadata', mask'); % apply mask
shading('flat');
plot(x_stat, y_stat, 'vw', 'linewidth', 0.5, 'markersize', 2, 'MarkerFaceColor', 'k')

cb = colorbar;
ylabel(cb, 'Velocity [km/s]', 'FontSize', 14)
title('Recovered model', 'FontSize', 14)

%sgtitle(['Checkerboard test at T = ' num2str(T) ' s'], 'FontSize', 18)i
sgtitle(sprintf('Checkerboard test at T = %ds\n', T));

% Save with enhanced resolution
fig_fname = [folder_name, '/checkerboard_T' num2str(T) '_mul' num2str(mul_lambda_thresh) '_' num2str(nb_el) 'x' num2str(nb_el) '_res' num2str(dx_grid) '.png'];
print(gcf, fig_fname, '-dpng', '-r300') % Save with 300 dpi resolution
disp(['Saved figure to: ' fig_fname])
close(gcf);

figure;
set(gcf,'color','w'); % Set background color to white
set(gcf, 'Position', [100, 100, 800, 800]); % Set figure size for better resolution (square size)

% Adjust axes properties
set(gca,'linewidth',1.5,'fontsize',14,'layer','top','CLim',[v_down, v_up])
hold on ; box on
imagesc(x_map,y_map,map) % Plot background
colormap('jet'); % Use the same colormap

% Ensure equal axis scaling and adjust limits
axis equal
axis tight
set(gca,'xlim',[min(x_grid), max(x_grid)],'ylim',[min(y_grid), max(y_grid)]);

% Plot data
im = pcolor(x_grid_eff,y_grid_eff,V_map'); % plot data
set(im,'facealpha','flat','alphadata',mask') % apply mask
shading('flat');
plot(x_stat,y_stat,'vw','linewidth',0.5,'markersize',5, 'MarkerFaceColor','k')

% Add colorbar
cb = colorbar;
ylabel(cb, 'Velocity [km/s]')

% Add title
%title('Recovered model')
filename = fullfile(folder_name, ['recovered_model_T_', num2str(T), '_res', num2str(dx_grid), '.png']);
exportgraphics(gcf, filename, 'Resolution', 300);  % Save as PNG or other desired format
disp(['Saved figure to: ', filename]);
% Save the figure in high resolution
%print('Recovered_model_high_res', '-dpng', '-r300') % Save as PNG with 300 dpi resolution
close(gcf);

% After plotting the synthetic model
% Export the synthetic model matrix
fileName1 = ['synthetic_model', '_T_', num2str(T), '_res', num2str(dx_grid)];
fullPath1 = fullfile(folder_name, [fileName1, '.mat']);
disp(['Saved mat to: ' fullPath1])
% or save as .mat file
save(fullPath1, 'V_2D', 'x_grid_eff', 'y_grid_eff');

% After plotting the recovered model
% Export the recovered model matrix
fileName2 = ['recovered_model', '_T_', num2str(T), '_res', num2str(dx_grid)];
fullPath2 = fullfile(folder_name, [fileName2, '.mat']);
disp(['Saved mat to: ' fullPath2])
% or save as .mat file
save(fullPath2, 'V_map', 'x_grid_eff', 'y_grid_eff');
