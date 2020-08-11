%% Specify the file names

% path to files
pname = '/Users/mbcx9cd4/Documents/MATLAB/ebsd/MTEX-texture-slice-analysis/';

% where are the FE results located
mech_tester = 'Dilatometer'
data_folder = '/Example Data/'
sample_name = 'D5' % edit this line to name of file
FE_results_file = 'FE Results D5.txt'
FE_results_path = strcat(mech_tester,data_folder,sample_name,'/',FE_results_file)
fname_FE_results = [pname FE_results_path]

% where should the analysis be saved
analysis_folder = '/Example Analysis/'
analysis_path = strcat(mech_tester,analysis_folder,sample_name,'/') % path for saving the data

% where is the texture data from the slice analysis saved
TEXTURE_results_file = 'D5_texture_strength_20.txt'
TEXTURE_results_path = strcat(mech_tester,analysis_folder,sample_name,'/',TEXTURE_results_file)
fname_TEXTURE_results = [pname TEXTURE_results_path]

%% Load FE model results

% import the FE data
FE_model = importdata(fname_FE_results)

% write the data to new arrays
FE_distance = flip(FE_model.data(:,1)*1000) % flip only distance to swop
FE_temperature = FE_model.data(:,2)
FE_strain_radial = FE_model.data(:,3)
FE_strain_axial = FE_model.data(:,4)
FE_strain_equiv = FE_model.data(:,5)
FE_strain_rate = FE_model.data(:,6)

% calculate a positive value for the strain rate
FE_strain_rate_positive = sqrt(FE_strain_rate.*FE_strain_rate)

% give the nominal values
nominal_temperature = 850
nominal_strain = 0.69
nominal_strain_rate = -0.1

% calculate temperature difference
FE_temperature_diff = ( FE_temperature - min(FE_temperature) )+( max(FE_temperature) - min(FE_temperature) )

%% Load TEXTURE results

% import the TEXTURE data
TEXTURE = importdata(fname_TEXTURE_results)

% write the data to new arrays
TEXTURE_strip_index = flip(TEXTURE.data(:,1)*1000) % flip only distance to swop
TEXTURE_TI = TEXTURE.data(:,2)
TEXTURE_ODF_max = TEXTURE.data(:,3)
TEXTURE_PHI = TEXTURE.data(:,4).*( 180/pi() ) %for degrees
TEXTURE_misori = TEXTURE.data(:,5)
TEXTURE_basal_ND = TEXTURE.data(:,6).*100 %for percent
TEXTURE_tilted_ND = TEXTURE.data(:,7).*100 %for percent

%% Get the maxima and minima values

% FE maxima
max_FE_temperature = max(FE_temperature)
max_FE_strain_equiv = max(FE_strain_equiv)
max_FE_strain_rate = max(FE_strain_rate)
max_FE_strain_rate_positive = max(FE_strain_rate_positive)

% FE minima
min_FE_temperature = min(FE_temperature)
min_FE_strain_equiv = min(FE_strain_equiv)
min_FE_strain_rate = min(FE_strain_rate)
min_FE_strain_rate_positive = min(FE_strain_rate_positive)

% TEXTURE maxima
max_TEXTURE_TI = max(TEXTURE_TI)
max_TEXTURE_PHI = max(TEXTURE_PHI)
max_TEXTURE_basal_ND = max(TEXTURE_basal_ND)

%TEXTURE minima
min_TEXTURE_TI = min(TEXTURE_TI)
min_TEXTURE_PHI = min(TEXTURE_PHI)
min_TEXTURE_basal_ND = min(TEXTURE_basal_ND)

%% Or, give global maxima and minima to adjust the values

% FE maxima
max_FE_temperature = 1
max_FE_strain_equiv = 1
max_FE_strain_rate = 1
max_FE_strain_rate_positive = 1

% FE minima
min_FE_temperature = 1
min_FE_strain_equiv = 1
min_FE_strain_rate = 1
min_FE_strain_rate_positive = 1

% TEXTURE maxima
max_TEXTURE_TI = 1
max_TEXTURE_PHI = 1*( 180/pi() )
max_TEXTURE_basal_ND = 1*100 %for percent

%TEXTURE minima
min_TEXTURE_TI = 1
min_TEXTURE_PHI = 1*( 180/pi() )
min_TEXTURE_basal_ND = 1*100 %for percent

%% Calculate normalised values

% calculate the normalised values for texture
TEXTURE_TI_norm = ( TEXTURE_TI - min_TEXTURE_TI ) / ( max_TEXTURE_TI -  min_TEXTURE_TI )
TEXTURE_basal_ND_norm = ( TEXTURE_basal_ND - min_TEXTURE_basal_ND ) / ( max_TEXTURE_basal_ND -  min_TEXTURE_basal_ND )
TEXTURE_PHI_norm = ( TEXTURE_PHI - min_TEXTURE_PHI ) / ( max_TEXTURE_PHI -  min_TEXTURE_PHI )

% calculate normalised values
FE_temperature_norm = ( FE_temperature - min_FE_temperature ) / ( max_FE_temperature -  min_FE_temperature )
FE_strain_equiv_norm = ( FE_strain_equiv - min_FE_strain_equiv ) / ( max_FE_strain_equiv -  min_FE_strain_equiv )
FE_strain_rate_norm = ( FE_strain_rate_positive - min_FE_strain_rate_positive ) / ( max_FE_strain_rate_positive -  min_FE_strain_rate_positive )

%% Create an array of positions to plot the strips against

% Give the maximum sample length in mm
% max_sample_length = max(FE_distance)
% max_sample_length = 5.2 %dilatometer 5.2205
% max_sample_length = 7.8 %servotest 7.8303
max_sample_length = 8.2 %hydrawedge 8.2323

num_strips = length(TEXTURE_strip_index)
strip_width = max_sample_length / num_strips
for strip_index = 1:num_strips
    strip_position(strip_index) = 0.5*strip_width + (strip_index - 1)*strip_width
end

%% Fit texture variation using polyfit or using weighted polyfit

% polyfitweighted.m downloaded from - https://uk.mathworks.com/matlabcentral/fileexchange/13520-polyfitweighted

% weights for 20 positions, put 0 on excluded points, 1 on used points
weight_TI=[1,1,0,0,1,0,1,0,1,0,1,0,1,1,1,0,0,0,1,0] % weight of points from top of plot down i.e. position 0 -> x mm
weight_basal_ND = [0,1,0,1,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,0]
weight_PHI = [1,1,1,1,1,0,1,0,10,1,1,1,0,1,1,1,1,1,0,1]

% Fit the texture variation (note, have to swop the x and y to fit as polynomial)
% poly_texture_TI = polyfit(transpose(strip_position),TEXTURE_TI_norm,3)
poly_texture_TI = polyfitweighted(transpose(strip_position),TEXTURE_TI_norm,3,transpose(weight_TI)) % use this to weight fit
y2_texture_TI = min(strip_position):.01:max(strip_position)
x2_texture_TI = polyval(poly_texture_TI,y2_texture_TI);

% Fit the texture variation (note, have to swop the x and y to fit as polynomial)
% poly_texture_basal_ND = polyfit(transpose(strip_position),TEXTURE_basal_ND_norm,3)
poly_texture_basal_ND = polyfitweighted(transpose(strip_position),TEXTURE_basal_ND_norm,3,transpose(weight_basal_ND)) % use this to weight fit
y2_basal_ND = min(strip_position):.01:max(strip_position)
x2_basal_ND = polyval(poly_texture_basal_ND,y2_basal_ND);

% Fit the texture variation (note, have to swop the x and y to fit as polynomial)
% poly_texture_PHI = polyfit(transpose(strip_position),TEXTURE_PHI_norm,3)
poly_texture_PHI = polyfitweighted(transpose(strip_position),TEXTURE_PHI_norm,3,transpose(weight_PHI)) % use this to weight fit
y2_texture_PHI = min(strip_position):.01:max(strip_position)
x2_texture_PHI = polyval(poly_texture_PHI,y2_texture_PHI);

%% Plot the texture variation and the FE results

% setup the figure
texture_variation_FE_results = figure();
set(gcf, 'Position', [10 10 900 1000]) % format is [left bottom width height]
set(gca,'XAxisLocation','top','YAxisLocation','left','yDir','reverse', 'Fontsize', 16, 'lineWidth',2);

% plot the FE results
line(FE_temperature_norm,FE_distance, 'lineWidth',3, 'Color', [1 0 0]);
xlabel('Normalised Strain, Strain Rate and Temperature');
ylabel('Position (mm)');
hold on
line(FE_strain_equiv_norm,FE_distance,'lineWidth',3, 'Color', [1 0.5 0]);
hold on
line(FE_strain_rate_norm,FE_distance,'lineWidth',3,'Color', [1 0.8 0] );
hold on

% plot the texture variation
line(TEXTURE_TI_norm,strip_position,'Marker','o','Markersize', 10,'lineWidth',3,'lineStyle','none', 'Color', [0 0 1]);
hold on
xlabel('Normalised Texture Values');
ylabel('Position (mm)');
line(TEXTURE_basal_ND_norm,strip_position,'Marker','x','Markersize', 10,'lineWidth',3,'lineStyle','none', 'Color', [0.75 0 1]);
hold on
line(TEXTURE_PHI_norm,strip_position,'Marker','+','Markersize', 10,'lineWidth',3,'lineStyle','none', 'Color', [0 0.75 1]);
hold on

% Plot the fit or weighted fit of the texture variation
% line(x2_texture_TI,y2_texture_TI,'lineWidth',1,'lineStyle','-', 'Color', [0 0 1]);
% hold on
% line(x2_basal_ND,y2_basal_ND,'lineWidth',1,'lineStyle','-', 'Color', [0.75 0 1]);
% hold on
% line(x2_texture_PHI,y2_texture_PHI,'lineWidth',1,'lineStyle','-', 'Color', [0 0.75 1]);
% hold on

% add a legend
legend({'Temperature ({\circC})','Strain','Strain Rate ({s^-^1})','Texture Index','0{\circ},0{\circ},0{\circ} Component (%)','Phi Angle ({\circ})'}, 'Location', [0.75,0.45,0,0])

% adding axes for normalised parameters
ax1 = gca; % current axes
ax1.XColor = 'none'; % clear previous x-axis
ax1_pos = ax1.Position; % position of first axes
set(ax1, 'Position', ax1_pos + [0 0.1 -0.3 -0.3]); % move existing axis up a bit and reduce height, format is [left bottom width height]
ax2 = axes('position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.10 -0.3 0], 'color', 'none', 'linewidth', 2);
ax3 = axes('position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.025 -0.3 0], 'color', 'none', 'linewidth', 2);
ax4 = axes('position', (ax1_pos .* [1 1 1 1e-3]) + [0 -0.05 -0.3 0], 'color', 'none', 'linewidth', 2);
ax5 = axes('XAxisLocation','top','position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.615 -0.3 0], 'color', 'none', 'linewidth', 2, 'xDir', 'reverse');
ax6 = axes('XAxisLocation','top','position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.69 -0.3 0], 'color', 'none', 'linewidth', 2, 'xDir', 'reverse');
ax7 = axes('XAxisLocation','top','position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.765 -0.3 0], 'color', 'none', 'linewidth', 2);
set(ax1, 'ylim', [0 max_sample_length]);
set(ax2, 'xlim', [min_TEXTURE_TI max_TEXTURE_TI], 'Fontsize', 16);
set(ax3, 'xlim', [min_TEXTURE_basal_ND max_TEXTURE_basal_ND], 'Fontsize', 16);
set(ax4, 'xlim', [min_TEXTURE_PHI max_TEXTURE_PHI], 'Fontsize', 16);
set(ax5, 'xlim', [min_FE_strain_rate max_FE_strain_rate], 'Fontsize', 16);
set(ax6, 'xlim', [-max_FE_strain_equiv -min_FE_strain_equiv], 'Fontsize', 16);
set(ax7, 'xlim', [min_FE_temperature max_FE_temperature], 'Fontsize', 16);
axes(ax2); xlabel('Texture Index');
axes(ax3); xlabel('0{\circ},0{\circ},0{\circ} Component (%)');
axes(ax4); xlabel('Phi Angle ({\circ})');
axes(ax5); xlabel('Strain Rate ({s^-^1})');
axes(ax6); xlabel('Strain');
axes(ax7); xlabel('Temperature ({\circC})');

hold off

saveas (texture_variation_FE_results, strcat(analysis_path,sample_name,'_texture_variation_FE_results_scatter_',num2str(num_strips), '.bmp'));

%% Notes

% adapt this for plotting another dataset on different axes (not used currently)
% ax1 = gca; % current axes
% ax1.XColor = 'k';
% ax1.YColor = 'k';
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,'XAxisLocation','bottom','Color','none', 'yDir', 'reverse');
% xlabel('Texture Values')
% hold on
% line(FE_strain_equiv_norm,FE_distance,'Parent',ax2,'lineWidth',2);

max(FE_distance)