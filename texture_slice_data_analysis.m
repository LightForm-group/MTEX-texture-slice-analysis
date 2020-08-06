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

% give the nominal values
nominal_temperature = 850
nominal_strain = 0.69
nominal_strain_rate = -0.1

% calculate normalised values
FE_temperature_norm = ( FE_temperature - min(FE_temperature) ) / ( max(FE_temperature) -  min(FE_temperature) )
FE_strain_radial_norm = ( FE_strain_radial  - min(FE_strain_radial) ) / ( max(FE_strain_radial) -  min(FE_strain_radial) )
FE_strain_axial_norm = ( FE_strain_axial  - min(FE_strain_axial) ) / ( max(FE_strain_axial) -  min(FE_strain_axial) )
FE_strain_equiv_norm = ( FE_strain_equiv - min(FE_strain_equiv) ) / ( max(FE_strain_equiv) -  min(FE_strain_equiv) )

FE_strain_rate_positive = sqrt(FE_strain_rate.*FE_strain_rate)
FE_strain_rate_norm = ( FE_strain_rate_positive - min(FE_strain_rate_positive) ) / ( max(FE_strain_rate_positive) -  min(FE_strain_rate_positive) )

% calculate temperature difference
FE_temperature_diff = ( FE_temperature - min(FE_temperature) )+( max(FE_temperature) - min(FE_temperature) )

%% Load TEXTURE results

% import the TEXTURE data
TEXTURE = importdata(fname_TEXTURE_results)

% write the data to new arrays
TEXTURE_strip_index = flip(TEXTURE.data(:,1)*1000) % flip only distance to swop
TEXTURE_TI = TEXTURE.data(:,2)
TEXTURE_ODF_max = TEXTURE.data(:,3)
TEXTURE_PHI = TEXTURE.data(:,4)
TEXTURE_misori = TEXTURE.data(:,5)
TEXTURE_basal_ND = TEXTURE.data(:,6)
TEXTURE_tilted_ND = TEXTURE.data(:,7)

% calculate the normalised values
TEXTURE_TI_norm = ( TEXTURE_TI - min(TEXTURE_TI) ) / ( max(TEXTURE_TI) -  min(TEXTURE_TI) )
TEXTURE_basal_ND_norm = ( TEXTURE_basal_ND - min(TEXTURE_basal_ND) ) / ( max(TEXTURE_basal_ND) -  min(TEXTURE_basal_ND) )
TEXTURE_PHI_norm = ( TEXTURE_PHI - min(TEXTURE_PHI) ) / ( max(TEXTURE_PHI) -  min(TEXTURE_PHI) )

%% Get the maxima and minima values

% FE maxima
max_FE_temperature = max(FE_temperature)
max_FE_strain_equiv = max(FE_strain_equiv)
max_FE_strain_rate = max(FE_strain_rate)

% FE minima
min_FE_temperature = min(FE_temperature)
min_FE_strain_equiv = min(FE_strain_equiv)
min_FE_strain_rate = min(FE_strain_rate)

% TEXTURE maxima
max_TEXTURE_TI = max(TEXTURE_TI)
max_TEXTURE_PHI = max(TEXTURE_PHI)*( 180/pi() )
max_TEXTURE_basal_ND = max(TEXTURE_basal_ND)*100 %for percent

%TEXTURE minima
min_TEXTURE_TI = min(TEXTURE_TI)
min_TEXTURE_PHI = min(TEXTURE_PHI)*( 180/pi() )
min_TEXTURE_basal_ND = min(TEXTURE_basal_ND)*100 %for percent

%% [TODO] - Get global maxima and minima and adjust the values

%% Create an array of positions to plot the strips against

% Give the maximum sample length in mm
max_sample_length = 5.2

num_strips = length(TEXTURE_strip_index)
strip_width = max_sample_length / num_strips
for strip_index = 1:num_strips
    strip_position(strip_index) = 0.5*strip_width + (strip_index - 1)*strip_width
end

%% Plot the texture variation and the FE results

% setup the figure
texture_variation_FE_results = figure();
set(gcf, 'Position', [10 10 900 1000]) % format is [left bottom width height]
set(gca,'XAxisLocation','top','YAxisLocation','left','yDir','reverse', 'Fontsize', 16, 'lineWidth',2);

% plot the FE results
line(FE_strain_equiv_norm,FE_distance,'lineWidth',2);
xlabel('Normalised Strain, Strain Rate and Temperature');
hold on
line(FE_strain_rate_norm,FE_distance,'lineWidth',2);
hold on
line(FE_temperature_norm,FE_distance, 'lineWidth',2);
hold on

% plot the texture variation
line(TEXTURE_TI_norm,strip_position,'Marker','o','lineWidth',2,'lineStyle','none');
hold on
xlabel('Normalised Texture Values');
ylabel('Position (mm)');
plot(TEXTURE_basal_ND_norm,strip_position,'Marker','o','lineWidth',2,'lineStyle','none');
hold on
plot(TEXTURE_PHI_norm,strip_position,'Marker','o','lineWidth',2,'lineStyle','none');
hold on

% adapt this for plotting another dataset on different axes
% ax1 = gca; % current axes
% ax1.XColor = 'k';
% ax1.YColor = 'k';
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,'XAxisLocation','bottom','Color','none', 'yDir', 'reverse');
% xlabel('Texture Values')
% hold on
% line(FE_strain_equiv_norm,FE_distance,'Parent',ax2,'lineWidth',2);

% add a legend
legend({'Equivalent Strain','Strain Rate','Temperature','Texture Index','Basal ND','Phi Angle of ODF Maxima'}, 'Location', [0.75,0.45,0,0])

% adding axes for normalised parameters
ax1 = gca; % current axes
ax1.XColor = 'none'; % clear previous x-axis
ax1_pos = ax1.Position; % position of first axes
set(ax1, 'Position', ax1_pos + [0 0.1 -0.3 -0.3]); % move existing axis up a bit and reduce height, format is [left bottom width height]
ax2 = axes('position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.10 -0.3 0], 'color', 'none', 'linewidth', 2);
ax3 = axes('position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.025 -0.3 0], 'color', 'none', 'linewidth', 2);
ax4 = axes('position', (ax1_pos .* [1 1 1 1e-3]) + [0 -0.05 -0.3 0], 'color', 'none', 'linewidth', 2);
ax5 = axes('XAxisLocation','top','position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.615 -0.3 0], 'color', 'none', 'linewidth', 2);
ax6 = axes('XAxisLocation','top','position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.69 -0.3 0], 'color', 'none', 'linewidth', 2);
ax7 = axes('XAxisLocation','top','position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.765 -0.3 0], 'color', 'none', 'linewidth', 2);
set(ax1, 'ylim', [0 max_sample_length]);
set(ax2, 'xlim', [min_TEXTURE_TI max_TEXTURE_TI], 'Fontsize', 16);
set(ax3, 'xlim', [min_TEXTURE_basal_ND max_TEXTURE_basal_ND], 'Fontsize', 16);
set(ax4, 'xlim', [min_TEXTURE_PHI max_TEXTURE_PHI], 'Fontsize', 16);
set(ax5, 'xlim', [min_FE_strain_rate max_FE_strain_rate], 'Fontsize', 16);
set(ax6, 'xlim', [min_FE_strain_equiv max_FE_strain_equiv], 'Fontsize', 16);
set(ax7, 'xlim', [min_FE_temperature max_FE_temperature], 'Fontsize', 16);
axes(ax2); xlabel('Texture Index');
axes(ax3); xlabel('0,0,0 Component (%)');
axes(ax4); xlabel('Phi Angle (degree)');
axes(ax5); xlabel('Strain Rate');
axes(ax6); xlabel('Strain');
axes(ax7); xlabel('Temperature');

hold off

saveas (texture_variation_FE_results, strcat(analysis_path,sample_name,'_texture_variation_FE_results_',num2str(num_strips), '.bmp'));