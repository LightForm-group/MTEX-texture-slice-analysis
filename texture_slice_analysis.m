%% Specify Crystal and Specimen Symmetries

% crystal symmetry with alpha first
% CS = {... 
%    'notIndexed',...
%    crystalSymmetry('6/mmm', [2.954 2.954 4.729], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ti-Hex', 'color', 'light blue'),...
%    crystalSymmetry('m-3m', [3.192 3.192 3.192], 'mineral', 'Titanium cubic', 'color', 'light green')};
% SS = specimenSymmetry('triclinic')

% crystal symmetry with beta first
CS = {... 
    'notIndexed',...
    crystalSymmetry('m-3m', [3.192 3.192 3.192], 'mineral', 'Titanium cubic', 'color', 'light green'),...
    crystalSymmetry('6/mmm', [2.954 2.954 4.729], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ti-Hex', 'color', 'light blue')};

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','intoPlane');

%% Specify File Names

% path to files
pname = '/Users/mbcx9cd4/Documents/MATLAB/ebsd/MTEX-texture-slice-analysis/';

% which files to be imported
mech_tester = 'Dilatometer'
data_file = '/Example Data/'
sample_name = 'D5' % edit this line to name of file
ctf_file = '/sampleD5_repeat Data.ctf'
data_path = strcat(mech_tester,data_file,sample_name,ctf_file)

analysis_folder = '/Example Analysis/'
analysis_path = strcat(mech_tester,analysis_folder,sample_name,'/') % path for saving the data

% which files to be imported
fname = [pname data_path]

%% Specify File Names

% path to files
pname = '/Users/mbcx9cd4/Documents/MATLAB/ebsd/MTEX-texture-slice-analysis/';

% which files to be imported
mech_tester = 'Dilatometer'
data_file = '/Data/'
sample_name = 'Chris D2' % edit this line to name of file
ctf_file = '/Chris D2 High Res Line Scan.ctf'
data_path = strcat(mech_tester,data_file,sample_name,ctf_file)

analysis_folder = '/Analysis/'
analysis_path = strcat(mech_tester,analysis_folder,sample_name,'/') % path for saving the data

% which files to be imported
fname = [pname data_path]

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame')

%% Reduce the Data

ebsd = reduce(ebsd,5)

%% Rotate the data

rot = rotation('Euler', 90*degree, 90*degree, 0*degree);
% rot = rotation('Euler', 0*degree, 90*degree, 90*degree);

ebsd = rotate(ebsd,rot,'keepXY'); % rotate the orientation data
ebsd = rotate(ebsd,90*degree,'keepEuler') % rotate the spatial data

% ebsd = rotate(ebsd,rot); % or we could just use this to rotate the map as well

fprintf('Note, the (x,y) origin on the map will have changed and x or y could be negative!')

%% Set whether figures are visible

visible = 'on'

%% Plot the IPF colour map

phase = 'alpha'
outputFileName = strcat(analysis_path,sample_name,'_IPF_map_entire_region')
IPF_map_plot(phase, ebsd, outputFileName, visible)

%% Calculate Sample Dimensions

ebsd_grid = ebsd.gridify;
ebsd_shape = size(ebsd_grid.id);
original_y = ebsd_shape(1);
original_x = ebsd_shape(2);
stepSize = ebsd_grid.dx;
size_y = ebsd_shape(1)*stepSize
size_x = ebsd_shape(2)*stepSize

fprintf('Sample size limits are, y limit =', size_y,' and x limit = ', size_x)

%% Plot the pole figures for the whole compression sample

phase = 'alpha'
ori = ebsd('Ti-Hex').orientations
contour_step = 0.1
pf_max = 2.5
outputFileName = strcat(analysis_path,sample_name,'_PF_entire_region')
pole_figure_plot(phase, ori, CS, contour_step, pf_max, outputFileName, visible);

%% Calculate an automatic halfwidth for the ODF 

% if the orientations are spatially independant...
psi=calcKernel(ebsd('Ti-Hex').orientations);

% if the EBSD measurements are not rough texture measurements and are spatially dependant (more than one point per grain), then perform grain reconstruction and estimate the halfwidth from the grains...
% grain reconstruction (default is 10 degrees so could just use calcGrains(ebsd))...
% grains = calcGrains(ebsd,'angle',10*degree);
% correct for small grains...
% grains = grains(grains.grainSize>5);
% compute optimal halfwidth from the meanorientations of grains...
% psi = calcKernel(grains('Zirc-alloy4').meanOrientation)

HALF_WIDTH = psi

%% Calculate the ODF using the optimised halfwidth

odf = calcDensity(ebsd('Ti-Hex').orientations,'kernel',psi);

%% Calculate the Texture Index for the whole compresion sample

TEXTURE_INDEX = textureindex(odf)

%% Plot the ODF slices without contouring for the whole compression sample

odf_max = 3.0
outputFileName = strcat(analysis_path,sample_name,'_ODF_entire_region')
specSym = 'triclinic'
ODF_plot(phase, odf, odf_max, outputFileName, specSym, visible)

%% Crop the map into a rectangle (avoiding the barrelled edges of the compression sample)
% note, in this case x is vertical starting at 0 and runs downwards as
% negative values, y is horizontal starting at 0 and runs left-to-right as
% positive values. 

x_top = -50
y_left = -216
x_bottom = -5100
y_right = 2868

x_width = x_bottom-x_top
y_width = y_right-y_left

region = [x_top, y_left, x_width, y_width]; % note, region is defined as x,y origin and an x,y width which is added onto the origin
condition = inpolygon(ebsd,region); % points located within region
ebsd_cropped = ebsd(condition);
ori_cropped = ebsd_cropped('Ti-Hex').orientations

% plot the IPF map for the cropped compression sample
outputFileName = strcat(analysis_path,sample_name,'_IPF_map_cropped')
IPF_map_plot(phase, ebsd_cropped, outputFileName, visible)

% plot the pole figures for the cropped compression sample
outputFileName = strcat(analysis_path,sample_name,'_PF_cropped')
pole_figure_plot(phase, ori_cropped, CS, contour_step, pf_max, outputFileName, visible);

% calculate an automatic halfwidth for the ODF for the cropped compression sample
psi=calcKernel(ori_cropped);
HALF_WIDTH = psi

% calculate the ODF using the optimised halfwidth for the cropped compression sample
odf_cropped = calcDensity(ori_cropped,'kernel',psi);

% calculate the Texture Index for the cropped compression sample
TEXTURE_INDEX = textureindex(odf_cropped)

% plot the ODF slices without contouring for the cropped compression sample
outputFileName = strcat(analysis_path,sample_name,'_ODF_cropped')
specSym = 'triclinic'
ODF_plot(phase, odf_cropped, odf_max, outputFileName, specSym, visible)

%% Cropping dimensions for different samples

% sample D2 on Dilatometer
%x_top = -20
%y_left = 2700
%x_bottom = -3900
%y_right = 5200

% sample D3 on Dilatometer
%x_top = -200
%y_left = 2460
%x_bottom = -5300
%y_right = 6140

% sample D5 on Dilatometer
%x_top = -250
%y_left = 2400
%x_bottom = -5100
%y_right = 5600

% sample D9 on Dilatometer
%x_top = -140
%y_left = 2340
%x_bottom = -4770
%y_right = 5570

% sample S1 on Servotester
%x_top = -180
%y_left = 2100
%x_bottom = -7000
%y_right = 6330

% sample S2 on Servotester
%x_top = -40
%y_left = 1600
%x_bottom = -7000
%y_right = 6330

% sample S3 on Servotester
%x_top = -50
%y_left = 2600
%x_bottom = -7190
%y_right = 7600

% sample S4 on Servotester
%x_top = -20
%y_left = 1355
%x_bottom = -6800
%y_right = 6090

% sample S5 on Servotester
%x_top = -40
%y_left = 3600
%x_bottom = -7300
%y_right = 8700

% sample S6 on Servotest
%x_top = -40
%y_left = 2000
%x_bottom = -7200
%y_right = 7000

% sample S7 on Servotest
%x_top = -30
%y_left = 2430
%x_bottom = -6420
%y_right = 6360

% sample S8 on Servotest
%x_top = -80
%y_left = 2040
%x_bottom = -7730
%y_right = 6540

% sample S9 on Servotest
%x_top = -30
%y_left = 1900
%x_bottom = -6960
%y_right = 6830

% sample H1 on Hydrawedge
%x_top = -40
%y_left = 140
%x_bottom = -7800
%y_right = 7100

% sample H2 on Hydrawedge
%x_top = -330
%y_left = 4800
%x_bottom = -7700
%y_right = 11300

% sample H3 on Hydrawedge
%x_top = -60
%y_left = 580
%x_bottom = -8100
%y_right = 6950

% sample H4 on Hydrawedge
%x_top = -50
%y_left = 1800
%x_bottom = -7500
%y_right = 6200

% sample H5 on Hydrawedge
%x_top = -30
%y_left = 800
%x_bottom = -6730
%y_right = 6520

% sample H6 on Hydrawedge
%x_top = -70
%y_left = 1655
%x_bottom = -7820
%y_right = 7300

%% Choose whether to slice the full map or the cropped map here

% use for cropped map
ebsd = ebsd_cropped
x_origin = x_top
y_origin = y_left

% use for entire map
% ebsd = ebsd
% x_origin = 0
% y_origin = 0

%% Slice the entire map or the cropped map

num_strips = 30; % number of strips to cut the map into (resolution)

% define the size of the EBSD map
ebsd_grid = ebsd.gridify;
ebsd_shape = size(ebsd_grid.id);
original_y = ebsd_shape(1);
original_x = ebsd_shape(2);
stepSize = ebsd_grid.dx;

x_min = (sqrt(x_origin * x_origin)/stepSize)
x_max = original_x + (sqrt(x_origin * x_origin)/stepSize)
x_length = x_max - x_min

y_min = (sqrt(y_origin * y_origin)/stepSize)
y_max = original_y + (sqrt(y_origin * y_origin)/stepSize)
y_length = y_max - y_min

% used if splitting into strips along y
y_width = floor(y_length / num_strips) % round to nearest integer
y_axis = (1:num_strips)

% used if splitting into strips along x
x_width = floor(x_length / num_strips) % round to nearest integer
x_axis = (1:num_strips)

cutmap = containers.Map('KeyType', 'int32', 'ValueType', 'any'); % creates an empty Map object

for strip_index = 0:num_strips-1
    % separate the map section
    
    % set out the coordinates for the edge of the region
    % note, region is defined as x,y origin and an x,y width which is added onto the origin
    
    % if splitting into strips along y (breaking up y)
    % y_min_strip = strip_index * y_width + y_min;
    % region = [x_min*stepSize, y_min_strip*stepSize, x_length*stepSize, y_width*stepSize];
    
    % if splitting into strips along y (breaking up y) and x is negative
    % y_min_strip = strip_index * y_width + y_min;
    % region = [-x_min*stepSize, y_min_strip*stepSize, -x_length*stepSize, y_width*stepSize];
    
    % if splitting into strips along x (breaking up x)
    % x_min_strip = strip_index * x_width + x_min;
    % region = [x_min_strip*stepSize, y_min*stepSize, x_width*stepSize, y_length*stepSize];
    
    % if splitting into strips along x (breaking up x) and x is negative
    % and x is moving in negative direction
    % x_min_strip = strip_index * x_width + x_min;
    % region = [-x_min_strip*stepSize, y_min*stepSize, -x_width*stepSize, y_length*stepSize];
    
    % if splitting into strips along x (breaking up x) and x and y are negative
    % and x is moving in negative direction, 
    x_min_strip = strip_index * x_width + x_min;
    region = [-x_min_strip*stepSize, -y_min*stepSize, -x_width*stepSize, y_length*stepSize];
    
    % Cut the EBSD map
    condition = inpolygon(ebsd,region); % points located within region
    ebsd_strip = ebsd(condition); % create ebsd map for region
    cutmap(strip_index) = ebsd_strip; % store strip in Map object with index
    ebsd_cutmap = cutmap(strip_index); % read out ebsd_cutmap from the Map object
    
    % plot the IPF map to check the slices
    outputFileName = strcat(analysis_path,sample_name,'_IPF_map_strip_',num2str(strip_index))
    IPF_map_plot(phase, ebsd_cutmap, outputFileName, visible)
    
end 

%% Analyse and plot the sliced data to see how the texture components change along the length

% define the crystal system for the texture components
cs = ebsd('Ti-Hex').CS;

% define the maximum possible misorientation
misorientation = 10

% Define a texture component for the hexagonal phase
% Define a texture component for the hexagonal phase
% basal_ND = symmetrise(orientation.byMiller([0 0 0 1],[1 0 -1 0],cs),'unique'); % define component with directions
basal_ND = symmetrise(orientation.byEuler(0*degree,0*degree,0*degree,cs),'unique') % define component with Euler angles
tilted_ND = symmetrise(orientation.byEuler(0*degree,30*degree,0*degree,cs),'unique') % define component with Euler angles
basal_TD = symmetrise(orientation.byEuler(0*degree,90*degree,0*degree,cs),'unique') % define component with Euler angles
basal_RD = symmetrise(orientation.byEuler(90*degree,90*degree,0*degree,cs),'unique') % define component with Euler angles

% Define a texture component for the cubic phase
rotated_cube = orientation.byMiller([0 0 1],[0 1 1],cs); % define component with directions
% Define a texture fibre for the cubic phase
gamma_fibre = fibre(Miller(1,1,1,cs),xvector);
alpha_fibre = fibre(Miller(1,1,0,cs),yvector);

for strip_index = 0:num_strips-1
    
    ebsd_cutmap = cutmap(strip_index); % read out ebsd_cutmap from the Map object
    
    % plot the IPF map, pole figures and odf slices
    outputFileName = strcat(analysis_path,sample_name,'_IPF_map_strip_',num2str(strip_index))
    IPF_map_plot(phase, ebsd_cutmap, outputFileName, visible)
    
    ori_strip = ebsd_cutmap('Ti-Hex').orientations
    outputFileName = strcat(analysis_path,sample_name,'_PF_strip_',num2str(strip_index))
    [maxval] = pole_figure_plot(phase, ori_strip, CS, contour_step, pf_max, outputFileName, visible);
    PF_basal_max(strip_index+1) = maxval(1);
    PF_prismatic1_max(strip_index+1) = maxval(2);
    PF_prismatic2_max(strip_index+1) = maxval(3);
    
    outputFileName = strcat(analysis_path,sample_name,'_ODF_strip_',num2str(strip_index))
    psi=calcKernel(ori_strip);
    HALF_WIDTH = psi
    odf_strip = calcDensity(ori_strip,'kernel',psi);
    ODF_plot(phase, odf_strip, odf_max, outputFileName, specSym, visible)
     
    % caclulate texture index
    TEXTURE_INDEX_strip(strip_index+1) = textureindex(odf_strip)
    
    % calculate strength of ODF maxima
    [odf_strip_max(strip_index+1),ori_strip_max(strip_index+1)]= max(odf_strip)
    
    % calculate misorientation of ODF maxima wrt 0002 in CD
    mori = (orientation.byEuler(0*degree,0*degree,0*degree,cs))*ori_strip_max
    misorientation_ODF_max = angle(mori)/degree

    % calculate PHI angle of ODF maxima
    phi1 = ori_strip_max.phi1
    PHI = ori_strip_max.Phi
    phi2 = ori_strip_max.phi2
    
    % seperate the texture components and calculate the volume fractions
    total_volume = length(ebsd_cutmap) % calculate the total volume as the number of points in the map

    % seperate a texture component and calculate the volume fraction
    ebsd_basal_ND = ebsd_cutmap('Ti-Hex').findByOrientation(basal_ND, misorientation*degree);
    basal_ND_volume = length(ebsd_basal_ND);
    basal_ND_volume_fraction(strip_index+1) = (basal_ND_volume/total_volume)
    
    % seperate a texture component and calculate the volume fraction
    ebsd_tilted_ND = ebsd_cutmap('Ti-Hex').findByOrientation(tilted_ND, misorientation*degree);
    tilted_ND_volume = length(ebsd_tilted_ND);
    tilted_ND_volume_fraction(strip_index+1) = (tilted_ND_volume/total_volume)
    
    % seperate a texture component and calculate the volume fraction
    ebsd_basal_TD = ebsd_cutmap('Ti-Hex').findByOrientation(basal_TD, misorientation*degree);
    basal_TD_volume = length(ebsd_basal_TD);
    basal_TD_volume_fraction(strip_index+1) = (basal_TD_volume/total_volume)
    
    % seperate a texture component and calculate the volume fraction
    ebsd_basal_RD = ebsd_cutmap('Ti-Hex').findByOrientation(basal_RD, misorientation*degree);
    basal_RD_volume = length(ebsd_basal_RD);
    basal_RD_volume_fraction(strip_index+1) = (basal_RD_volume/total_volume)
    
end

%% Open and write to file to save the different texture strength values

fileTS = fopen(fullfile(analysis_path, strcat(sample_name,'_texture_strength_',num2str(num_strips),'.txt')),'w');
fprintf(fileTS, 'Strip Index \t Texture Index \t ODF Max \t Misorientation of ODF Max \t phi1 Angle of ODF Max \t PHI Angle of ODF Max \t phi2 Angle of ODF Max \t {0002} PF Max \t {10-10} PF Max \t {11-20} PF Max \t Basal ND Volume Fraction\t Tilted ND Volume Fraction \t Basal TD Volume Fraction \t Basal RD Volume Fraction \n');

for strip_index = 0:num_strips-1    
    % write the texture strength values to file
    fprintf(fileTS, '%f\t%f\t%f\t%f\t%f\t%f\t%f%f\t%f\t%f\t%f\t%f\t%f\t%f\n', strip_index, TEXTURE_INDEX_strip(strip_index+1), odf_strip_max(strip_index+1), misorientation_ODF_max(strip_index+1), rad2deg(phi1(strip_index+1)), rad2deg(PHI(strip_index+1)), rad2deg(phi2(strip_index+1)), PF_basal_max(strip_index+1), PF_prismatic1_max(strip_index+1), PF_prismatic2_max(strip_index+1), basal_ND_volume_fraction(strip_index+1), tilted_ND_volume_fraction(strip_index+1), basal_TD_volume_fraction(strip_index+1), basal_RD_volume_fraction(strip_index+1))
end

% close any open files
fclose(fileTS);

%% Plot the texture variation

strip_index = 0:num_strips-1

TEXTURE_INDEX_figure = figure();
hold on
plot(strip_index+1,TEXTURE_INDEX_square,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('Texture Index')
hold off
legend('Texture Index')
saveas (TEXTURE_INDEX_figure, strcat(analysis_path,sample_name,'_TI_line_plot_',num2str(num_squares), '.png'));

odf_max_figure = figure();
hold on
plot(strip_index+1,odf_square_max,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('ODF Maximum')
hold off
legend('ODF Maximum')
saveas (odf_max_figure, strcat(analysis_path,sample_name,'_odf_max_line_plot_',num2str(num_squares), '.png'));

phi1_figure = figure();
hold on
plot(strip_index+1,rad2deg(phi1),'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('Phi1 Angle')
hold off
legend('Phi1 Angle')
saveas (phi1_figure, strcat(analysis_path,sample_name,'_phi1_line_plot_',num2str(num_squares), '.png'));

PHI_figure = figure();
hold on
plot(strip_index+1,rad2deg(PHI),'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('PHI Angle')
hold off
legend('PHI Angle')
saveas (PHI_figure, strcat(analysis_path,sample_name,'_PHI_line_plot_',num2str(num_squares), '.png'));

phi2_figure = figure();
hold on
plot(strip_index+1,rad2deg(phi2),'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('Phi2 Angle')
hold off
legend('Phi2 Angle')
saveas (phi2_figure, strcat(analysis_path,sample_name,'_phi2_line_plot_',num2str(num_squares), '.png'));

PF_basal_figure = figure();
hold on
plot(strip_index+1,PF_basal_max,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('{0002} Pole Figure Max.')
hold off
legend('{0002} Pole Figure Max.')
saveas (PF_basal_figure, strcat(analysis_path,sample_name,'_PF_basal_line_plot_',num2str(num_squares), '.png'));

PF_prismatic1_figure = figure();
hold on
plot(strip_index+1,PF_prismatic1_max,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('{10-10} Pole Figure Max.')
hold off
legend('{10-10} Pole Figure Max.')
saveas (PF_prismatic1_figure, strcat(analysis_path,sample_name,'_PF_prismatic1_line_plot_',num2str(num_squares), '.png'));

PF_prismatic2_figure = figure();
hold on
plot(strip_index+1,PF_prismatic2_max,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('{11-20} Pole Figure Max.')
hold off
legend('{11-20} Pole Figure Max.')
saveas (PF_prismatic2_figure, strcat(analysis_path,sample_name,'_PF_prismatic2_line_plot_',num2str(num_squares), '.png'));