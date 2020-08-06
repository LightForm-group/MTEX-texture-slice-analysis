%% Specify Crystal and Specimen Symmetries

% crystal symmetry
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
ctf_file = '/SampleD5_repeat Data.ctf'
data_path = strcat(mech_tester,data_file,sample_name,ctf_file)

analysis_folder = '/Example Analysis/'
analysis_path = strcat(mech_tester,analysis_folder,sample_name,'/') % path for saving the data

% which files to be imported
fname = [pname data_path]

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame')

%% Rotate the data

rot = rotation('Euler', 90*degree, 90*degree, 0*degree);

ebsd = rotate(ebsd,rot,'keepXY'); % rotate the orientation data
ebsd = rotate(ebsd,90*degree,'keepEuler') % rotate the spatial data

% ebsd = rotate(ebsd,rot); % or we could just use this to rotate the map as well

fprintf('Note, the (x,y) origin on the map will have changed and x or y could be negative!')

%% Plot the IPF colour map

phase = 'alpha'
visible = 'off'
outputFileName = strcat(analysis_path,sample_name,'_IPF_map_entire_region')
IPF_map_plot(phase, ebsd, outputFileName, visible)

%% Plot the pole figures for the whole compression sample

phase = 'alpha'
ori = ebsd('Ti-Hex').orientations
contour_step = 0.1
pf_max = 2.5
outputFileName = strcat(analysis_path,sample_name,'_PF_entire_region')
visible = 'off'
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
x_top = -250
y_left = 2400
x_bottom = -5100
y_right = 5600

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

num_strips = 20; % number of strips to cut the map into (resolution)

% define the size of the EBSD map
ebsd_grid = ebsd.gridify;
ebsd_shape = size(ebsd_grid.id);
original_y = ebsd_shape(1);
original_x = ebsd_shape(2);
stepSize = ebsd_grid.dx;

x_min = (sqrt(x_origin * x_origin)/stepSize);
x_max = original_x + (sqrt(x_origin * x_origin)/stepSize);
x_length = x_max - x_min;

y_min = (sqrt(y_origin * y_origin)/stepSize);
y_max = original_y + (sqrt(y_origin * y_origin)/stepSize);
y_length = y_max - y_min;

% used if splitting into strips along y
y_width = floor(y_length / num_strips); % round to nearest integer
y_axis = (1:num_strips);

% used if splitting into strips along x
x_width = floor(x_length / num_strips); % round to nearest integer
x_axis = (1:num_strips);

cutmap = containers.Map('KeyType', 'int32', 'ValueType', 'any'); % creates an empty Map object

for strip_index = 0:num_strips-1
    % separate the map section
    
    % set out the coordinates for the edge of the region
    % note, region is defined as x,y origin and an x,y width which is added onto the origin
    
    % if splitting into strips along y (breaking up y)
    % y_min_strip = strip_index * y_width;
    % region = [x_min*stepSize, y_min_strip*stepSize, x_length*stepSize, y_width*stepSize];
    
    % if splitting into strips along y (breaking up y) and x is negative
    % y_min_strip = strip_index * y_width;
    % region = [-x_min*stepSize, y_min_strip*stepSize, -x_length*stepSize, y_width*stepSize];
    
    % if splitting into strips along x (breaking up x)
    % x_min_strip = strip_index * x_width;
    % region = [x_min_strip*stepSize, y_min*stepSize, x_width*stepSize, y_length*stepSize];
    
    % if splitting into strips along x (breaking up x) and x is negative
    x_min_strip = strip_index * x_width + x_min;
    region = [-x_min_strip*stepSize, y_min*stepSize, -x_width*stepSize, y_length*stepSize];
    
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
% basal_ND = symmetrise(orientation.byMiller([0 0 0 1],[1 0 -1 0],cs),'unique'); % define component with directions
basal_ND = symmetrise(orientation.byEuler(0*degree,0*degree,0*degree,cs),'unique') % define component with Euler angles
tilted_ND = symmetrise(orientation.byEuler(0*degree,30*degree,0*degree,cs),'unique') % define component with Euler angles

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
    pole_figure_plot(phase, ori_strip, CS, contour_step, pf_max, outputFileName, visible);
    
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
    PHI=ori_strip_max.Phi
    
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
    
end

%% Open and write to file to save the different texture strength values

fileTS = fopen(fullfile(analysis_path, strcat(sample_name,'_texture_strength_',num2str(num_strips),'.txt')),'w');
fprintf(fileTS, 'Strip Index\tTexture Index\tODF Max\tPHI Angle of ODF Max.\tMisorientation of ODF Max.\tBasal ND Volume Fraction\tTilted ND Volume Fraction\n');

for strip_index = 0:num_strips-1    
    % write the texture strength values to file
    fprintf(fileTS, '%f\t%f\t%f\t%f\t%f\t%f\t%f\n', strip_index, TEXTURE_INDEX_strip(strip_index+1), odf_strip_max(strip_index+1), PHI(strip_index+1), misorientation_ODF_max(strip_index+1), basal_ND_volume_fraction(strip_index+1), tilted_ND_volume_fraction(strip_index+1))
end

% close any open files
fclose(fileTS);

%% Plot the texture variation and the FE results

strip_index = 1:num_strips

vol_frac_line_figure = figure();
hold on
plot(strip_index,basal_ND_volume_fraction*100,'Color',[1,0,0],'lineWidth',2) % red);
xlabel('Slice Number')
ylabel('Volume Fraction (%)')
hold on
plot(strip_index,tilted_ND_volume_fraction*100,'Color',[0,1,0],'lineWidth',2) % green);
hold off
legend('Basal ND','Tilted ND')
saveas (vol_frac_line_figure, strcat(analysis_path,sample_name,'_vol_frac_line_plot_',num2str(num_strips), '.bmp'));
 
text_ind_line_figure = figure();
hold on
plot(strip_index,TEXTURE_INDEX_strip,'Color',[1,0,0],'lineWidth',2) % red);
xlabel('Slice Number')
ylabel('Texture Index or ODF Max')
hold on
plot(strip_index,odf_strip_max,'Color',[0,1,0],'lineWidth',2) % green);
hold off
legend('Texture Index','ODF Maximum')
saveas (text_ind_line_figure, strcat(analysis_path,sample_name,'_texture_index_line_plot_',num2str(num_strips), '.bmp'));

text_ODF_ori_figure = figure();
hold on
plot(strip_index, PHI,'Color',[0,0,1],'lineWidth',2); % blue
xlabel('Slice Number')
ylabel('PHI Angle or Misorientation')
hold on
plot(strip_index, misorientation_ODF_max,'Color',[0,1,0],'lineWidth',2); % green
legend('PHI Angle of ODF Max.','Misorientation of ODF Max.')
hold off
saveas (text_ODF_ori_figure, strcat(analysis_path,sample_name,'_texture_ODF_orientation_plot_',num2str(num_strips), '.bmp'));