%%%%%%%%%%%%%% COASTCAMS %%%%%%%%%%%%%
% Please cite when using CoastCams: 

% Nuyts, S., Almar, R., Morichon, D., Dealbera, S., Abalia, A., Muñoz, J. M., Abessolo, G. O., & Regard, V. (2023). 
% CoastCams: A MATLAB toolbox making accessible estimations of nearshore processes, mean water levels, and morphology from timestack images. 
% Environmental Modelling & Software, 168, 105800. https://doi.org/https://doi.org/10.1016/j.envsoft.2023.105800 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script will compute nearshore processes from timestack images (e.g. wave
% height, wave period, wave celerity), mean water levels, 
% and morphology (e.g. shoreline positions).

%% A: Housekeeping
close all
clearvars 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START ESSENTIAL USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% B: Set up paths
% Check if CoastCams directory exists, if not, clone it
repo_path = fullfile(pwd, 'CoastCams');
if ~exist(repo_path, 'dir')
    !git clone https://github.com/NuytsSiegmund/CoastCams.git
else
    fprintf('CoastCams directory already exists. Using existing directory.\n');
end

% Set paths
addpath(genpath(repo_path));

% B1: Path to User Scripts (now in the repo)
user_scripts_path = fullfile(repo_path, 'UsersScripts');
addpath(genpath(user_scripts_path));

% B2: Path to Timestack Images (you may need to adjust this)
img_path = fullfile(repo_path, 'Timestacks');

% B3: Output path: Location where outputs will be stored
out_path = fullfile(repo_path, 'Output');
if ~exist(out_path, 'dir')
    mkdir(out_path);
end

%% C: Select Timestack Images
Img = dir(fullfile(img_path, 'S_1_*.jpeg')); % Name should include datetime as S_1_yyyymmddHHMM e.g. "S_1_202212310700.jpeg"

%% D: Parameters that require user input
% D1: Parameters for image processing
dt          = 1/2; % dt = 1/freq with freq = 2 (frequency acquisition of the camera: 2 images per second)
H_camera    = 27.240; % Camera height above MSL [m] 
res         = 0.1; % Resolution: Size of each pixel [m]
rotation    = 270; % Waves in timestack image should come from top-left corner - rotate accordingly 

% D2: Reading camera parameters 
CoordCam    = [0,0,H_camera]; 
dur         = 14; % Duration of the Timestack in minutes

% D3: Parameters for cross-correlation computation
Nlim        = 1600; % Width of Timestack Image calculation (can decrease to only calculate portion of transcet and to speed-up calculations)
dpha        = 1; % Time lag parametes used for cross-correlation calculations (should be smaller than smallest wave period)
icmin       = 1; % Minimum value on the x-axis for pre-processing
icmax       = 680; % Maximum value on the x-axis for pre-processing
dc          = 100; % see Thuan et al., 2019 and Abessolo et al., 2020 (Figure 2 above)
resc        = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START OPTIONAL USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% D4: Selection of shoreline method
ShoreMethod = 1; % 1 = grayscale (specifically designed for Socoa rocky shore);
                 % 2 = Red- minus Blue- channel method;
                 % 3 = Colour convergence method;
threshold = 30;  % threshold between peaks and troughs (only used in ShoreMethod = 1)
                 
% D5: Plot option: If user wants to plot figures for each individual timestack
plotoption_shore = 1; % plots for shoreline position
plotoption_imfor = 1; % plots for image formating
plotoption_crosscor = 1; % plots for cross-correlation matrix
compare_bathymetry_option = true;

% D6: Camera number: Only used for naming output files
camera = 1;

%% E: Variables to be calculated from Timestack images
Img_date = []; WaveCelerity = []; df = []; WaveLength = []; Hs_TS = []; Tp_TS = []; Tm_TS = []; WaveEnergy = []; RollerLength = [];  BreakpointDepth = [];
BreakpointLocation = []; Cf1 = []; WLe1 = []; Depth_S = []; ShorePosition = []; Stack_av = []; BL_Coordinates = []; WaterDepth_L = []; Time_TS = [];

%% F: Coordinate System
% (Not an essential step)
% Transfer pixel coordinates to real life coordinates 
% Create a .txt file with real life xy coordinates for pixel coordinates
coordinate_option = 0; % 0 if no tranformation is needed, 1 if user wants to transform pixel coordinates to real life coordinates
coordinates = []; % Initialize as empty array
if coordinate_option == 1
    coordinate_file = fullfile(user_scripts_path, 'Coordinates.txt');
    if exist(coordinate_file, 'file')
        coordinates = load(coordinate_file);
    else
        warning('Coordinates.txt file not found. Proceeding without coordinate transformation.');
        coordinate_option = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END ALL USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% G: Computing Parameters from Timestack Images 
for i = 1:length(Img)
    % G1: Extract datetime from Timestack Images 
    Img_date = datenum(Img(i).name(5:end-5), 'yyyymmddHHMM');
    Time_TS = [Time_TS; datetime(Img_date, 'ConvertFrom', 'datenum')];
    
    % G2: Load and read Timestack Images
    fprintf('Analysing Timestack from %s\n', datestr(Img_date, 'mmmm dd, yyyy HH:MM AM'));
    try
        Timestack = imread(fullfile(img_path, Img(i).name));
    catch
        fprintf("Error reading image file %s\n", Img(i).name);
        continue
    end
    
    Timestack_rot = imrotate (Timestack, rotation); % Rotate Timestack Image so that waves arrive from top left corner     
    [nt,nc,ncol] = size(Timestack_rot); % Get dimensions from Timestack Image   

    % G3: Extract Shoreline Position
    try
        [shoreline_positions] = shoreline_position(Timestack_rot, threshold, res, ShoreMethod, plotoption_shore);
        ShorePosition = [ShorePosition; shoreline_positions]; 
    catch ME
        warning('Error in shoreline position extraction: %s', ME.message);
        ShorePosition = [ShorePosition; NaN];
    end
      
    % G4: Image Formatting   
    dx = 0:0.1:nc/10;
    dxXX  = (dx(1:nc)).';
    sLimit = [1:1:icmax - dc]';
    
    if Nlim > nt
        error('Error: Nlim must be smaller than the time dimension of the Timestack image. Please decrease Nlim');
    end
    if icmax > nc
        error('Error: icmax must be smaller than the cross-shore dimension of the Timestack image. Please decrease icmax');
    end
     
    if plotoption_imfor == 1
        try
            figure(1)
            set(gcf, 'Position', [10 10 600 1200], 'Color', [1, 1, 1])
            subplot(311)
            image(Timestack_rot)
            xlabel('X[Pixels]')
            ylabel('Frames')
            title(sprintf('%s',datestr(datenum(Img(i).name(5:16),'yyyymmddHHMM'),'mmmm dd, yyyy HH:MM AM')))
            set(gca, 'fontsize',14)
                               
            subplot(312)
            image(Timestack_rot)
            xlabel('X [m]')
            ylabel('Frames')
            set(gca, 'fontsize', 14)
            set(gca,'XTick',0:200:nc)
            set(gca,'XTickLabel',0:20:nc/10)   
        catch ME
            warning('Error in plotting image formatting: %s', ME.message);
        end
    end

    So = double(Timestack_rot);
    S_std = nanstd(So);
    
    S_std0 = smooth(S_std(:,:,3)); % Only the blue band is considered 
    iS = find(S_std0>0.5*max(S_std0)); 
    
    S0 = double(So(:,:,3)); 
    S1 = S0(1:min([Nlim size(S0,1)]),:);
    
     % G5: Pre-processing and check
    try
        [S2] = ImagePreProcessing_CoastCams(S1,icmin,icmax,dt,resc,1); 
    catch ME
        warning('Error in image pre-processing: %s', ME.message);
        S2 = NaN(size(S1));
    end

    if plotoption_imfor == 1
        try
            subplot(313)
            image(S0)
            box on
            hold on
            plot(smooth(S_std(:,:,1),30).*40,'r')
            plot(smooth(S_std(:,:,2),30).*40,'g')
            plot(smooth(S_std(:,:,3),30).*40,'b')
            xlim([1 size(S0,2)])
            patch([min(iS) min(iS)],[0 size(S0,1)],'k','Linewidth',2)
            patch([max(iS) max(iS)],[0 size(S0,1)],'k','Linewidth',2)
            hold off
            box on
            xlabel('X [Pixels]')
            ylabel('Frames')
            set(gca, 'fontsize',14)
        catch ME
            warning('Error in plotting image formatting: %s', ME.message);
        end
    end
    
    % G6: Compute Wave Parameters
    try 
        [C,depth_s,depth_l,hs,hm,Tm,Tp,Er,rolL,nbwave,Breakstd,Breakmean1,Breakmean2,BreakLocs,BreakDepth] = WaveParameters_CoastCams(fullfile(img_path, Img(i).name),dt,CoordCam,dxXX./100,dur);
    catch ME
        warning('Error in WaveParameters_CoastCams: %s', ME.message);
        C = NaN; depth_s = NaN; depth_l = NaN; hs = NaN; hm = NaN; Tm = NaN; Tp = NaN; Er = NaN; rolL = NaN; nbwave = NaN;
        Breakstd = NaN; Breakmean1 = NaN; Breakmean2 = NaN; BreakLocs = NaN; BreakDepth = NaN;
    end
   
    if ~isempty(depth_s) && ~all(isnan(depth_s))
        WaterDepth(i, 1:numel(depth_s)) = movmean(depth_s./10, 10);
    else
        WaterDepth(i, :) = NaN;
    end
    
    Hs_TS = [Hs_TS; hs];
    Tp_TS = [Tp_TS; Tp];
    Tm_TS = [Tm_TS; Tm];
    WaveEnergy(i,1:size(Er,2)) = nanmean(Er);
    RollerLength(i,1:size(rolL,2)) = nanmean(rolL);
    BreakDepthf(i, 1:size(depth_s,2)) = nanmean(depth_s); 
    BreakpointLocation = [BreakpointLocation; (BreakLocs+(dc/2))];  
    BreakpointDepth = [BreakpointDepth; BreakDepth];
        
    % G7: Cross-correlation Calculations   
    try
        [Cf1, WLe1, Tp1, Hs1, RM] = CrossCorrelation_CoastCams(S2, dpha, dt, dc);
        disp('CrossCorrelation_Coastcams executed successfully');      
    catch ME
        warning('Error in CrossCorrelation_Coastcams: %s', ME.message);
        Cf1 = NaN(size(S2, 2), 1);
        WLe1 = NaN(size(S2, 2), 1);
    end

    WaveCelerity(i,1:size(Cf1,2)) = movmean(Cf1./10, 10); %get moving average
    WaveLength(i,1:size(WLe1,2)) = movmean(WLe1,10);

    % G8: Calculate depth using linear wave theory
    try 
        [df] = LinearC(Tp, WaveCelerity(i,:), 0.01);
    catch ME
        warning('Error in LinearC: %s', ME.message);
        df = NaN(1,numel(sLimit));
    end
    WaterDepth_L = [WaterDepth_L; movmean(df(:,1:numel(sLimit)), 10)];   

    % G9 : Calculating Additional Parameters
    % Sea level Anomaly
    try
        [rows, cols] = size(WaveCelerity);
        Nr = min(rows, 1);  % Use 1 as the minimum smoothing window for rows
        Nc = min(cols, 30); % Use the smaller of 30 or the number of columns
        Csmooth_S = smooth2(WaveCelerity, Nr, Nc); % smooth signal
        [rows, cols] = size(WaterDepth_L);
        Nr = min(rows, 1);  % Use 1 as the minimum smoothing window for rows
        Nc = min(cols, 30); % Use the smaller of 30 or the number of columns
        Csmooth_L = smooth2(WaterDepth_L, Nr, Nc);
        sLimit_size = min([size(Csmooth_S, 2), size(Csmooth_L, 2), numel(sLimit)]);
        SLA_S = Csmooth_S(:,1:sLimit_size) - nanmean(Csmooth_S(:,1:sLimit_size)); % calculate sea level anomaly using shallow water
        SLA_L = Csmooth_L(:,1:sLimit_size) - nanmean(Csmooth_L(:,1:sLimit_size)); % calculate sea level anomaly using linear wave
    catch ME
        warning('Error in calculating sea level anomaly: %s', ME.message);
        SLA_S = NaN(size(WaveCelerity, 1), numel(sLimit));
        SLA_L = NaN(size(WaterDepth_L, 1), numel(sLimit));
    end
               
    % Water Level and Relative Tidal Range
    Level_TS = nanmean(WaveCelerity, 2);
    nRTR_thresh = Level_TS - min(Level_TS);
    nRTR_thresh(nRTR_thresh < 0.2) = NaN; % threshold to avoid dividing by 0
    RTR = hs./nRTR_thresh;  
    
    % G10: Produce an average Timestack Image to see change over time
    Stack_av = [Stack_av; mean(Timestack_rot, 1, 'omitnan')];
        
    % G11: Convert to real life coordinates (if option is enabled)
    if coordinate_option == 1
        try
            for j = 1:length(BreakLocs)
                if BreakLocs(j) <= size(coordinates, 1)
                    BL_Lam_Mid = coordinates(BreakLocs(j), :);
                    BL_Coordinates = [BL_Coordinates; BL_Lam_Mid];
                else
                    warning('BreakLocs exceeds the size of coordinates. Skipping coordinate conversion for iteration %d.', j);
                    BL_Coordinates = [BL_Coordinates; NaN(1, size(coordinates, 2))];
                end
            end
        catch ME
            warning('Error in coordinate conversion: %s', ME.message);
            disp(getReport(ME, 'extended'));
        end
   end
end

%% H: Pixels coordinates to real life coordinates
if coordinate_option == 1
    BL_Coordinates_x = BL_Coordinates(:, 1);
    BL_Coordinates_y = BL_Coordinates(:, 2);
    BL_Coordinates = timetable(Time_TS, BL_Coordinates_x, BL_Coordinates_y);
    writetimetable(BL_Coordinates, fullfile(out_path, 'BL_Coordinates.txt'));
end

%% I: Bathymetry Calculation 
% Debug output
disp(['Size of WaterDepth: ', mat2str(size(WaterDepth))]);
disp(['Size of SLA_S: ', mat2str(size(SLA_S))]);

if ~isempty(WaterDepth) && ~isempty(SLA_S)
    % Transpose SLA_S to match WaterDepth orientation
    SLA_S_transposed = SLA_S';
    disp(['Size of SLA_S after transposition: ', mat2str(size(SLA_S_transposed))]);
    
    % Determine the common size
    [rows_wd, cols_wd] = size(WaterDepth);
    [rows_sla, cols_sla] = size(SLA_S_transposed);
    common_rows = min(rows_wd, rows_sla);
    common_cols = min(cols_wd, cols_sla);
    
    % Adjust both WaterDepth and SLA_S
    WaterDepth_adjusted = WaterDepth(1:common_rows, 1:common_cols);
    SLA_S_adjusted = SLA_S_transposed(1:common_rows, 1:common_cols);
    
    disp(['Size of WaterDepth_adjusted: ', mat2str(size(WaterDepth_adjusted))]);
    disp(['Size of SLA_S_adjusted: ', mat2str(size(SLA_S_adjusted))]);
    
    if ~isequal(size(WaterDepth_adjusted), size(SLA_S_adjusted))
        error('Failed to adjust WaterDepth and SLA_S to compatible sizes.');
    end
    
    % Calculate bathymetry for each timestep
    Bathymetry_full = WaterDepth_adjusted - SLA_S_adjusted;
    
    % Calculate mean bathymetry across all timesteps
    Bathymetry = nanmean(Bathymetry_full, 1);
    
    % Debug output
    disp(['Size of Bathymetry_full: ', mat2str(size(Bathymetry_full))]);
    disp(['Size of Bathymetry (mean): ', mat2str(size(Bathymetry))]);
    disp(['Range of Bathymetry values: ', num2str(min(Bathymetry, [], 'omitnan')), ' to ', num2str(max(Bathymetry, [], 'omitnan'))]);
    
    % Create cross-shore distance vector
    x = (0:length(Bathymetry)-1) * res;
    disp(['Cross-shore distance range: 0 to ', num2str(max(x)), ' meters']);
    disp(['Resolution (res): ', num2str(res), ' meters']);
    
    % Plot the bathymetry profile
    figure;
    plot(x, Bathymetry, 'b-', 'LineWidth', 2);
    xlabel('Cross-shore distance [m]', 'FontSize', 14);
    ylabel('Depth [m]', 'FontSize', 14);
    title('Mean Bathymetry Profile', 'FontSize', 16);
    grid on;
    set(gca, 'YDir', 'reverse');  % Invert y-axis to show depth properly
    set(gcf, 'Color', 'w');
else
    warning('WaterDepth or SLA_S is empty. Unable to calculate Bathymetry.');
    Bathymetry = [];
    x = [];
end

% Save Bathymetry to output structure or variable
if exist('Output', 'var')
    Output.Bathymetry = Bathymetry;
    Output.BathymetryDistance = x;
else
    % If 'Output' doesn't exist, create a new variable
    BathymetryResult.profile = Bathymetry;
    BathymetryResult.x = x;
end

% Optionally, save the bathymetry profile to a file
save(fullfile(out_path, 'bathymetry_profile.mat'), 'Bathymetry', 'x');

%% J: Collecting all parameters from Timestack Image
num_rows = length(Time_TS);
variables = {BreakpointLocation, BreakpointDepth, Hs_TS, WaveEnergy, RollerLength, WaveCelerity, ...
             Tp_TS, Tm_TS, WaveLength, WaterDepth_L, ShorePosition, SLA_S, SLA_L, RTR, Bathymetry};

% Ensure all variables have the same number of rows
for i = 1:length(variables)
    if size(variables{i}, 1) < num_rows
        variables{i} = padarray(variables{i}, [num_rows - size(variables{i}, 1), 0], NaN, 'post');
    elseif size(variables{i}, 1) > num_rows
        variables{i} = variables{i}(1:num_rows, :);
    end
end

% Create timetable
WP = timetable(Time_TS, variables{:});
WP.Properties.VariableNames = {'BreakPointLocation', 'BreakpointDepth', 'Hs', 'WaveEnergy', 'RollerLength', ...
                               'WaveCelerity', 'Tp', 'Tm', 'WaveLength', 'WaterDepth', 'ShorelinePosition', ...
                               'SLA_S', 'SLA_L', 'RTR', 'Bathymetry'};

% Resample to regular 15-minute intervals
Data_TS = retime(WP, 'regular', 'mean', 'TimeStep', minutes(15));

% Show average daily values
WP_daily = retime(Data_TS, 'daily', 'mean');
disp(WP_daily)

%% K: Save Outputs
% K1: Save Workspace
date = datetime('now', 'Format', 'yyyyMMdd');
Ft = char(date);
Fn = sprintf('Camera_%d_', camera);
FileName = ['Output_S01_', Fn, Ft];
save(fullfile(out_path, FileName));

% K2: Export to .txt
txt = fullfile(out_path, [FileName '.txt']);
writetimetable(WP, txt);

%% L: Error Handling and Logging
logFile = fullfile(out_path, 'error_log.txt');
if ~exist(logFile, 'file')
    fid = fopen(logFile, 'w');
    fprintf(fid, 'CoastCams Error Log\n\n');
    fclose(fid);
end

% Function to log errors
function logError(errorMsg)
    global logFile
    fid = fopen(logFile, 'a');
    fprintf(fid, '%s: %s\n', datestr(now), errorMsg);
    fclose(fid);
end

%% M: Plots
plot_coastcams_main(Time_TS, Stack_av, SLA_S, Hs_TS, Tp_TS, rotation)

%% N: Final Output and Cleanup
% Save all figures
figHandles = findall(0, 'Type', 'figure');
for i = 1:length(figHandles)
    saveas(figHandles(i), fullfile(out_path, sprintf('Figure_%d.png', i)));
end

% Generate a summary report
summaryFile = fullfile(out_path, 'analysis_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'CoastCams Analysis Summary\n\n');
fprintf(fid, 'Analysis Date: %s\n', datestr(now));
fprintf(fid, 'Number of Images Processed: %d\n', length(Img));
fprintf(fid, 'Time Range: %s to %s\n', datestr(Time_TS(1)), datestr(Time_TS(end)));
fprintf(fid, 'Mean Significant Wave Height: %.2f m\n', mean(Hs_TS, 'omitnan'));
fprintf(fid, 'Mean Peak Wave Period: %.2f s\n', mean(Tp_TS, 'omitnan'));
fclose(fid);

% Display completion message
disp('CoastCams analysis completed successfully.');
disp(['Results saved in: ' out_path]);
disp(['Summary report: ' summaryFile]);

%% End of script

