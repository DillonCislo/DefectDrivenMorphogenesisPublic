%% Order Parameter Correlation Validation =================================
% This is a script to validate the significance of the orientational order
% observed in a growing Parhyale embryo relative to synthetically generated
% totally disordered data sets. It will produced the analysis displayed in
% Fig. S9
%
% by Dillon Cislo
%==========================================================================

% Add relevant file structure to MATLAB path
[DDMDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
DDMDir = [DDMDir '/..'];
addpath(genpath(DDMDir));
clear DDMDir

%% Set Up Pipeline ========================================================
clear; close all; clc;

% Load the appropriate geometric characteristics of the embryo featured in
% 'Dataset_20180603' at the reference time points.
% NOTE: 'tissueLengths' and 'tissueWidths' are in units of the average cell
% size
load('Correlation_Function_Validation_20210303.mat', 'tissueLengths');
load('Correlation_Function_Validation_20210303.mat', 'tissueWidths');
load('Correlation_Function_Validation_20210303.mat', 'cellLengths');

%--------------------------------------------------------------------------
% Time Point Validation Choices
%--------------------------------------------------------------------------

% The ID of the time point of interest
tidx = 15;
% tidx = 161;
% tidx = 234;

% A rough indication of the average length scale of the cells
cellSize = cellLengths(tidx);

% The measured length of the included region
measL = cellSize * tissueLengths(tidx);

% The measured width of the included region
measW = cellSize * tissueWidths(tidx);

clear cellLengths tissueLengths tissueWidths

%--------------------------------------------------------------------------
% Point Set Generation Options
%--------------------------------------------------------------------------

% Points are generated inside a large box. Only points within a smaller
% box, itself contained within the large box, will be included in
% calculations of the correlation functions
incBoxX = measW .* [-1 1] / 2;
incBoxY = measL .* [-1 1] / 2;
bigBoxX = incBoxX * 2;
bigBoxY = incBoxY * 2;

% Lattice Options Point Generation Instructions:---------------------------
% 'square' - number set by (xLim, yLim, sideLength)
% 'noisySquare' - number set by (xLim, yLim, sideLength)
% 'rectangular' - number set by (xLim, yLim, numRows, numCols)
% 'noisyRectangular' - number set by (xLim, yLim, numRows, numCols)
% 'equilateral' - number set by (xLim, yLim, sideLength)
% 'noisyEquilateral' - number set by (xLim, yLim, sideLength)
% 'triangular' - number set by (xLim, yLim, numRows, numCols)
% 'noisyTriangular' - number set by (xLim, yLim, numRows, numCols)
% 'random' - number set by (numPoints)
% 'randomMinDist' - number set by (xLim, yLim, minDist, numPoints)

latticeOptions = struct();
latticeOptions.xLim = bigBoxX;
latticeOptions.yLim = bigBoxY;

% Random Point Locations --------------------------------------------------
latticeOptions.pointSetType = 'randomMinDist';
latticeOptions.minDist = 0.8 * cellSize;
latticeOptions.maxRandSamples = 30;
latticeOptions.numPoints = Inf;

%% Generate Correlation Functions =========================================
% This section generates the data contained in the following files:
%
%   'Random_Correlations_T015.mat' 'Random_Correlations_T161.mat'
%   'Random_Correlations_T234.mat'
%
% found on the associated Dryad data repository. This can take quite some
% time. If you download the data from the Dryad repository you can just
% skip to the figure making sections.
clc;

% The number of random ensemble realizations to generate
numSamples = 1000;

allXY = cell(numSamples, 1);
allIncNodeIDx = cell(numSamples, 1);

allPsi4 = cell(numSamples, 1);
allPsi6 = cell(numSamples, 1);

allPsi4Corr = cell(numSamples, 1);
allPsi6Corr = cell(numSamples, 1);
allNumCellsR = cell(numSamples, 1);

allGRMax = cell(numSamples, 1);
allGRX = cell(numSamples, 1);
allGRY = cell(numSamples, 1);
allGRR = cell(numSamples, 1);

normType = 'none';
useEffectiveVolume = true;
effectiveVolumeMethod = 'polygon';

RMax = 50;
numRBins = 200;
normDists = true;
r = linspace(0, RMax, numRBins);

% For reproducible random numbers
rng(1, 'twister');

for i = 1:numSamples
    
    fprintf('NOW PROCESSING ENSEMBLE #%d\n', i);
    
    % Generate the point set
    XY = CELTIGS.generate_points( latticeOptions );
    
    % Determine which nodes to include in the calculation
    incNodeIDx = (incBoxX(1) <= XY(:,1)) & (XY(:,1) <= incBoxX(2)) ...
        & (incBoxY(1) <= XY(:,2)) & (XY(:,2) <= incBoxY(2));
    incNodeIDx = find(incNodeIDx);
    
    % Calcualte the bond orientational order parameters
    psi4 = calculateBondOrientationalOrder(XY, 4, 2, false, false);
    psi6 = calculateBondOrientationalOrder(XY, 6, 2, false, false);
    
    % Calculate the bond orientational correlation functions
    [psi4Corr, ~, numCellsR, ~] = ...
        calculateBondOrientationalCorrelations( XY, psi4, ...
        'MaxR', RMax, 'NumRBins', numRBins, ...
        'IncludedNodes', incNodeIDx, ...
        'NormalizeDistances', normDists, 'PlotResults', false );
    
    [psi6Corr, ~, ~, ~] = ...
        calculateBondOrientationalCorrelations( XY, psi6, ...
        'MaxR', RMax, 'NumRBins', numRBins, ...
        'IncludedNodes', incNodeIDx, ...
        'NormalizeDistances', normDists, 'PlotResults', false );
    
    [GRMax, GRX, GRY, GRR, ~] = calculatePositionalCorrelations( ...
        XY, 'NumRBins', numRBins, 'EffectiveVolume', useEffectiveVolume, ...
        'EffectiveVolumeMethod', effectiveVolumeMethod, ...
        'IncludedNodes', incNodeIDx, 'Normalization', normType, ...
        'NormalizeDistances', normDists, 'PlotResults', false, ...
        'MaxR', RMax );
    
    % Update ensemble arrays
    allXY{i} = XY;
    allIncNodeIDx{i} = incNodeIDx;
    allPsi4{i} = psi4;
    allPsi6{i} = psi6;
    allPsi4Corr{i} = psi4Corr;
    allPsi6Corr{i} = psi6Corr;
    allNumCellsR{i} = numCellsR;
    allGRMax{i} = GRMax;
    allGRX{i} = GRX;
    allGRY{i} = GRY;
    allGRR{i} = GRR;
    
end
 
save( sprintf('Random_Correlations_T%03d.mat', tidx) );

%% Generate 4-Fold Orientational Order Comparisons ========================
clear; close all; clc;

% Load Synthetic Results for Early Times ----------------------------------
tidx = 15;

% Radial bins should be the same for all time points
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'r');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'RMax');

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measL');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measW');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'cellSize');

measL_T15 = measL;
measW_T15 = measW;
cellSize_T15 = cellSize;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allPsi4Corr');

psi4Corr_T15 = mean(cell2mat(allPsi4Corr), 1);

clear allPsi4Corr measL measW cellSize

% Load Synthetic Results for Intermediate Times ---------------------------
tidx = 161;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measL');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measW');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'cellSize');

measL_T161 = measL;
measW_T161 = measW;
cellSize_T161 = cellSize;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allPsi4Corr');

psi4Corr_T161 = mean(cell2mat(allPsi4Corr), 1);

clear allPsi4Corr measL measW cellSize

% Load Synthetic Results for Late Times -----------------------------------
tidx = 234;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measL');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measW');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'cellSize');

measL_T234 = measL;
measW_T234 = measW;
cellSize_T234 = cellSize;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allPsi4Corr');

psi4Corr_T234 = mean(cell2mat(allPsi4Corr), 1);

clear allPsi4Corr measL measW cellSize

%--------------------------------------------------------------------------
% Generate Figures
%--------------------------------------------------------------------------

fig = figure('Color', [1 1 1], 'visible', 'on');

% Axis bounds
xLim = [0.5 RMax];
yLim = [1e-4, 0.5];

hold on

% Plot correlation functions
plot(r, psi4Corr_T15, 'LineWidth', 1.5);
plot(r, psi4Corr_T161, 'LineWidth', 1.5);
plot(r, psi4Corr_T234, 'LineWidth', 1.5);

% plot(r, abs(psi4Corr_T15), 'LineWidth', 1.5);
% plot(r, abs(psi4Corr_T161), 'LineWidth', 1.5);
% plot(r, abs(psi4Corr_T234), 'LineWidth', 1.5);

% Plot tissue widths
plot( ...
    mean([measW_T161/cellSize_T161, measW_T234/cellSize_T234]) .* [1 1], ...
    [yLim(1) 0.8*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );

% Plot algebraic decay guide line
plot(r, 0.1 .* r.^(-1/4), '--k', 'LineWidth', 1.5);

hold off

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

xlim(xLim);
ylim(yLim);

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

xticks([1 10]);

% ylabel('\langle\psi_4({\itr/L})\psi_4(0)\rangle', ...
%     'FontWeight', 'bold', 'FontSize', 7);
% xlabel('{\it r/L} (Mean cell lengths)', ...
%     'FontWeight', 'bold', 'FontSize', 7);

yLabelStr = sprintf(['4-fold correlations\n' ...
    char(12296) char(968) '_4({\\itr}) '...
    char(968) '_4^*(0)' char(12297) ]);

ylabel(yLabelStr, ...
    'FontWeight', 'bold', 'FontSize', 5);
xlabel(sprintf('Distance\n(Mean cell lengths)'), ...
    'FontWeight', 'bold', 'FontSize', 7);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

% ratio = fig.Position(4) ./ fig.Position(3);
% fig.Position(3) = (4/3)*5;
% fig.Position(4) = 5; %ratio * fig.Position(3);

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = 4.5; %4;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

%% Generate 6-Fold Orientational Order Comparisons ========================
clear; close all; clc;

% Load Synthetic Results for Early Times ----------------------------------
tidx = 15;

% Radial bins should be the same for all time points
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'r');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'RMax');

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measL');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measW');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'cellSize');

measL_T15 = measL;
measW_T15 = measW;
cellSize_T15 = cellSize;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allPsi6Corr');

psi6Corr_T15 = mean(cell2mat(allPsi6Corr), 1);

clear allPsi6Corr measL measW cellSize

% Load Synthetic Results for Intermediate Times ---------------------------
tidx = 161;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measL');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measW');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'cellSize');

measL_T161 = measL;
measW_T161 = measW;
cellSize_T161 = cellSize;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allPsi6Corr');

psi6Corr_T161 = mean(cell2mat(allPsi6Corr), 1);

clear allPsi6Corr measL measW cellSize

% Load Synthetic Results for Late Times -----------------------------------
tidx = 234;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measL');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measW');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'cellSize');

measL_T234 = measL;
measW_T234 = measW;
cellSize_T234 = cellSize;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allPsi6Corr');

psi6Corr_T234 = mean(cell2mat(allPsi6Corr), 1);

clear allPsi6Corr measL measW cellSize

%--------------------------------------------------------------------------
% Generate Figures
%--------------------------------------------------------------------------

fig = figure('Color', [1 1 1], 'visible', 'on');

% Axis bounds
xLim = [0.5 RMax];
yLim = [1e-4, 0.5];

hold on

% Plot correlation functions
plot(r, psi6Corr_T15, 'LineWidth', 1.5);
plot(r, psi6Corr_T161, 'LineWidth', 1.5);
plot(r, psi6Corr_T234, 'LineWidth', 1.5);

% plot(r, abs(psi6Corr_T15), 'LineWidth', 1.5);
% plot(r, abs(psi6Corr_T161), 'LineWidth', 1.5);
% plot(r, abs(psi6Corr_T234), 'LineWidth', 1.5);

% Plot tissue widths
plot( ...
    mean([measW_T161/cellSize_T161, measW_T234/cellSize_T234]) .* [1 1], ...
    [yLim(1) 0.8*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );

% Plot algebraic decay guide line
plot(r, 0.1 .* r.^(-1/4), '--k', 'LineWidth', 1.5);

hold off

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

xlim(xLim);
ylim(yLim);

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

% ylabel('\langle\psi_6({\itr/L})\psi_6(0)\rangle', ...
%     'FontWeight', 'bold', 'FontSize', 7);
% xlabel('{\it r/L} (Mean cell lengths)', ...
%     'FontWeight', 'bold', 'FontSize', 7);

yLabelStr = sprintf(['6-fold correlations\n' ...
    char(12296) char(968) '_6({\\itr}) '...
    char(968) '_6^*(0)' char(12297) ]);

ylabel(yLabelStr, ...
    'FontWeight', 'bold', 'FontSize', 5);
xlabel(sprintf('Distance\n(Mean cell lengths)'), ...
    'FontWeight', 'bold', 'FontSize', 7);

xticks([1 10]);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

% ratio = fig.Position(4) ./ fig.Position(3);
% fig.Position(3) = (4/3)*5;
% fig.Position(4) = 5; %ratio * fig.Position(3);

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = 4.5; % 4;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

%% Generate Postional Correlation Function Comparisons ====================
clear; close all; clc;

% Load Synthetic Results for Early Times ----------------------------------
tidx = 15;

% Radial bins should be the same for all time points
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'r');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'RMax');

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measL');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measW');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'cellSize');

measL_T15 = measL;
measW_T15 = measW;
cellSize_T15 = cellSize;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allGRX');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allGRY');

GRX_T15 = mean(cell2mat(allGRX), 1);
GRY_T15 = mean(cell2mat(allGRY), 1);

clear allGRX allGRY measL measW cellSize

% Load Synthetic Results for Intermediate Times ---------------------------
tidx = 161;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measL');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measW');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'cellSize');

measL_T161 = measL;
measW_T161 = measW;
cellSize_T161 = cellSize;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allGRX');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allGRY');

GRX_T161 = mean(cell2mat(allGRX), 1);
GRY_T161 = mean(cell2mat(allGRY), 1);

clear allGRX allGRY measL measW cellSize

% Load Synthetic Results for Late Times -----------------------------------
tidx = 234;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measL');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'measW');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'cellSize');

measL_T234 = measL;
measW_T234 = measW;
cellSize_T234 = cellSize;

load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allGRX');
load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allGRY');

GRX_T234 = mean(cell2mat(allGRX), 1);
GRY_T234 = mean(cell2mat(allGRY), 1);

clear allGRX allGRY measL measW cellSize

%--------------------------------------------------------------------------
% Generate Figures
%--------------------------------------------------------------------------
plotDir = 'X';

fig = figure('Color', [1 1 1],  'visible', 'on');

hold on

% Plot averaged correlations
if strcmp(plotDir, 'X')

    plotGRX_T15 = GRX_T15-1; plotGRX_T15(plotGRX_T15 <= 0) = 1e-10;
    plotGRX_T161 = GRX_T161-1; plotGRX_T161(plotGRX_T161 <= 0) = 1e-10;
    plotGRX_T234 = GRX_T234-1; plotGRX_T234(plotGRX_T234 <= 0) = 1e-10;
    
    plot(r, plotGRX_T15, 'LineWidth', 1.5)
    plot(r, plotGRX_T161, 'LineWidth', 1.5)
    plot(r, plotGRX_T234, 'LineWidth', 1.5)
    
    % Y-axis label
    % yStr = ['g({\itr/L}, ' char(952) '_{DV})'];
    % yStr = sprintf('Radial distribution\nalong DV-axis');
    yStr = sprintf('DV pair correlations\n g_{DV}(r)-1');
    
    % Axis bounds
    xLim = [0.6 RMax];
    yLim = [5e-3 1];
    
    % Coefficient of algebraic decay guideline
    algCoeff = 0.21;
    
elseif strcmp(plotDir, 'Y')

    plotGRY_T15 = GRY_T15-1; plotGRY_T15(plotGRY_T15 <= 0) = 1e-10;
    plotGRY_T161 = GRY_T161-1; plotGRY_T161(plotGRY_T161 <= 0) = 1e-10;
    plotGRY_T234 = GRY_T234-1; plotGRY_T234(plotGRY_T234 <= 0) = 1e-10;
    
    plot(r, plotGRY_T15, 'LineWidth', 1.5)
    plot(r, plotGRY_T161, 'LineWidth', 1.5)
    plot(r, plotGRY_T234, 'LineWidth', 1.5)
    
    % Y-axis label
    % yStr = ['g({\itr/L}, ' char(952) '_{AP})'];
    % yStr = sprintf('Radial distribution\nalong AP-axis');
    yStr = sprintf('AP pair correlations\n g_{AP}(r)-1');
    
    % Axis bounds
    xLim = [0.6 RMax];
    yLim = [5e-3 1];
    
    algCoeff = 0.2;
    
else
    
    error('Invalid plot direction');
    
end

% Plot tissue lengths
plot( ...
    mean([measL_T161/cellSize_T161, measL_T234/cellSize_T234]) .* [1 1], ...
    [yLim(1) 0.8*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );

% Plot tissue widths
plot( ...
    mean([measW_T161/cellSize_T161, measW_T234/cellSize_T234]) .* [1 1], ...
    [yLim(1) 0.8*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );

% Plot algebraic decay guide line
% plot(r, algCoeff .* r.^(-1/3), '--k', 'LineWidth', 1.5);
% plot(r, algCoeff .* r.^(-1/5), '--k', 'LineWidth', 1.5);

hold off

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

xlim(xLim);
ylim(yLim);

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

xticks([1 10]);

ylabel(yStr, 'FontWeight', 'bold', 'FontSize', 5);
% xlabel('{\it r/L} (Mean cell lengths)', ...
%     'FontWeight', 'bold', 'FontSize', 7);
xlabel(sprintf('Distance\n(Mean cell lengths)'), ...
    'FontWeight', 'bold', 'FontSize', 7);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

% ratio = fig.Position(4) ./ fig.Position(3);
% fig.Position(3) = (4/3)*5;
% fig.Position(4) = 5; %ratio * fig.Position(3);

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = 4.5;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

%% Generate Order Parameter Magnitdue Comparisons =========================
clear; close all; clc;

tidx = 15;
orderType = 4;
displayType = 'magnitude';

if strcmpi(displayType, 'magnitude')
    
    savename = sprintf(['ValidationFigures/Hist_T%d_Abs_Psi%d-' ...
        date '.png'], tidx, orderType);
    
elseif strcmpi(displayType, 'angle')
    
    savename = sprintf(['ValidationFigures/Hist_T%d_Ang_Psi%d-' ...
        date '.png'], tidx, orderType);
    
else
    
    error('Invalid display type');
    
end

if (orderType == 4)
    
    % Load synthetic results ----------------------------------------------
    load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allIncNodeIDx');
    load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allPsi4');
    
    synthPsi = [];
    for i = 1:length(allIncNodeIDx)
        synthPsi = [ synthPsi; allPsi4{i}(allIncNodeIDx{i}) ];
    end
    
    clear allIncNodeIDx allPsi4 i
    
    % Load measured results -----------------------------------------------
    load('Correlation_Function_Validation_20210303.mat', 'allIncNodeIDx');
    load('Correlation_Function_Validation_20210303.mat', 'smoothPsi4');
    
    switch tidx
        
        case 15
            
            measPsi = smoothPsi4(vertcat(allIncNodeIDx{10:20}));
            
        case 161
            
            measPsi = smoothPsi4(vertcat(allIncNodeIDx{159:163}));
            
        case 234
            
            measPsi = smoothPsi4(vertcat(allIncNodeIDx{231:234}));
            
        otherwise
            
            error('Invalid time point');
    end
    
    clear allIncNodeIDx smoothPsi4
    
elseif (orderType == 6)
    
    % Load synthetic results ----------------------------------------------
    load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allIncNodeIDx');
    load(sprintf('Random_Correlations_T%03d.mat', tidx), 'allPsi6');
    
    synthPsi = [];
    for i = 1:length(allIncNodeIDx)
        synthPsi = [ synthPsi; allPsi6{i}(allIncNodeIDx{i}) ];
    end
    
    clear allIncNodeIDx allPsi6 i
    
    % Load measured results -----------------------------------------------
    load('Correlation_Function_Validation_20210303.mat', 'allIncNodeIDx');
    load('Correlation_Function_Validation_20210303.mat', 'smoothPsi6');
    
    switch tidx
        
        case 15
            
            measPsi = smoothPsi6(vertcat(allIncNodeIDx{10:20}));
            
        case 161
            
            measPsi = smoothPsi6(vertcat(allIncNodeIDx{159:163}));
            
        case 234
            
            measPsi = smoothPsi6(vertcat(allIncNodeIDx{231:234}));
            
        otherwise
            
            error('Invalid time point');
    end
    
    clear allIncNodeIDx smoothPsi6
    
else
    
    error('Invalid order type');
    
end


%--------------------------------------------------------------------------
% Generate Figures
%--------------------------------------------------------------------------

fig = figure('Color', [1 1 1]);

hold on

if strcmpi(displayType, 'magnitude')
    
    binEdges = linspace(0, 1, 25);
    histogram( abs(synthPsi), binEdges, 'Normalization', 'probability' );
    histogram( abs(measPsi), binEdges, 'Normalization', 'probability' );
    
    xlabel(sprintf('\\mid \\psi_%d \\mid', orderType));
    ylabel('Probability');
    
    xTicks = [0 0.5 1];
    xTickStr = {'0', '0.5', '1'};
    % yTick = [0 0.5 1];
    % yTickStr = {'0', '0.5', '1'};
    
    
elseif strcmpi(displayType, 'angle')
    
    binEdges = linspace(-pi, pi, 25);
    histogram( angle(synthPsi), binEdges, 'Normalization', 'probability' );
    histogram( angle(measPsi), binEdges, 'Normalization', 'probability' );
    
    xlabel(sprintf('arg[\\psi_%d]', orderType));
    ylabel('Probability');
    
    xTicks = [-pi 0 pi];
    xTickStr = {'-\pi', '0', '\pi'};
    % yTicks = [0 0.5 1];
    % yTickStr = {'0', '0.5', '1'};
    
else
    
    error('Invalid display type');
    
end

hold off

histAx = gca;

xticks(xTicks);
xticklabels(xTickStr);
% yticks(yTicks);
% yticklabels(yTickStr);

set(histAx, 'FontSize', 6);
set(histAx, 'FontWeight', 'bold');
set(histAx, 'LineWidth', 1);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = 4;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

% Save figure
% print( savename );


