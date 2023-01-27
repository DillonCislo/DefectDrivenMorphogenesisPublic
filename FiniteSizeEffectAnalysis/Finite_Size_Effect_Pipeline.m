%% Finite Size Effect Analysis ==========================================
% This is a script to analyze finite size effects on our definition of
% the orientational order parameter, and the orientational/translational
% correlation functions. It can be used to generate the analysis/figure
% panels found in Fig. S8
%
% by Dillon Cislo 2022/06/08
% =======================================================================

% Add relevant file structure to MATLAB path
[DDMDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
DDMDir = [DDMDir '/..'];
addpath(genpath(DDMDir));
clear DDMDir

%% Prepare Simulations ==================================================
clear; close all; clc;

% Generate Output File Structure ----------------------------------------

% The directory where project files will be generated and saved
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

cd(projectDir);

saveDir = 'Finite_Size_Effect_Analysis_08-Jun-2022'; % Use with data from Dryad
% saveDir = ['Finite_Size_Effect_Analysis_' date];
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

RFileBase = fullfile(saveDir, 'R_SF%d_SNR%d.csv');

XYFileBase = fullfile(saveDir, 'XY_SF%d_SNR%d_N%d.csv');
IncNodeFileBase = fullfile(saveDir, 'IncNodeIDx_SF%d_SNR%d_N%d.csv');

Psi4FileBase = fullfile(saveDir, 'Psi4_SF%d_SNR%d_N%d.csv');
Psi6FileBase = fullfile(saveDir, 'Psi6_SF%d_SNR%d_N%d.csv');

Psi4CorrFileBase = fullfile(saveDir, 'Psi4Corr_SF%d_SNR%d_N%d.csv');
Psi6CorrFileBase = fullfile(saveDir, 'Psi6Corr_SF%d_SNR%d_N%d.csv');
NumRCellsFileBase = fullfile(saveDir, 'NumRCells_SF%d_SNR%d_N%d.csv');

GRXFileBase = fullfile(saveDir, 'GRX_SF%d_SNR%d_N%d.csv');
GRYFileBase = fullfile(saveDir, 'GRY_SF%d_SNR%d_N%d.csv');
GRRFileBase = fullfile(saveDir, 'GRR_SF%d_SNR%d_N%d.csv');

figDir = fullfile(projectDir, 'Finite_Size_Analysis_Figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end


% Set Simulation Metadata -----------------------------------------------

cellSize = 1; % Average length scale of cells

% Distances for correlation functions will be normalized by cell sizes 
normDists = true;

% The dimensions of the basic domain used to analyze the finite size
% effects
basicIncBoxX = cellSize * 5 * [-1 1];
basicIncBoxY = cellSize * 5 * [-1 1];

% The scale factors used to generate the different domains
scaleFactors = [1 2 4 8];

% The number of samples to generate per category
numSamples = 25;

% The SNR ratios used to generate the noisy configurations
allSNR = [100 25 15 0];

%% Run Simulations ======================================================

% For reproducible random numbers
rng(1, 'twister');

for simID = 1:(numel(scaleFactors) * numel(allSNR) * numSamples)

    [i, j, k] = ...
        ind2sub([numel(scaleFactors), numel(allSNR), numSamples], simID);

    SF = scaleFactors(i); SNR = allSNR(j);

    fprintf('Running Simulation: Scale = %d, SNR = %d, N = %d\n', ...
        SF, SNR, k);

    %--------------------------------------------------------------------
    % Handle Output File Structure
    %--------------------------------------------------------------------

    RFile = sprintf(RFileBase, SF, SNR);

    XYFile = sprintf(XYFileBase, SF, SNR, k);
    IncNodeFile = sprintf(IncNodeFileBase, SF, SNR, k);

    Psi4File = sprintf(Psi4FileBase, SF, SNR, k);
    Psi6File = sprintf(Psi6FileBase, SF, SNR, k);

    Psi4CorrFile = sprintf(Psi4CorrFileBase, SF, SNR, k);
    Psi6CorrFile = sprintf(Psi6CorrFileBase, SF, SNR, k);
    NumRCellsFile = sprintf(NumRCellsFileBase, SF, SNR, k);

    GRXFile = sprintf(GRXFileBase, SF, SNR, k);
    GRYFile = sprintf(GRYFileBase, SF, SNR, k);
    GRRFile = sprintf(GRRFileBase, SF, SNR, k);

    allFilesExist = exist(RFile, 'file');

    allFilesExist = allFilesExist && exist(XYFile, 'file');
    allFilesExist = allFilesExist && exist(IncNodeFile, 'file');

    allFilesExist = allFilesExist && exist(Psi4File, 'file');
    allFilesExist = allFilesExist && exist(Psi6File, 'file');
    allFilesExist = allFilesExist && exist(NumRCellsFile, 'file');

    allFilesExist = allFilesExist && exist(GRXFile, 'file');
    allFilesExist = allFilesExist && exist(GRYFile, 'file');
    allFilesExist = allFilesExist && exist(GRRFile, 'file');

    if allFilesExist, continue; end

    %--------------------------------------------------------------------
    % Generate Point Set
    %--------------------------------------------------------------------

    % Points are generated inside a large box. Only points within a smaller
    % box, itself contained within the large box, will be included in
    % calculations of the correlation functions
    incBoxX = SF * basicIncBoxX;
    incBoxY = SF * basicIncBoxY;
    bigBoxX = incBoxX + cellSize * 5 * [-1 1];
    bigBoxY = incBoxY + cellSize * 5 * [-1 1];

    % Generate distance bin structure for cell domain
    RMax = ceil((sqrt(2) * diff(bigBoxX)) / cellSize);
    r = 0:0.25:RMax;
    nBins = numel(r);

    % if (k == 1), allR{simID} = r; end
    if (k == 1), writematrix(r, RFile); end

    latticeOptions = struct();
    latticeOptions.xLim = bigBoxX;
    latticeOptions.yLim = bigBoxY;

    if (SNR == 0)

        latticeOptions.pointSetType = 'randomMinDist';
        latticeOptions.minDist = 0.8 * cellSize;
        latticeOptions.maxRandSamples = 30;
        latticeOptions.numPoints = Inf;

    else

        latticeOptions.pointSetType = 'noisySquare';
        latticeOptions.SNR = SNR;
        latticeOptions.sideLength = cellSize;

    end

    XY = CELTIGS.generate_points( latticeOptions );

    % Determine which points are included in the small box
    incNodeIDx = (incBoxX(1) <= XY(:,1)) & (XY(:,1) <= incBoxX(2)) ...
        & (incBoxY(1) <= XY(:,2)) & (XY(:,2) <= incBoxY(2));
    incNodeIDx = find(incNodeIDx);

    fprintf('N = %d cells generated\n', numel(incNodeIDx));

    writematrix(XY, XYFile);
    writematrix(incNodeIDx, IncNodeFile);

    %--------------------------------------------------------------------
    % Calculate Bond Orientational Order/Correlations
    %--------------------------------------------------------------------

    % NOTE: I AM CHOOSING TO IMPOSE SOME SPATIAL SMOOTHING ON THE
    % ORDER PARAMETERS!!!!
    psi4 = calculateBondOrientationalOrder(XY, 4, 2, true, false);
    psi6 = calculateBondOrientationalOrder(XY, 6, 2, true, false);

    [psi4Corr, ~, numCellsR, ~] = ...
        calculateBondOrientationalCorrelations( XY, psi4, ...
        'CellSize', cellSize, 'MaxR', RMax, 'NumRBins', nBins, ...
        'IncludedNodes', incNodeIDx, 'OrderType', 4, ...
        'NormalizeDistances', normDists, 'PlotResults', false );

    [psi6Corr, ~, ~, ~] = ...
        calculateBondOrientationalCorrelations( XY, psi6, ...
        'CellSize', cellSize, 'MaxR', RMax, 'NumRBins', nBins, ...
        'IncludedNodes', incNodeIDx, 'OrderType', 6, ...
        'NormalizeDistances', normDists, 'PlotResults', false );

    writematrix(psi4, Psi4File);
    writematrix(psi6, Psi6File);
    writematrix(psi4Corr, Psi4CorrFile);
    writematrix(psi6Corr, Psi6CorrFile);
    writematrix(numCellsR, NumRCellsFile);

    %--------------------------------------------------------------------
    % Calculate Pair Correlation Functions
    %--------------------------------------------------------------------

    useEffectiveVolume = true;
    effectiveVolumeMethod = 'polygon';

    [~, GR_X, GR_Y, GR_R, ~] = calculatePositionalCorrelations( ...
        XY, 'CellSize', cellSize, 'MaxR', RMax, 'NumRBins', nBins, ...
        'IncludedNodes', incNodeIDx, 'NormalizeDistances', normDists, ...
        'EffectiveVolume', useEffectiveVolume, ...
        'EffectiveVolumeMethod', effectiveVolumeMethod );

    writematrix(GR_X, GRXFile);
    writematrix(GR_Y, GRYFile);
    writematrix(GR_R, GRRFile);

end

clear simID


%% **********************************************************************
% ***********************************************************************
%                       GENERATE FIGURES
% ***********************************************************************
% ***********************************************************************

% Running the above pipeline with the default supplied setting will
% generate the figures from the paper (modulo some very small differences
% perhaps since we use the MATLAB built-in method for calculating polygon
% intersection areas, rather than the Boost method), but will take a long
% time.
%
% If you want to generate the figures directly, download the folder
% 'Finite_Size_Effect_Analysis_08-Jun-2022' from the Dryad host repository
% located at: 
%
% In the first cell block make sure that  you set
%
%   saveDir = 'Finite_Size_Effect_Analysis_08-Jun-2022';
%
% and run accordingly.

%% Generate Mean Orientational Order Figure =============================
close all; clc;

dispSNR = 15;
SNRID = find(allSNR == dispSNR);
orderType = 4;

% Load Average Order Parameters -----------------------------------------

allMeanPsi = nan(numel(scaleFactors), numSamples);
for i = 1:numel(scaleFactors)
    for k = 1:numSamples

        if orderType == 4
            psiFile = sprintf(Psi4FileBase, scaleFactors(i), dispSNR, k);
        else
            psiFile = sprintf(Psi6FileBase, scaleFactors(i), dispSNR, k);
        end

        incNodeFile = ...
            sprintf(IncNodeFileBase, scaleFactors(i), dispSNR, k);

        if exist(psiFile, 'file')
            psi = readmatrix(psiFile);
            incNodeIDx = readmatrix(incNodeFile);
            allMeanPsi(i,k) = mean(psi(incNodeIDx));
        end

    end
end

clear i k psiFile incNodeFile psi incNodeIDx

meanPsi = mean(allMeanPsi, 2, 'omitnan');
stdPsi = std(allMeanPsi, 0, 2, 'omitnan');

% Generate Figure -------------------------------------------------------
close all; clc;
scaleColors = brewermap(4, 'Dark2');

sysSize = 10 * scaleFactors;

fig = figure('Color', [1 1 1], 'visible', 'on');

violinplot(abs(allMeanPsi).', sysSize, 'ViolinColor', scaleColors, ...
    'ShowMean', true, 'ShowData', false);

% hold on
% for i = 1:numel(scaleFactors)
% 
%     scatter( repmat(sysSize(i), sum(~isnan(allMeanPsi(i,:))), 1), ...
%         abs(allMeanPsi(i, ~isnan(allMeanPsi(i,:)))), ...
%         'x', 'Color', scaleColors(i,:));
% 
%     errorbar( sysSize(i), abs(meanPsi(i)), abs(stdPsi(i)), 'o', ...
%         'Color', 'k', 'MarkerFaceColor', 'k');
% 
% 
% end
% hold off

%xlim([0 100]);

if dispSNR == 15
    yticks(0.3:0.1:0.6);
    ylim([0.3 0.6]);
elseif dispSNR == 0
    yticks(0:0.1:0.2);
    ylim([0 0.2]);
end

% ylabel(sprintf('\\mid \\langle \\psi_%d \\rangle \\mid', orderType));
ylabel(['\mid \langle' char(968) '_' num2str(orderType) '\rangle \mid']);
xlabel('System size (Mean cell lengths)')

% Re-size Figure for Paper ----------------------------------------------

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = 5; 
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

% Save Results ----------------------------------------------------------

saveFileName = fullfile(figDir, ['Psi%d_SNR%d_' date '.png']);
saveFileName = sprintf(saveFileName, orderType, dispSNR);
print(saveFileName);

%% Generate Pair Correlation Function Comparison Figure =================
close all; clc;

dispSNR = 0;
% SNRID = find(allSNR == dispSNR);
orderType = 'X';

switch orderType
    case 'R'
        
        switch dispSNR
            case 0, yLim = [0 1];
            case 15, yLim = [0 1];
            case 25, yLim = [0 2];
            case 100, yLim = [0 2];
        end

        yLabelStr = sprintf('Pair correlation function\n g(r)-1');

    case 'X'

        switch dispSNR
            case 0, yLim = [0 1];
            case 15, yLim = [0 2];
            case 25, yLim = [0 5];
            case 100, yLim = [0 5];
        end

        yLabelStr = sprintf('Pair correlation function\n g(r,x)-1');

    case 'Y'

        switch dispSNR
            case 0, yLim = [0 2];
            case 15, yLim = [0 2];
            case 25, yLim = [0 5];
            case 100, yLim = [0 5];
        end

        yLabelStr = sprintf('Pair correlation function\n g(r,y)-1');

    otherwise
        error('Invalid order type')
end

% Load Correlation Functions --------------------------------------------

% Load distance bin positions
allR = cell(numel(scaleFactors), 1);
for i = 1:numel(scaleFactors)
    RFile = sprintf(RFileBase, scaleFactors(i), dispSNR);
    allR{i} = readmatrix(RFile);
end

allG = cell(numel(scaleFactors), 1);
for i = 1:numel(scaleFactors)

    count = 0;
    curG = zeros(1, numel(allR{i}));
    for k = 1:numSamples

        switch orderType
            case 'R'
                GFile = sprintf(GRRFileBase, scaleFactors(i), dispSNR, k);
            case 'X'
                GFile = sprintf(GRXFileBase, scaleFactors(i), dispSNR, k);
            case 'Y'
                GFile = sprintf(GRYFileBase, scaleFactors(i), dispSNR, k);
            otherwise
                error('Invalid order type');
        end

        if exist(GFile, 'file')
            curG = curG + readmatrix(GFile);
            count = count+1;
        end

    end

    curG = curG ./ count;
    allG{i} = curG;

end

clear i k RFile count curG 

% Generate Figure -------------------------------------------------------
close all; clc;
scaleColors = flipud(brewermap(4, 'Dark2'));
sysSize = 10 * scaleFactors;

corrLineWidth = 1;
sysLineWidth = 2;

generateLogPlot = true;
if generateLogPlot, yLim(1) = 1e-2; end

fig = figure('Color', [1 1 1], 'visible', 'on');

hold on

for i = numel(scaleFactors):-1:1

    if generateLogPlot
        plotG = allG{i}-1;
        plotG(plotG <= 0) = 1e-10;
    else
        plotG = allG{i};
    end

    plot(allR{i}, plotG, ...
        'LineWidth', corrLineWidth, 'Color', scaleColors(i, :));

    plot(sysSize(i) * [1 1], yLim, ':', ...
        'LineWidth', sysLineWidth, 'Color', scaleColors(i, :));

end

hold off

if generateLogPlot
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xticks([1 10 100]);
end

ylim(yLim)
xlim([0.5 100]);

ylabel(yLabelStr);
xlabel('Distance (Mean cell lengths)');

% Re-size Figure for Paper ----------------------------------------------

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = 5;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

% Save Results ----------------------------------------------------------

saveFileName = fullfile(figDir, ['G' orderType '_SNR%d_' date '.eps']);
saveFileName = sprintf(saveFileName, dispSNR);
print(saveFileName);

%% Generate Orientational Order Correlation Figures =====================
close all; clc;

dispSNR = 0;
% SNRID = find(allSNR == dispSNR);

orderType = 4;
yLabelStr = ['Correlation function\n' char(12296) char(968) '_%d({\\itr}) '...
    char(968) '_%d^*(0)' char(12297) ];
yLabelStr = sprintf(yLabelStr, orderType, orderType);

switch orderType

    case 4
        
        switch dispSNR
            case 0, yLim = [0 0.1];
            case 15, yLim = [0 0.3];
            case 25, yLim = [0 1];
            case 100, yLim = [0 1];
        end

    case 6

        switch dispSNR
            case 0, yLim = [0 0.1];
            case 15, yLim = [0 0.1];
            case 25, yLim = [0 0.1];
            case 100, yLim = [0 0.1];
        end

    otherwise

        error('Invalid order type');

end

% Load Correlation Functions --------------------------------------------

% Load distance bin positions
allR = cell(numel(scaleFactors), 1);
for i = 1:numel(scaleFactors)
    RFile = sprintf(RFileBase, scaleFactors(i), dispSNR);
    allR{i} = readmatrix(RFile);
end

allPsiCorr = cell(numel(scaleFactors), 1);
for i = 1:numel(scaleFactors)

    count = 0;
    curPsiCorr = zeros(1, numel(allR{i}));
    for k = 1:numSamples

        switch orderType
            case 4
                psiCorrFile = sprintf(Psi4CorrFileBase, scaleFactors(i), dispSNR, k);
            case 6
                psiCorrFile = sprintf(Psi6CorrFileBase, scaleFactors(i), dispSNR, k);
            otherwise
                error('Invalid order type');
        end

        if exist(psiCorrFile, 'file')
            curPsiCorr = curPsiCorr + readmatrix(psiCorrFile);
            count = count+1;
        end

    end

    curPsiCorr = curPsiCorr ./ count;
    allPsiCorr{i} = curPsiCorr;

end

clear i k RFile count curPsiCorr  psiCorrFile 

% Generate Figure -------------------------------------------------------
close all; clc;
scaleColors = flipud(brewermap(4, 'Dark2'));
sysSize = 10 * scaleFactors;

corrLineWidth = 1;
sysLineWidth = 2;

generateLogPlot = true;
if generateLogPlot, yLim(1) = 1e-2; end

fig = figure('Color', [1 1 1], 'visible', 'on');

hold on

for i = numel(scaleFactors):-1:1

    plotPsiCorr = allPsiCorr{i};
    if generateLogPlot, plotPsiCorr(plotPsiCorr <= 0) = 1e-10; end

    plot(allR{i}, plotPsiCorr, ...
        'LineWidth', corrLineWidth, 'Color', scaleColors(i, :));

    plot(sysSize(i) * [1 1], yLim, ':', ...
        'LineWidth', sysLineWidth, 'Color', scaleColors(i, :));

end

hold off

if generateLogPlot
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xticks([1 10 100]);
end

ylim(yLim)
xlim([0.5 120]);

ylabel(yLabelStr);
xlabel('Distance (Mean cell lengths)');

% Re-size Figure for Paper ----------------------------------------------

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = 5;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

% Save Results ----------------------------------------------------------

saveFileName = fullfile(figDir, ['Psi%dCorr_SNR%d_' date '.eps']);
saveFileName = sprintf(saveFileName, orderType, dispSNR);
print(saveFileName);