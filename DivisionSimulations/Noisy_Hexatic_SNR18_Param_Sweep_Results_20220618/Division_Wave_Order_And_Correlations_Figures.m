%% Division Wave Order and Correlations Figures Pipeline ==================
%
%   This is a pipeline to generate a set of figures showing the
%   orientational order, orientational order correlations, and positional
%   order functions of the final configuration of a set of simulations
%   after a single division wave
%
%   In particular, this script will generate Fig. 5D, Fig. 5E, Fig. S19B
%   and Fig. S19C
%
%   by Dillon Cislo 05/17/2022
%==========================================================================

% Add relevant file structure to MATLAB path
[DDMDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
DDMDir = [DDMDir '/..'];
addpath(genpath(DDMDir));
clear DDMDir

%% Set Up Pipeline ========================================================
clear; close all; clc;

% Set Division Distribution Paramters -------------------------------------
numSimK = 5;
allK = [0 2 6 20 Inf].';
cellSize = 1;

% Set Directory Data ------------------------------------------------------

% The working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% The directory holding the division wave simulation results
simResultDir = projectDir;

% The simulation file name template
simResultFileBase = 'Simulation_Result_K%d_N%d.mat';
simResultFileBase = fullfile(simResultDir, simResultFileBase);

% The output directory
saveDir = projectDir;

%% Calculate Parameters for All Simulation Runs ===========================
close all; clc;

% The total number of simulations to run
numSim = numSimK * numel(allK);

allPsi4 = cell(numSimK, numel(allK));
allPsi6 = cell(numSimK, numel(allK));

allPsi4Corr = cell(numSimK, numel(allK));
allPsi6Corr = cell(numSimK, numel(allK));
allNumCellsR = cell(numSimK, numel(allK));

allGRMax = cell(numSimK, numel(allK));
allGRX = cell(numSimK, numel(allK));
allGRY = cell(numSimK, numel(allK));
allGRR = cell(numSimK, numel(allK));

allCellCentroids = cell(numSimK, numel(allK));

RMax = 50;
numRBins = 100;
normDists = true;
useEffectiveVolume = true;
r = linspace(0, RMax, numRBins);

parfor simID = 1:numSim
    
    % Convert linear simulation index to (simK, K) format
    [kSimID, kIDx] = ind2sub([numSimK, numel(allK)], simID);
    fprintf('Running simulation %d for k = %d\n', kSimID, allK(kIDx));
    
    %======================================================================
    % Load Data For All Time Points of Current Simulation Run
    %======================================================================
    
    k = allK(kIDx);
    
    simResultFile = sprintf(simResultFileBase, k, kSimID);
    simResult = load(simResultFile);
    allGLattices = simResult.allGLattices;
    divOrder = simResult.divOrder;
    g0 = simResult.g0;
    latticeOptions = simResult.latticeOptions;
    
    %======================================================================
    % Extract Fields for Each Time Point
    %======================================================================
    
    simPsi4 = zeros(numel(allGLattices), 2);
    simPsi6 = zeros(numel(allGLattices), 2);
    
    simPsi4Corr = zeros(numel(allGLattices), numRBins);
    simPsi6Corr = zeros(numel(allGLattices), numRBins);
    simNumCellsR = zeros(numel(allGLattices), numRBins);
    
    simGRMax = zeros(numel(allGLattices), numRBins);
    simGRX = zeros(numel(allGLattices), numRBins);
    simGRY = zeros(numel(allGLattices), numRBins);
    simGRR = zeros(numel(allGLattices), numRBins);
    
    simCellCentroids = cell(numel(allGLattices), 1);
    
    for tidx = 1:numel(allGLattices)
        
        % The current lattice
        g = allGLattices{tidx};
        
        warning('off','all')
        cellCentroids = zeros(numel(g.cells), 2);
        for i = 1:numel(g.cells)
            
            X = g.verts(g.bonds(g.cells{i},1),1);
            Y = g.verts(g.bonds(g.cells{i},1),2);
            [cellCentroids(i,1), cellCentroids(i,2)] = ...
                centroid(polyshape(X,Y));
            
        end
        warning('on','all')
        
        simCellCentroids{tidx} = cellCentroids;
        
        % clear X Y simResultFile
        
        %------------------------------------------------------------------
        % Generate a set of boundary points around the current
        % configuration
        %------------------------------------------------------------------
        
        % METHOD 1: Simple Bounding Box -----------------------------------
        
        W = range(g.verts(:,1));
        H = range(g.verts(:,2));
        % boxScale = 1.1;
        % WW = boxScale * W; HH = boxScale * H;
        WW = W + cellSize; HH = H + 0.35 * cellSize;
        P = 2 * (WW+HH);
        
        numBdyPoints = ceil(P / cellSize);
        
        linS = linspace(0, P, numBdyPoints+1).';
        linS(end) = [];
        
        topIDx = linS < WW;
        rightIDx = (WW <= linS) & (linS < (WW+HH));
        bottomIDx = ((WW+HH) <= linS) & (linS < (2*WW+HH));
        leftIDx = (2*WW+HH) <= linS;
        
        linS(rightIDx) = linS(rightIDx) - WW;
        linS(bottomIDx) = linS(bottomIDx) - (WW+HH);
        linS(leftIDx) = linS(leftIDx) - (2*WW+HH);
        
        bdyXY = zeros(numBdyPoints, 2);
        bdyXY(topIDx, :) = [(-WW/2) + linS(topIDx), ...
            (HH/2) * ones(sum(topIDx), 1)];
        bdyXY(rightIDx, :) = [(WW/2) * ones(sum(rightIDx), 1), ...
            (HH/2) - linS(rightIDx)];
        bdyXY(bottomIDx, :) = [(WW/2) - linS(bottomIDx), ...
            (-HH/2) * ones(sum(bottomIDx), 1)];
        bdyXY(leftIDx, :) = [(-WW/2) * ones(sum(leftIDx), 1), ...
            (-HH/2) + linS(leftIDx)];
        
        bdyXY(:,1) = bdyXY(:,1) + (min(g.verts(:,1)) + W/2);
        bdyXY(:,2) = bdyXY(:,2) + (min(g.verts(:,2)) + H/2);
        
        % clear W H boxScale WW HH P numBdyPoints linS
        % clear topIDx rightIDx bottomIDx leftIDx
        
        % METHOD 2: Convex Hull -------------------------------------------
        
        % hullScale = 1.00;
        % hullIDx = convhull(g.verts(:,1), g.verts(:,2));
        % hullIDx(end) = [];
        % hullIDx = [hullIDx, circshift(hullIDx, [-1 0])];
        % hullXY = hullScale * g.verts(hullIDx(:,1), 1:2);
        %
        % % Calculate polygon edge lengths
        % hullL = g.verts(hullIDx(:,2), :) - g.verts(hullIDx(:,1), :);
        % hullL = sqrt(sum(hullL.^2, 2));
        %
        % % Interpolate to find the new points
        % bdyX = interp1([0; cumsum(hullL)], [hullXY(:,1); hullXY(1,1)], ...
        %     (0:cellSize:sum(hullL)).');
        % bdyY = interp1([0; cumsum(hullL)], [hullXY(:,2); hullXY(1,2)], ...
        %     (0:cellSize:sum(hullL)).');
        % bdyX(1) = []; bdyY(1) = [];
        % bdyXY = [bdyX, bdyY];
        %
        % % Calculate outward pointing normal of each new boundary point
        % eN = [circshift(bdyXY, [-1 0]) - bdyXY, zeros(size(bdyXY,1),1)];
        % eN = cross(eN, repmat([0 0 1], size(eN,1), 1), 2);
        % eN = eN ./ sqrt(sum(eN.^2, 2));
        % eN = eN(:, 1:2);
        % eN2 = circshift(eN, [-1 0]);
        % eNAngles = atan2(eN(:,2), eN(:,1));
        % bendingAngles = 2 * atan2( eN(:,1) .* eN2(:,2) - eN(:,2) .* eN2(:,1), ...
        %     1 + dot(eN, eN2, 2) );
        % vNAngles = eNAngles + bendingAngles/2;
        % vN = [cos(vNAngles), sin(vNAngles)];
        %
        % % Displace points outward
        % bdyXY = bdyXY + 0.5 .* cellSize .* vN;
        %
        % % clear hullScale hullIDx hullXY hullL bdyX bdyY eN eN2
        % % clear eNAngles bendingAngles vNAngles vN
        
        %------------------------------------------------------------------
        % Calculate Order Parameter Correlations
        %------------------------------------------------------------------
        allXY = [cellCentroids; bdyXY];
        
        % IDs of the virtual boundary cells added to the tissue
        virtualBdyIDx = (1:size(bdyXY,1)).' + size(cellCentroids,1);
        
        % IDs of the original boundary cells of the tissue
        % bdyBondIDx = find((g.bonds(:,3) == 0) | (g.bonds(:,4) == 0));
        % originalBdyIDx = find(cellfun(@(x) any(ismember(x, bdyBondIDx)), g.cells));
        
        divOrderIDx = [divOrder; ...
            numel(allGLattices{1}.cells)+ ...
            (1:(numel(g.cells)-numel(allGLattices{1}.cells))).'];
        
        psi4 = calculateBondOrientationalOrder(allXY, 4, 2, false, false);
        psi6 = calculateBondOrientationalOrder(allXY, 6, 2, false, false);
        
        dispCells = ~isnan(psi4) & ~isnan(psi6);
        dispCells = dispCells & ismember((1:numel(dispCells)).', divOrderIDx);
        % dispCells = dispCells & ~ismember((1:numel(dispCells)).', virtualBdyIDx);
        % dispCells = dispCells & ~ismember((1:numel(dispCells)).', originalBdyIDx);

        incNodeIDx = find(dispCells);
        
        % Calculate the bond orientational correlation functions
        [psi4Corr, ~, numCellsR, ~] = ...
            calculateBondOrientationalCorrelations( allXY, psi4, ...
            'MaxR', RMax, 'NumRBins', numRBins, 'OrderType', 4, ...
            'IncludedNodes', incNodeIDx, 'CellSize', cellSize, ...
            'NormalizeDistances', normDists, 'PlotResults', false );
        
        [psi6Corr, ~, ~, ~] = ...
            calculateBondOrientationalCorrelations( allXY, psi6, ...
            'MaxR', RMax, 'NumRBins', numRBins, 'OrderType', 6, ...
            'IncludedNodes', incNodeIDx, 'CellSize', cellSize, ...
            'NormalizeDistances', normDists, 'PlotResults', false );
        
        % Calculate positional correlation functions
        [GRMax, GRX, GRY, GRR, ~] = calculatePositionalCorrelations( ...
            allXY, 'MaxR', RMax, 'NumRBins', numRBins, ...
            'IncludedNodes', incNodeIDx, ...
            'EffectiveVolume', useEffectiveVolume, ...
            'CellSize', cellSize, 'ProgressBar', false, ...
            'NormalizeDistances', normDists, 'PlotResults', false );
        
        % Update ensemble arrays
        simPsi4(tidx,:) = [mean(psi4(dispCells)), ...
            std(psi4(dispCells)) ./ sqrt(sum(dispCells))];
        simPsi6(tidx,:) = [mean(psi6(dispCells)), ...
            std(psi6(dispCells)) ./ sqrt(sum(dispCells))];
        simPsi4Corr(tidx,:) = psi4Corr;
        simPsi6Corr(tidx,:) = psi6Corr;
        simNumCellsR(tidx,:) = numCellsR;
        simGRMax(tidx,:) = GRMax;
        simGRX(tidx,:) = GRX;
        simGRY(tidx,:) = GRY;
        simGRR(tidx,:) = GRR;
        simCellCentroids{simID} = cellCentroids;
        
    end
    
    % Update ensemble arrays
    allPsi4{simID} = simPsi4;
    allPsi6{simID} = simPsi6;
    allPsi4Corr{simID} = simPsi4Corr;
    allPsi6Corr{simID} = simPsi6Corr;
    allNumCellsR{simID} = simNumCellsR;
    allGRMax{simID} = simGRMax;
    allGRX{simID} = simGRX;
    allGRY{simID} = simGRY;
    allGRR{simID} = simGRR;
    allCellCentroids{simID} = simCellCentroids;
    
end

saveFileName = fullfile(saveDir, ['Simulation_Result_Fields_' date '.mat']);
save(saveFileName, 'RMax', 'numRBins', 'allPsi4', 'allPsi6', ...
    'allPsi4Corr', 'allPsi6Corr', 'allNumCellsR', 'allGRMax', ...
    'allGRX', 'allGRY', 'allGRR', 'allCellCentroids' );

clear psiIDx kidx initPsi4 k g cellCentroids
clear bdyXY allXY virtualBdyIDx bdyBondIDx originalBdyIDx
clear psi4 psi5 dispCells incNodeIDx psi4Corr psi6Corrr
clear numCellsR GRMax GRX GRY GRAll;

%% ************************************************************************
% *************************************************************************
%               GENERATE FIGURES
% *************************************************************************
% *************************************************************************

%% Generate Violin Plots of Final Order ===================================
clear; close all; clc;

numSimK = 5; % Number of simulations for each value of k
allK = [0 2 6 20 Inf].'; % All simulated values of k
cellSize = 1; % Length scale of cells

%--------------------------------------------------------------------------
% Load Data
%--------------------------------------------------------------------------

% The working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% Load results
dataFile = fullfile(projectDir, 'Simulation_Result_Fields_23-Jun-2022.mat');
load(dataFile);

% Set up bin structure
r = linspace(0, RMax, numRBins);

%--------------------------------------------------------------------------
% Extract Initial/Final Order Parameters
%--------------------------------------------------------------------------

orderType = 4;
plotType = 'final';
if orderType == 4 
    
    allPsi = allPsi4;
    
    if strcmpi(plotType, 'Initial')
        yLim = [0 0.5];
        yTicks = 0:0.1:0.5;
    elseif strcmpi(plotType, 'Final')
        yLim = [0 0.5];
        yTicks = 0:0.1:0.5;
    else
        error('Invalid plot type');
    end
    
elseif orderType == 6
    
    allPsi = allPsi6;
    
    if strcmpi(plotType, 'Initial')
        yLim = [0.4 0.7];
        yTicks = 0.4:0.1:0.7;
    elseif strcmpi(plotType, 'Final')
        yLim = [0 0.1];
        yTicks = 0:0.05:0.1;
    else
        error('Invalid plot type');
    end
    
else
    error('Invalid order type');
end

labelStr = ['\mid \langle' char(968) '_' num2str(orderType) '\rangle \mid'];

initPsi = zeros(size(allPsi));
finalPsi = zeros(size(allPsi));
for i = 1:size(initPsi,1)
    for j = 1:size(initPsi,2)
        initPsi(i,j) = allPsi{i,j}(1,1);
        finalPsi(i,j) = allPsi{i,j}(end,1);
    end
end

initPsi = initPsi.';
finalPsi = finalPsi.';

%--------------------------------------------------------------------------
% Generate Figures
%--------------------------------------------------------------------------

% Violin plot parameters
kColors = brewermap(5, 'Set2');
showMean = true;
showMedian = false;
medianColor = [0 0 0];
showData = false;
violinAlpha = 0.6;
meanColor = [0 0 0];

fig = figure('Color', [1 1 1]);

if strcmpi(plotType, 'Initial')
    
    violinplot( abs(initPsi).', allK, 'ViolinColor', kColors, ...
        'ShowMean', showMean, 'ShowMedian', showMedian, ...
        'ShowData', showData, 'ViolinAlpha', violinAlpha, ...
        'MedianColor', medianColor );
    
elseif strcmpi(plotType, 'Final')
    
    violinplot(abs(finalPsi).', allK, 'ViolinColor', kColors, ...
        'ShowMean', showMean, 'ShowMedian', showMedian, ...
        'ShowData', showData, 'ViolinAlpha', violinAlpha, ...
        'MedianColor', medianColor );
    
else
    error('Invalid plot type');
end


ylim(yLim);
yticks(yTicks);



% ylabel(sprintf('\\mid \\langle \\psi_%d \\rangle \\mid', orderType));
ylabel(labelStr);
xlabel(sprintf('Division orientation \nconcentration \\itk'));

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


%% Generate Time Course of Orientational Order ============================
clear; close all; clc;

numSimK = 5; % Number of simulations for each value of k
allK = [0 2 6 20 Inf].'; % All simulated values of k
cellSize = 1; % Length scale of cells

%--------------------------------------------------------------------------
% Load Data
%--------------------------------------------------------------------------

% The working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% Load results
dataFile = fullfile(projectDir, 'Simulation_Result_Fields_23-Jun-2022.mat');
load(dataFile);

% Set up bin structure
r = linspace(0, RMax, numRBins);

%--------------------------------------------------------------------------
% Extract Time Course of Orientational Order Parameters
%--------------------------------------------------------------------------

orderType = 4;
if orderType == 4 
    allPsi = allPsi4;
    yLim = [0 0.45];
    yTicks = 0:0.1:0.5;
elseif orderType == 6
    allPsi = allPsi6;
    yLim = [0 0.65];
    yTicks = 0:0.1:0.7;
else
    error('Invalid order type');
end

labelStr = ['\mid \langle' char(968) '_' num2str(orderType) '\rangle \mid'];

psi = zeros(numel(allK), numel(allPsi{1}(:,1)));
psiErr = zeros(numel(allK), numel(allPsi{1}(:,1)));
for j = 1:numel(allK)
    
    curPsi = zeros(1, numel(allPsi{1}(:,1)));
    curPsiErr = zeros(1, numel(allPsi{1}(:,1)));
    
    for i = 1:numSimK
    
    curPsi = curPsi + allPsi{i,j}(:,1).';
    curPsiErr = curPsiErr + (allPsi{i,j}(:,2).').^2;
    
    end
    
    curPsi = curPsi ./ numSimK;
    curPsiErr = sqrt(curPsiErr) ./ numSimK;
    
    psi(j,:) = curPsi;
    psiErr(j,:) = curPsiErr;
    
end

psi = abs(psi);
psiErr = abs(psiErr);

%--------------------------------------------------------------------------
% Generate Figure
%--------------------------------------------------------------------------

dispKIDx = [1 2 3 4 5];
allColors = brewermap(numel(dispKIDx), 'Dark2');

fig = figure('Color', [1 1 1],  'visible', 'on');

% Change figure color axes
% set(fig, 'defaultAxesColorOrder', allColors(2:3, :));

hold on

for i = 1:numel(dispKIDx)
    shadedErrorBar( 1:numel(psi(1,:)), ...
        psi(dispKIDx(i), :), psiErr(dispKIDx(i), :), ...
        'lineprops', {'Color', allColors(i,:), 'LineWidth', 1.5} );
end

hold off

xlim([1 numel(psi(1,:))]);
ylim(yLim);

yticks(yTicks);

ylabel(labelStr)
xlabel('Time (au)');

% Re-size Figure for Paper ----------------------------------------------

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = 4.5; %5; 
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

%% Generate Psi4 Time Comparison Figure ===================================
clear; close all; clc;

numSimK = 5; % Number of simulations for each value of k
allK = [0 2 6 20 Inf].'; % All simulated values of k
cellSize = 1; % Length scale of cells

%--------------------------------------------------------------------------
% Load Data
%--------------------------------------------------------------------------

% The working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% Load results
dataFile = fullfile(projectDir, 'Simulation_Result_Fields_23-Jun-2022.mat');
load(dataFile);

% Set up bin structure
r = linspace(0, RMax, numRBins);

% Tissue geometry parameters
tissueLength = 10/cellSize;
tissueWidth = 20/cellSize;

%--------------------------------------------------------------------------
% Format Correlations For Visualization
%--------------------------------------------------------------------------

% Choose which division concentration to visualize
dispK = Inf;
dispKIDx = find(allK == dispK);

dispTIDx = [1 54 108];
psiCorr = zeros(numel(dispTIDx), numel(r));

psi4_T1 = []; % zeros(1, numel(r));
psi4_T54 = []; % zeros(1, numel(r));
psi4_T108 = []; % zeros(1, numel(r));

numCells_T1 = [];
numCells_T54 = [];
numCells_T108 = [];

for i = 1:numSimK
    
    psi4_T1 = [psi4_T1; allPsi4Corr{i, dispKIDx}(1:5, :)];
    numCells_T1 = [numCells_T1; allNumCellsR{i, dispKIDx}(1:5, :)];
    % badR = sum((allNumCellsR(10:20, :) < 50), 1) > 0;
    % psi4_T1(badR) = 0;
    
    psi4_T54 = [psi4_T54; allPsi4Corr{i, dispKIDx}(52:56, :)];
    numCells_T54 = [numCells_T54; allNumCellsR{i, dispKIDx}(52:56, :)];
    % badR = sum((allNumCellsR(51:55, :) < 50), 1) > 0;
    % psi4_T54(badR) = 0;
    
    psi4_T108 = [psi4_T108; allPsi4Corr{i, dispKIDx}(105:108, :)];
    numCells_T108 = [numCells_T108; allNumCellsR{i, dispKIDx}(105:108, :)];
    % badR = sum((allNumCellsR(99:103, :) < 50), 1) > 0;
    % psi4_T101(badR) = 0;
    
end

psi4_T1 = mean(psi4_T1, 1);
psi4_T54 = mean(psi4_T54, 1);
psi4_T108 = mean(psi4_T108, 1);

numCells_T1 = mean(numCells_T1, 1);
numCells_T54 = mean(numCells_T54, 1);
numCells_T108 = mean(numCells_T108, 1);

minNum = 50;
psi4_T1(numCells_T1 < minNum) = 0;
psi4_T54(numCells_T54 < minNum) = 0;
psi4_T108(numCells_T108 < minNum) = 0;

%--------------------------------------------------------------------------
% Generate Figure
%--------------------------------------------------------------------------

% Axis bounds
xLim = [0.5 RMax];
yLim = [1e-2, 0.3];

fig = figure('Color', [1 1 1],  'visible', 'on');
% set(fig, 'defaultAxesColorOrder', ...
%     [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250]);

hold on

% Plot averaged correlations
t1Plot = plot(r, psi4_T1, 'LineWidth', 1.5);
t54Plot = plot(r, psi4_T54, 'LineWidth', 1.5);
t108Plot = plot(r, psi4_T108, 'LineWidth', 1.5);

plot( tissueLength .* [1 1], [yLim(1) 1*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );
plot( tissueWidth .* [1 1], [yLim(1) 1*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );

% Plot algebraic decay guide line
plot(r, 0.175 .* r.^(-1/4), '--k', 'LineWidth', 1.5);

hold off

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

xlim(xLim);
ylim(yLim);

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

yLabelStr = sprintf(['4-fold correlations ' ...
    char(12296) char(968) '_4({\\itr}) '...
    char(968) '_4^*(0)' char(12297) ]);

ylabel(yLabelStr, ...
    'FontWeight', 'bold', 'FontSize', 6);
xlabel('Distance (Mean cell lengths)', ...
    'FontWeight', 'bold', 'FontSize', 7);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = (4/3)*5;
fig.Position(4) = 5; % ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

%% Generate Psi6 Time Comparison Figure ===================================
clear; close all; clc;

numSimK = 5; % Number of simulations for each value of k
allK = [0 2 6 20 Inf].'; % All simulated values of k
cellSize = 1; % Length scale of cells

%--------------------------------------------------------------------------
% Load Data
%--------------------------------------------------------------------------

% The working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% Load results
dataFile = fullfile(projectDir, 'Simulation_Result_Fields_23-Jun-2022.mat');
load(dataFile);

% Set up bin structure
r = linspace(0, RMax, numRBins);

% Tissue geometry parameters
tissueLength = 10/cellSize;
tissueWidth = 20/cellSize;

%--------------------------------------------------------------------------
% Format Correlations For Visualization
%--------------------------------------------------------------------------

% Choose which division concentration to visualize
dispK = Inf;
dispKIDx = find(allK == dispK);

dispTIDx = [1 54 108];
psiCorr = zeros(numel(dispTIDx), numel(r));

psi6_T1 = []; % zeros(1, numel(r));
psi6_T54 = []; % zeros(1, numel(r));
psi6_T108 = []; % zeros(1, numel(r));

numCells_T1 = [];
numCells_T54 = [];
numCells_T108 = [];

for i = 1:numSimK
    
    psi6_T1 = [psi6_T1; allPsi6Corr{i, dispKIDx}(1:5, :)];
    numCells_T1 = [numCells_T1; allNumCellsR{i, dispKIDx}(1:5, :)];
    % badR = sum((allNumCellsR(10:20, :) < 50), 1) > 0;
    % psi4_T1(badR) = 0;
    
    psi6_T54 = [psi6_T54; allPsi6Corr{i, dispKIDx}(52:56, :)];
    numCells_T54 = [numCells_T54; allNumCellsR{i, dispKIDx}(52:56, :)];
    % badR = sum((allNumCellsR(51:55, :) < 50), 1) > 0;
    % psi4_T54(badR) = 0;
    
    psi6_T108 = [psi6_T108; allPsi6Corr{i, dispKIDx}(105:108, :)];
    numCells_T108 = [numCells_T108; allNumCellsR{i, dispKIDx}(105:108, :)];
    % badR = sum((allNumCellsR(99:103, :) < 50), 1) > 0;
    % psi4_T101(badR) = 0;
    
end

psi6_T1 = mean(psi6_T1, 1);
psi6_T54 = mean(psi6_T54, 1);
psi6_T108 = mean(psi6_T108, 1);

numCells_T1 = mean(numCells_T1, 1);
numCells_T54 = mean(numCells_T54, 1);
numCells_T108 = mean(numCells_T108, 1);

minNum = 50;
psi6_T1(numCells_T1 < minNum) = 0;
psi6_T54(numCells_T54 < minNum) = 0;
psi6_T108(numCells_T108 < minNum) = 0;

%--------------------------------------------------------------------------
% Generate Figure
%--------------------------------------------------------------------------

% Axis bounds
xLim = [0.5 RMax];
yLim = [5e-3, 0.5];

fig = figure('Color', [1 1 1],  'visible', 'on');
% set(fig, 'defaultAxesColorOrder', ...
%     [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250]);

hold on

% Plot averaged correlations
t1Plot = plot(r, psi6_T1, 'LineWidth', 1.5);
t54Plot = plot(r, psi6_T54, 'LineWidth', 1.5);
t108Plot = plot(r, psi6_T108, 'LineWidth', 1.5);

plot( tissueLength .* [1 1], [yLim(1) 1*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );
plot( tissueWidth .* [1 1], [yLim(1) 1*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );

% Plot algebraic decay guide line
plot(r, 0.275 .* r.^(-1/4), '--k', 'LineWidth', 1.5);

hold off

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

xlim(xLim);
ylim(yLim);

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

yLabelStr = sprintf(['6-fold correlations\n' ...
    char(12296) char(968) '_6({\\itr}) '...
    char(968) '_6^*(0)' char(12297) ]);

ylabel(yLabelStr, ...
    'FontWeight', 'bold', 'FontSize', 6);
xlabel('Distance (Mean cell lengths)', ...
    'FontWeight', 'bold', 'FontSize', 7);

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = (4/3)*5;
fig.Position(4) = 5; %ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

%% Positional Correlations Comparison Figure ==============================
clear; close all; clc;

numSimK = 5; % Number of simulations for each value of k
allK = [0 2 6 20 Inf].'; % All simulated values of k
cellSize = 1; % Length scale of cells

%--------------------------------------------------------------------------
% Load Data
%--------------------------------------------------------------------------

% The working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% Load results
dataFile = fullfile(projectDir, 'Simulation_Result_Fields_23-Jun-2022.mat');
load(dataFile);

% Set up bin structure
r = linspace(0, RMax, numRBins);

% Tissue geometry parameters
tissueLength = 10/cellSize;
tissueWidth = 20/cellSize;

%--------------------------------------------------------------------------
% Format Correlations For Visualization
%--------------------------------------------------------------------------

% Choose which division concentration to visualize
dispK = Inf;
dispKIDx = find(allK == dispK);

plotDir = 'R';
if strcmpi(plotDir, 'DV')
    
    allPlotGR = allGRX;
    
    % Y-axis label
    yStr = sprintf('DV pair correlations\n g_{DV}(r)-1');
    
    % Axis bounds
    xLim = [0.5 RMax];
    yLim = [5e-3 5];
    
    % Coefficient of algebraic decay guideline
    algCoeff = 0.075;
    plotAlg = false;
    
elseif strcmpi(plotDir, 'AP')
    
    allPlotGR = allGRY;
    
    % Y-axis label
    yStr = sprintf('AP pair correlations\n g_{AP}(r)-1');
    
    % Axis bounds
    xLim = [0.5 RMax];
    yLim = [5e-3 5];
    
    algCoeff = 0.075;
    plotAlg = false;
    
elseif strcmpi(plotDir, 'R')
    
    allPlotGR = allGRR;
    
    % Y-axis label
    yStr = sprintf('Pair correlations\n g(r)-1');
    
    % Axis bounds
    xLim = [0.5 RMax];
    yLim = [1e-2 1];
    
    % Coefficient of algebraic decay guideline
    algCoeff = 0.4;
    plotAlg = true;
    
else
    
    error('Invalid correlation type');
    
end

plotGR_T1 = [];
plotGR_T54 = [];
plotGR_T108 = [];

numCells_T1 = [];
numCells_T54 = [];
numCells_T108 = [];

for i = 1:numSimK
    
    plotGR_T1 = [plotGR_T1; allPlotGR{i, dispKIDx}(1:5, :)];
    numCells_T1 = [numCells_T1; allNumCellsR{i, dispKIDx}(1:5, :)];
    
    plotGR_T54 = [plotGR_T54; allPlotGR{i, dispKIDx}(52:56, :)];
    numCells_T54 = [numCells_T54; allNumCellsR{i, dispKIDx}(52:56, :)];
    
    plotGR_T108 = [plotGR_T108; allPlotGR{i, dispKIDx}(105:108, :)];
    numCells_T108 = [numCells_T108; allNumCellsR{i, dispKIDx}(105:108, :)];
    
end

plotGR_T1 = mean(plotGR_T1, 1);
plotGR_T54 = mean(plotGR_T54, 1);
plotGR_T108 = mean(plotGR_T108, 1);

numCells_T1 = mean(numCells_T1, 1);
numCells_T54 = mean(numCells_T54, 1);
numCells_T108 = mean(numCells_T108, 1);

minNum = 50;
plotGR_T1(numCells_T1 < minNum) = 0;
plotGR_T54(numCells_T54 < minNum) = 0;
plotGR_T108(numCells_T108 < minNum) = 0;

plotGR_T1 = plotGR_T1-1; plotGR_T1(plotGR_T1 <= 0) = 1e-10;
plotGR_T54 = plotGR_T54-1; plotGR_T54(plotGR_T54 <= 0) = 1e-10;
plotGR_T108 = plotGR_T108-1; plotGR_T108(plotGR_T108 <= 0) = 1e-10;

%--------------------------------------------------------------------------
% Generate Figure
%--------------------------------------------------------------------------

fig = figure('Color', [1 1 1],  'visible', 'on');
% set(fig, 'defaultAxesColorOrder', ...
%     [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250]);

hold on

% Plot averaged correlations
t1Plot = loglog(r, plotGR_T1, 'LineWidth', 1.5);
t54Plot = plot(r, plotGR_T54, 'LineWidth', 1.5);
t108Plot = plot(r, plotGR_T108, 'LineWidth', 1.5);

plot( tissueLength .* [1 1], [yLim(1) 1*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );
plot( tissueWidth .* [1 1], [yLim(1) 1*yLim(2)], ...
    ':k', 'LineWidth', 1.5 );

% Plot algebraic decay guide line
if plotAlg
    plot(r, algCoeff .* r.^(-1/3), '--k', 'LineWidth', 1.5);
end

hold off

xlim(xLim);
ylim(yLim);

set(gca, 'FontSize', 6);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1);

% grid on

ylh = ylabel(yStr, 'FontWeight', 'bold', 'FontSize', 6);
xlh = xlabel(sprintf('Distance (Mean cell lengths)'), ...
    'FontWeight', 'bold', 'FontSize', 7);

% legend( { ...
%     sprintf('\\boldmath{$T = %0.1f$} {\\bf hours}', timePoints(1)/12), ...
%     sprintf('\\boldmath{$T = %0.1f$} {\\bf hours}', timePoints(161)/12), ...
%     sprintf('\\boldmath{$T = %0.1f$} {\\bf hours}', timePoints(234)/12), ...
%     }, 'Interpreter', 'latex', 'FontSize', 20 );

% legend( { ...
%     sprintf('{\\it T} = %0.1f hours', timePoints(15)/12), ...
%     sprintf('{\\it T} = %0.1f hours', timePoints(161)/12), ...
%     sprintf('{\\it T} = %0.1f hours', timePoints(234)/12), ...
%     }, 'FontSize', 6, 'LineWidth', 3, 'EdgeColor', 'none');

% annotation( 'textarrow', [0.5 0.6], [0.805 0.705], ...
%     'String', '\sim {\itr}^{-1/4}', ...
%     'FontWeight', 'bold', 'FontSize', 6, 'LineWidth', 1, ...
%     'HeadLength', 5, 'HeadWidth', 5 );

% Re-scale Figure to Publication Size
set(fig, 'Units', 'centimeters');

ratio = fig.Position(4) ./ fig.Position(3);
fig.Position(3) = (4/3)*5;
fig.Position(4) = 5; %ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

% ylhPos = ylh.Position
% xlhPos = xlh.Position

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

% ylhPos
% xlhPos

% ylh.Position = ylhPos;
% xlh.Position = xlhPos;
