%% Hexatic Lattice Parameter Sweep ======================================
%
%   This is a script to determine the effect of division axis angles on the
%   orientational order of an epithelial tissue. Beginning with a cell
%   lattice with modest hexatic order, cells are divided in a wave through
%   a user specified number of divisions. Cell division axes are drawn from
%   a von Mises distribution. The concentration of the von Mises
%   distribution is the parameter being swept
%
%   by Dillon Cislo 2022/06/16
%========================================================================


% Add relevant file structure to MATLAB path
[DDMDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
DDMDir = [DDMDir '/..'];
addpath(genpath(DDMDir));
clear DDMDir

%% Set Parameter Sweep Metadata =========================================
clear; close all; clc;

% Set Division Distribution Parameters ----------------------------------
allK = [0 2 6 20 Inf];

% Set Directory Data ----------------------------------------------------

% The working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
% projectDir = pwd; % For use with headless mode
cd(projectDir);

% Directory that will hold the output results
% saveDir = fullfile(projectDir, ... % Use with Dryad data
%     'Noisy_Hexatic_SNR18_Param_Sweep_Results_20220618');
saveDir = fullfile(projectDir, ... % Make the data yourself
    ['Noisy_Hexatic_SNR18_Param_Sweep_Results_' date]);
if ~exist(saveDir,  'dir'), mkdir(saveDir); end

%% RUN PARAMETER SWEEP **************************************************
% ***********************************************************************

% The number of simulations to run for each value of k
numSimK = 5;

% The total number of simulations to run
numSim = numSimK * numel(allK);

% randSeeds = (1:5:(5*numSim))+2;
randSeeds = (1:5:(5*numSim))+3;
assert(numel(randSeeds) == numSim, ...
    'Number of seeds must equal the number of total simulations');

plotResults = false;
verbose = false;
parfor simID = 1:numSim

    warning('off', 'all');
    
    rng(randSeeds(simID), 'twister');

    % Convert linear simulation index to (simK, K) format
    [kSimID, kIDx] = ind2sub([numSimK, numel(allK)], simID);
    fprintf('Running simulation %d for k = %d\n', kSimID, allK(kIDx));

    saveFileName = fullfile(saveDir, ...
        sprintf('Simulation_Result_K%d_N%d.mat', allK(kIDx), kSimID));
    if exist(saveFileName, 'file'), continue; end

    %====================================================================
    % GENERATE INITIAL LATTICE
    %====================================================================

    cellSize = 1; % Length scale for cells

    % Points are generated inside a large box. Only points within a smaller
    % box, itself contained within the large box, will be included initial
    % cell lattice
    W = 2*10*cellSize; % Small bounding box width
    H = 9*cellSize; %2*4*cellSize; % Small bounding box height
    WW = W + 10*cellSize; % Big bounding box width
    HH = H + 10*cellSize; % Big bounding box height

    latticeOptions = struct();
    latticeOptions.xLim = WW * [-1 1] / 2;
    latticeOptions.yLim = HH * [-1 1] / 2;
    latticeOptions.pointSetType = 'equilateral';
    latticeOptions.SNR = 18;
    latticeOptions.sideLength = cellSize;
    latticeOptions.orientation = 'horizontal';

    % Generate intial equilateral point set
    allXY = CELTIGS.generate_points( latticeOptions );

    % Determine which points are included in the small box
    incBoxX = (W+1.0*cellSize) * [-1 1] / 2;
    incBoxY = (H+1.0*cellSize) * [-1 1] / 2;
    incNodeIDx = ...
        (incBoxX(1) <= allXY(:,1)) & (allXY(:,1) <= incBoxX(2)) ...
        & (incBoxY(1) <= allXY(:,2)) & (allXY(:,2) <= incBoxY(2));

    if verbose
        fprintf('N = %d interior cells generated\n', sum(incNodeIDx));
    end

    % Re-set the bounding box stored in the lattice options
    latticeOptions.xLim = W * [-1 1] / 2;
    latticeOptions.yLim = H * [-1 1] / 2;

    %--------------------------------------------------------------------
    % Establish Cell Row Identity
    %--------------------------------------------------------------------
    % WAIT UNTIL AFTER THIS SECTION TO ADD NOISE TO THE CELL POSITIONS!!!

    % Only establish row identities for included cells
    XY = allXY(incNodeIDx, :);

    % Estimate the number of rows
    % [rowIDx, rowY] = kmeans(XY(:,2), ceil(H/cellSize)+1);

    klist = 2:ceil(2*H/cellSize); % The number of clusters you want to try
    myfunc = @(X,K)(kmeans(X, K));
    eva = evalclusters(XY(:,2), myfunc, 'CalinskiHarabasz', 'klist', klist);

    numMeans = eva.OptimalK;
    % numMeans = 5;
    [rowIDx, rowY]  = kmeans(XY(:,2), numMeans);

    % Sort the rows by Y-position
    [~, sortIDx] = sort(rowY, 'descend');
    [~, newRowIDx] = ismember((1:numel(rowY)).', sortIDx);
    rowIDx = changem(rowIDx, newRowIDx, (1:numel(rowY)).');

    assert( max(rowIDx(:)) >= 11, ...
        ['Invalid number of rows generated in simulation ', ...
        '%d for k = %d\n'], kSimID, allK(kIDx));

    % clear klist myfunc eva sortIDx newRowIDx numMeans

    % Generate a Division Order From Row IDs ----------------------------
    divOrderFromRows = true;
    if divOrderFromRows

        [~, divOrder] = sortrows([rowIDx, abs(XY(:,1))]);
        [~, divOrder] = ismember((1:numel(rowIDx)).', divOrder);
        [~, divOrder] = sort(divOrder);

        incRowIDx = 4:8;
        divOrder(~ismember(rowIDx(divOrder), incRowIDx)) = [];

    end

    %--------------------------------------------------------------------
    % Calculate Voronoi Diagram Parameters
    %--------------------------------------------------------------------

    % OPTIONAL: Add noise at the end following lattice creation
    allXY = awgn(allXY, latticeOptions.SNR);

    % Determine Voronoi connectivity of current configuration
    delTri = delaunayTriangulation(allXY(:,1), allXY(:,2));
    [v0, c0] = voronoiDiagram(delTri);

    % Extract boundary cells
    bdyCells = cellfun( @(x) ismember(1, x), c0 );
    assert(~any(incNodeIDx & bdyCells), 'Invalid boundary produced');

    % Convert the voronoi cell connectivity list to a (NaN-padded)
    % matrix Useful for plotting purposes
    maxFaceSize = max(cellfun(@(x) numel(x), c0));
    voronoiFace = nan(size(c0,1), maxFaceSize);
    for i = 1:size(c0,1)
        voronoiFace(i, 1:numel(c0{i})) = c0{i};
    end

    % Prune Vertex List -------------------------------------------------

    v = v0; c = c0;
    c(~incNodeIDx) = []; % Remove exterior cells
    XY = allXY(incNodeIDx, :);

    incVIDx = unique([c{:}]); % The Voronoi vertices remaining
    oldVIDx = (1:size(v,1)).'; % Original vertex IDs
    rmVIDx = ~ismember(oldVIDx, incVIDx); % IDs of unused vertices
    v(rmVIDx, :) = []; % Prune vertex list

    % Prune individual cell vertex lists and renumber vertices
    oldVIDx(rmVIDx) = [];
    c = cellfun( @(x) changem(x, (1:size(v,1)).', oldVIDx), ...
        c, 'Uni', false);

    assert(all(ismember((1:size(v,1)).', unique([c{:}]))), ...
        'Voronoi vertex inclusion mismatch');

    %--------------------------------------------------------------------
    % Report Cell Geometry
    %--------------------------------------------------------------------

    % Find the Voronoi area/perimeter of each cell
    cellAreas = zeros(numel(c), 1);
    cellPerim = zeros(numel(c), 1);
    for cid = 1:numel(c)
        cellAreas(cid) = ...
            polyarea(v(c{cid}, 1), v(c{cid}, 2));
        cellPerim(cid) = ...
            perimeter(polyshape(v(c{cid}, 1), v(c{cid}, 2)));
    end

    if verbose

        fprintf('Mean cell area A = %0.5f\n', mean(cellAreas));
        fprintf('Mean cell perimeter P = %0.5f\n', mean(cellPerim));

    end

    % clear maxFaceSize voronoiFace
    % clear delTri bdyCells incBoxX incBoxY
    % clear xLim yLim

    %====================================================================
    % SET SIMULATION PARAMETERS
    %====================================================================

    % Generate 'mechanicalFeedback' style lattice structure
    if (size(v,2) ~= 3), v = [v zeros(size(v,1), 1)]; end
    g = GLattConversion2(c, v, false);

    % Set target area for each cell
    % g.A0 = cellSize.^2 * ones([numel(g.cells) 1]); % Square lattice
    g.A0 = mean(cellAreas) * ones([numel(g.cells) 1]);

    % Prefactor for area elasticity energy;
    g.kA0 = ones([numel(g.cells) 1]);

    % Set target perimeter for each cell
    % g.p0 = 4 * cellSize * ones([numel(g.cells) 1]); % Square lattice
    g.p0 = 4 * sqrt(g.A0);

    % Prefactor for perimeter tension energy
    g.pT0 = ones([numel(g.cells) 1]);

    % Set target edge lengths to be the current bond lengths...
    g.l0 = g.verts(g.bonds(:,2), :)-g.verts(g.bonds(:,1), :);
    g.l0 = sqrt(sum((g.l0).^2, 2));

    % Define edge tensions
    g.T0 = ones([size(g.bonds,1) 1]);
    g.T = g.T0;

    % Define Clonal Identities
    g.clones = nan([numel(g.cells) 1]);
    clsubset = [];
    g.clones(clsubset) = 1;

    % clear clsubset

    % Simulation Parameters ('g.param') -----------------------------------
    % g.param(1): Prefactor for cell area energy
    %       - p1 * kA0 * (A-A0).^2 / 2 OR
    %       - p1 * A / A0
    %
    % g.param(2): Prefactor for cell perimeter tension term
    %       - p2 * pT0 * p.^(1+p6)
    %
    % g.param(3): Prefactor for cell perimeter elasticity term
    %       - p3 * (p-p0).^2 / 2
    %
    % g.param(4): Prefactor for bond length tension term
    %       - p4 * T0 * l^(1+p6)
    %
    % g.param(5): Prefctor for bond length elasticity term
    %       - p5 * (l-l0)^2 / 2
    %
    % g.param(6): Used as exponent in tension terms
    %
    % ONLY CALCULATES Earea if param(1) > 0
    % ONLY CALCULATES Eperim if (param(2) > 0 || param(3) > 0)
    % ONLY CALCULATES Ebond if (param(4) > 0 || param(5) > 0)

    g.param = [1 0 1 0 0 0];

    % View Results ------------------------------------------------------

    if plotResults

        ect = sumTension(g);

        divOrderCT = nan(numel(g.cells), 1);
        [isDiv, divID] = ismember((1:numel(g.cells)).', divOrder);
        divOrderCT(isDiv) = divID(isDiv);
        divOrderCM = [parula(numel(divOrder)); 0 0 0];
        divOrderCT(isnan(divOrderCT)) = size(divOrderCM,1);

        LatticePresentation2(g, struct('cellIndex', false, ...
            'colorTable', divOrderCT, 'cellColormap', divOrderCM, ...
            'edgeColor', ect, 'edgeColorAbsolute', true, ...
            'includeFaceColorbar', true, ...
            'cellColorRange', [1 size(divOrderCM,1)]));

        % LatticePresentation2(g, struct('cellIndex', false, ...
        %     'colorTable', g.row, ...
        %     'cellColormap', parula(numel(unique(g.row))), ...
        %     'edgeColor', ect, 'edgeColorAbsolute', true, ...
        %     'includeFaceColorbar', true));

        vcounts = histcounts([c{:}], (1:(size(v,1)+1))-0.5);
        allSingleV = find(vcounts == 1);
        hold on
        scatter(v(allSingleV, 1), v(allSingleV, 2), 'filled', 'r');
        hold off

        title('Initial Lattice');

        clear ect vcounts allSingleV
        clear divOrderCT isDiv divID divOrderCM

        pause(5);
        close all;

    end

    %====================================================================
    % RELAX INITIAL LATTICE MECHANICALLY
    %====================================================================

    % Make a copy of the unrelaxed lattice
    g0 = g;

    maxiter = 5000;
    g = relaxEpithelium2D(g, maxiter);

    if plotResults

        ect = sumTension(g);
        ct = g.stress;
        ct(isnan(ct))=0;

        LatticePresentation2(g, struct('transparent', false, ...
            'cellIndex', false, 'includeFaceColorbar', true, ...
            'edgeColor', ect, 'colorTable', ct, 'lineWidth', 2))

        title('Initial Relaxed Lattice')

        clear ect ct

        pause(5);
        close all;

    end

    %====================================================================
    % SIMULATE DIVISION WAVES
    %====================================================================

    % The number of division waves
    numDivWaves = 1;
    assert(numDivWaves >= 1, 'Please simulate at least one wave');
    cell_cmap = brewermap(numDivWaves+1, 'Dark2');
    cell_crange = [0 numDivWaves];

    % T1 parameters
    T1_Lmin = 0; %0.05;
    allowAngleT1s = true;
    T1_AngMin = deg2rad(20);
    T1_freq = 25;

    % Edge merger parameters
    allowEdgeMergers = true;
    minEdgeAngle = T1_AngMin; % deg2rad(20);
    edgeMergeFreq = 25;

    % Whether or not to enforce horizontal confinement
    enforceConfinement = true;
    enforceConfinementRows = false;

    % Determine growth parameters (ASSUMES A0 IS UNIFORM FOR ALL CELLS)
    A0 = mean(g.A0); % The basic intrinsic area for cells
    P0 = mean(g.p0); % The basic intrinsic perimeter for cells
    Adiv = 1.5 * A0; % The intrinsic area at which cells will divide
    Pdiv = 4 * sqrt(Adiv); % sqrt(2) * P0;
    numGrowthSteps = 50;% The number of time steps over which the cell is grown

    % Cell growth rate (implicit assumption of unit time step dt = 1)
    growthRate = (Adiv-A0)/numGrowthSteps;
    perimGrowthRate = (Pdiv-P0)/numGrowthSteps;

    % The time point IDs for when each division wave is completed
    % divTimes = ...
    %     ((numGrowthSteps+1) * numel(g.cells) .* ...
    %     (2.^(0:numDivWaves)-1)).' + 1;

    % Division axis options
    divAxisOptions = struct();
    divAxisOptions.mu = pi/2;
    divAxisOptions.k = allK(kIDx);
    divAxisOptions.n = 2;
    divAxisOptions.truncateDomain = true;
    divAxisOptions.domain = 'Pi';
    divAxisOptions.plotSamples = false;

    assert(divAxisOptions.k >= 0, ...
        'Distribution concentration must be non-negative');

    % Add a field to reconstruct cell lineage tracknig
    g.parentID = zeros(numel(g.cells), 1);

    % Add a field to count divisions within each lineage
    g.numDiv = zeros(numel(g.cells), 1);

    % Store the cell lattice structures following division events
    allGLattices = {g};

    % A time point counter across all waves
    tidx = 2;

    if plotResults

        % Visualize initial configuration
        fig = figure('units', 'normalized', ...
            'outerposition', [0.5 0 0.5 1],  'color', [1 1 1]);

        L = 1.1 * max(abs(g.verts(:))); % Axis limits
        LatticePresentation2(g, struct('transparent', false, ...
            'cellIndex', false, ...
            'colorTable', g.numDiv, 'cellColormap', cell_cmap, ...
            'cellColorRange', cell_crange, 'lineWidth', 2));
        axis([-L L -L L]);
        % axis off
        drawnow

    end

    for nwIDx = 1:numDivWaves

        if verbose
            fprintf('Simulating Division Wave %d/%d\n', ...
                nwIDx, numDivWaves);
        end

        %----------------------------------------------------------------
        % Reset Division Order
        %----------------------------------------------------------------
        if divOrderFromRows

            if (nwIDx > 1)

                newDivOrder = numel(allGLattices{1}.cells) + ...
                    (1:numel(originalDivOrder)).';
                newDivOrder = repmat(newDivOrder, 1, nwIDx-1) + ...
                    repmat(numel(originalDivOrder) * ((1:(nwIDx-1))-1), ...
                    numel(originalDivOrder), 1);
                newDivOrder = [originalDivOrder, newDivOrder];

                originalRowIDx = g.row(originalDivOrder);
                [~,~,ix] = unique(originalRowIDx);
                numInRow = accumarray(ix,1);

                newDivOrder = mat2cell(newDivOrder, ...
                    numInRow, size(newDivOrder, 2));
                newDivOrder = cellfun(@(x) x(:), newDivOrder, 'Uni', false);
                newDivOrder = vertcat(newDivOrder{:});
                divOrder = newDivOrder;

                % clear newDivOrder originalRowIDx ix numInRow

            end

        else

            % Extract cell centroids
            cellCentroids = zeros(numel(g.cells), 2);
            for i = 1:numel(g.cells)
                [cellCentroids(i,1), cellCentroids(i,2)] = ...
                    centroid( polyshape( ...
                    g.verts(g.bonds(g.cells{i},1),1), ...
                    g.verts(g.bonds(g.cells{i},1),2) ) );
            end

            % Origin lies at the top of the patch
            wavefrontOrigin = ...
                [ min(cellCentroids(:,1)) + range(cellCentroids(:,1))/2, ...
                max(cellCentroids(:,2)) ];

            % Origin lies at center of the patch
            % wavefrontOrigin = ...
            %     [ min(cellCentroids(:,1)) + range(cellCentroids(:,1))/2, ...
            %     min(cellCentroids(:,2)) + range(cellCentroids(:,2))/2 ];

            divOrderOptions = struct();
            divOrderOptions.plotOrder = false;
            divOrderOptions.wavefrontOrigin = wavefrontOrigin;
            divOrderOptions.wavefrontVelocity = [15 1];

            divOrder = CELTIGS.determine_division_order( ...
                cellCentroids, divOrderOptions);
            [~, divOrder] = sort(divOrder);

        end

        %----------------------------------------------------------------
        % Generate Division Waves
        %----------------------------------------------------------------

        for didx = 1:numel(divOrder)
    
            if verbose
                try
                    progressbar(didx, numel(divOrder));
                catch
                    fprintf('Simulating the %d/%d for current wave\n', ...
                        didx, numel(divOrder));
                end
            end

            % The ID of the cell about to divide
            divID = divOrder(didx);

            % Choose the division axis
            if isinf(divAxisOptions.k)
                divAxis = divAxisOptions.mu;
            elseif (divAxisOptions.k == 0)
                divAxis = 2 * pi * rand(1);
            else
                divAxis = randAngleVonMises( 1, ...
                    'Mu', divAxisOptions.mu, ...
                    'Kappa', divAxisOptions.k, ...
                    'numPeaks', divAxisOptions.n, ...
                    'Domain', divAxisOptions.domain, ...
                    'TruncateDomain', divAxisOptions.truncateDomain, ...
                    'PlotSamples', divAxisOptions.plotSamples );
            end

            divAxis = [cos(divAxis), sin(divAxis)];

            % Run Growth Steps ------------------------------------------

            for i = 1:numGrowthSteps

                % Perform T1 transitions
                if (mod(tidx, T1_freq) == 0)
                    % g = removeDegenerateEdges(g);
                    t1idx = decideT1(g, T1_Lmin);
                    if allowAngleT1s
                        t1idx = cat(1, t1idx, ...
                            decideT1BoundaryAngle(g, T1_AngMin));
                    end
                    if ~isempty(t1idx)
                        % disp(['Performing T1 transitions for ' ...
                        %     num2str(t1idx.')]);
                        for j = 1:numel(t1idx)
                            [g, newBondIDx, oldBondIDx] = ...
                                T1transition(g, t1idx(j));
                            if ~isempty(newBondIDx)
                                t1idx = changem(t1idx, ...
                                    newBondIDx, oldBondIDx);
                            end
                        end
                        % g = removeDegenerateEdges(g);
                    end
                end

                % Perform degenerate edge mergers
                if (mod(tidx, edgeMergeFreq) == 0) && allowEdgeMergers
                    % g = removeDegenerateEdges(g);
                    g = mergeDegenerateBoundaryEdges(g, minEdgeAngle);
                    % g = removeDegenerateEdges(g);
                end

                % Increase target area of dividing cell
                g.A0(divID) = g.A0(divID) + growthRate; % *(1);
                g.p0(divID) = g.p0(divID) + perimGrowthRate; % *(1);

                % OPTIONAL CONFINEMENT STEP
                if enforceConfinementRows
                    uniqueRows = unique(g.row);
                    for ridx = 1:numel(uniqueRows)
                        maxXVal = maxXRows(ridx) + cellSize/2;
                        minXVal = -maxXVal;
                        curRCIDx = g.row == uniqueRows(ridx);
                        curRBIDx = unique([g.cells{curRCIDx}]);
                        curRVIDx = unique([g.bonds(curRBIDx,1); ...
                            g.bonds(curRBIDx,2)]);
                        curRVIDx = ismember((1:size(g.verts,1)).', curRVIDx);
                        g.verts((g.verts(:,1) > maxXVal) & curRVIDx, 1) = maxXVal;
                        g.verts((g.verts(:,1) < minXVal) & curRVIDx, 1) = minXVal;
                    end
                    %clear uniqueRows ridx maxXVal minXVal
                    %clear curRCIDx curRBIDx curRVIDx
                elseif enforceConfinement
                    maxXVal = latticeOptions.xLim(2) + cellSize/2;
                    minXVal = latticeOptions.xLim(1) - cellSize/2;
                    g.verts(g.verts(:,1) > maxXVal, 1) = maxXVal;
                    g.verts(g.verts(:,1) < minXVal, 1) = minXVal;
                end

                % Relax tissue mechanically
                g = relaxEpithelium2D(g, maxiter);

                % Visualize deformation
                if plotResults

                    clf
                    L = 1.1 * max(abs(g.verts(:))); % Axis limits
                    % visualizePressure(g, struct('L', L));
                    LatticePresentation2(g, ...
                        struct('transparent', false, ...
                        'colorTable', g.numDiv, ...
                        'cellColormap', cell_cmap, ...
                        'cellColorRange', cell_crange, ...
                        'lineWidth', 2, 'cellIndex', false));
                    axis([-L L -L L]);
                    % axis off
                    drawnow

                end

                tidx = tidx+1;

            end

            % Divide Cell -----------------------------------------------

            intersectMethod = 'lineIntersect';
            g = divideCell(g, divID, divAxis, intersectMethod);

            % Manually set all target areas/perimeters back to default
            % (This shouldn't be necessary, but I'm careful)
            g.A0 = A0 * ones([numel(g.cells) 1]);
            g.p0 = P0 * ones([numel(g.cells) 1]);

            % Update division count
            g.numDiv(divID) = g.numDiv(divID)+1;
            g.numDiv = [g.numDiv; g.numDiv(divID)];

            % Update lineage IDs
            g.parentID = [(1:numel(g.cells)).'; divID];

            % Perform T1 transitions
            if (mod(tidx, T1_freq) == 0)
                % g = removeDegenerateEdges(g);
                t1idx = decideT1(g, T1_Lmin);
                if allowAngleT1s
                    t1idx = cat(1, t1idx, ...
                        decideT1BoundaryAngle(g, T1_AngMin));
                end
                if ~isempty(t1idx)
                    % disp(['Performing T1 transitions for ' ...
                    %     num2str(t1idx.')]);
                    for j = 1:numel(t1idx)
                        [g, newBondIDx, oldBondIDx] = ...
                            T1transition(g, t1idx(j));
                        if ~isempty(newBondIDx)
                            t1idx = changem(t1idx, ...
                                newBondIDx, oldBondIDx);
                        end
                    end
                    % g = removeDegenerate Edges(g);
                end
            end

            % Perform degenerate edge mergers
            if (mod(tidx, edgeMergeFreq) == 0) && allowEdgeMergers
                % g = removeDegenerateEdges(g);
                g = mergeDegenerateBoundaryEdges(g, minEdgeAngle);
                % g = removeDegenerateEdges(g);
            end

            % OPTIONAL CONFINEMENT STEP
            if enforceConfinementRows
                uniqueRows = unique(g.row);
                for ridx = 1:numel(uniqueRows)
                    maxXVal = maxXRows(ridx) + cellSize/2;
                    minXVal = -maxXVal;
                    curRCIDx = g.row == uniqueRows(ridx);
                    curRBIDx = unique([g.cells{curRCIDx}]);
                    curRVIDx = unique([g.bonds(curRBIDx,1); ...
                        g.bonds(curRBIDx,2)]);
                    curRVIDx = ismember((1:size(g.verts,1)).', curRVIDx);
                    g.verts((g.verts(:,1) > maxXVal) & curRVIDx, 1) = maxXVal;
                    g.verts((g.verts(:,1) < minXVal) & curRVIDx, 1) = minXVal;
                end
                %clear uniqueRows ridx maxXVal minXVal
                %clear curRCIDx curRBIDx curRVIDx
            elseif enforceConfinement
                maxXVal = latticeOptions.xLim(2) + cellSize/2;
                minXVal = latticeOptions.xLim(1) - cellSize/2;
                g.verts(g.verts(:,1) > maxXVal, 1) = maxXVal;
                g.verts(g.verts(:,1) < minXVal, 1) = minXVal;
            end

            % Relax tissue mechanically
            g = relaxEpithelium2D(g, maxiter);

            allGLattices = cat(1, allGLattices, g);

            % Visualize deformation
            if plotResults

                clf
                L = 1.1 * max(abs(g.verts(:))); % Axis limits
                % visualizePressure(g, struct('L', L));
                LatticePresentation2(g, ...
                    struct('transparent', false, ...
                    'colorTable', g.numDiv, ...
                    'cellColormap', cell_cmap, ...
                    'lineWidth', 2, 'cellIndex', false, ...
                    'cellColorRange', cell_crange ));
                axis([-L L -L L]);
                % axis off
                drawnow

            end

            tidx = tidx+1;

        end

    end

    %====================================================================
    % SAVE RESULTS
    %====================================================================

    saveVarNames = {'latticeOptions', 'g0', 'allGLattices', 'divOrder'};
    saveVars = {latticeOptions, g0, allGLattices, divOrder};
    parsave(saveFileName, saveVarNames, saveVars);

    warning('on', 'all');

end


