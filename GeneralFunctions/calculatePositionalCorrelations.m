function [ GR_Max, GR_X, GR_Y, GR_All, r ] = ...
    calculatePositionalCorrelations( X, varargin )
%CALCULATEPOSITIONALCORRELATIONS Calculate the pairwise positional
%correlation function (g(r = |r1-r2|)). Correlation funciton is calculated
%by averaging over a cut in the 2D histogram of relative point placements
%
%   INPUT PARAMETERS:
%
%       - X:            #Px2 list of scattered point set (x,y) coordinates
%
%   OPTIONAL (Name, Value = Default)-PAIRS:
%
%       - ('CellSize', cellSize = []): A user supplied length scale in
%       terms of which pairwise distances may be normalized
%
%       - ('MaxR', RMax = []): The maximum distance over which to calculate
%       the correlation functions. Units of RMax depend on the value of
%       'normDists'
%
%       - ('NumRBins', numRBins = []): The number of bins with which to
%       calculate the correlations
%
%       - ('NumThetaBins', numThetaBins = 500): The number of bins used to
%       calculate the axis of symmetry of the sample
%
%       - ('NumGridBins', numGridBins = 500): The number of histogram bins
%       used to calculate the axis of symmetry of the sample
%
%       - ('SymmetryRingR', symR = []): The size of the ring over which to
%       calculate the axis of symmetry of the sample in units of cell
%       lengths
%
%       - ('CutWidth', cutWidth = []); The width of the radial cut used to
%       calculate the correlations
%
%       - ('IncludedNodes', incNodeIDx = []): #IPx1 vector of point IDs
%       indicated which nodes to include in the calculation
%
%       - ('NormalizeDistances', normDists = true): If true, distances
%       between points are normalized by the cellSize parameter
%
%       - ('Normalization', normType = 'none'): The type of normalization
%       for the distribution function
%
%       - ('CorrectFiniteSize', correctFiniteSize = false): If true, radial
%       distribution function is multiplied by a shape parameter to correct
%       for finite size effects. WARNING: This method does not produce good
%       results for the cut distributions and really only produces
%       reasonable results for the full radial distribution for disordered
%       configurations
%
%       - ('EffectiveVolume', useEffectiveVolume = false): If true, radial
%       distribution function is calculated using an effective volume
%       method to accound for finite sample size. This flag will override
%       both the 'Normalization' flag and the 'CorrectFiniteSize' flag
%
%       - ('EffectiveVolumeMethod', effectiveVolumeMethod = 'polygon'): The
%       method used to calculate the effective volume. If the boost library
%       is present on the machine the user should always select the
%       'polygon' method. The polygon method is still usable without it,
%       but is intolerably slow. The 'grid' method can be faster than the
%       MATLAB built-in 'polygon' method, but may produce erroneous results
%       if the grid size is too sparse (especially for the linear
%       distribution cuts)
%
%       - ('EffectiveVolmeGridSize', gridSize = [100 100]): The grid size
%       used to calculate the effective volume using the 'grid' method
%
%       - ('PlotResults', plotResults = false): If true, radial
%       correlations are plotted
%
%   OUTPUT PARRAMETERS:
%
%       - GR_Max:   1x#(numRBins) radial distribution function along the
%                   axis of symmetry of the system
%
%       - GR_X:     1x#(numRBins) radial distribution function along the
%                   x-axis
%
%       - GR_Y:     1x#(numRBins) radial distribution function along the
%                   y-axis
%
%       - GR_All:   1x#(numRBins) integrated radial distribution function
%                   for all theta
%
%       - r:        1x#(numBins)radial distance bin locations. Used for
%                   plotting results
%
%   by Dillon Cislo 03/02/2021

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if (nargin < 1), error('Please input point set coordinates'); end
validateattributes(X, {'numeric'}, {'2d', 'ncols', 2, 'finite', 'real'});

% Set default values
cellSize = [];
RMax = [];
numRBins = [];
numThetaBins = 500;
numGridBins = 500;
symR = 10;
cutWidth = [];
incNodeIDx = [];
normDists = true;
normType = 'none';
correctFiniteSize  = false;
useEffectiveVolume = false;
effectiveVolumeMethod = 'polygon';
gridSize = [100 100];
plotResults = false;

allNormTypes = {'none', 'probability', 'pairdist', 'radialdist'};
allEVMethods = {'grid', 'polygon'};

for i = 1:length(varargin)
    
    if isa(varargin{i}, 'double'), continue; end
    
    if isa(varargin{i}, 'logical'), continue; end
    
    if strcmpi(varargin{i}, 'CellSize')
        cellSize = varargin{i+1};
        if ~isempty(cellSize)
            validateattributes(cellSize, {'numeric'}, ...
                {'scalar', 'positive', 'finite', 'real'});
        end
    end
    
    if strcmpi(varargin{i}, 'MaxR')
        RMax = varargin{i+1};
        if ~isempty(RMax)
            validateattributes(RMax, {'numeric'}, ...
                {'scalar', 'positive', 'finite', 'real'});
        end
    end
    
    if strcmpi(varargin{i}, 'NumRBins')
        numRBins = varargin{i+1};
        if ~isempty(numRBins)
            validateattributes(numRBins, {'numeric'}, ...
                {'scalar', 'integer', 'positive', 'finite', 'real'});
        end
    end
    
    if strcmpi(varargin{i}, 'NumThetaBins')
        numThetaBins = varargin{i+1};
        validateattributes(numThetaBins, {'numeric'}, ...
            {'scalar', 'integer', 'positive', 'finite', 'real'});
    end
    
    if strcmpi(varargin{i}, 'NumGridBins')
        numGridBins = varargin{i+1};
        validateattributes(numGridBins, {'numeric'}, ...
            {'scalar', 'integer', 'positive', 'finite', 'real'});
    end
    
    if strcmpi(varargin{i}, 'SymmetryRingR')
        symR = varargin{i+1};
        validateattributes(symR, {'numeric'}, ...
            {'scalar', 'positive', 'finite', 'real'});
    end
    
    if strcmpi(varargin{i}, 'CutWidth')
        cutWidth = varargin{i+1};
        if ~isempty(cutWidth)
            validateattributes(cutWidth, {'numeric'}, ...
                {'scalar', 'positive', 'finite', 'real'});
        end
    end
    
    if strcmpi(varargin{i}, 'IncludedNodes')
        incNodeIDx = varargin{i+1};
        if ~isempty(incNodeIDx)
            validateattributes(incNodeIDx, {'numeric'}, ...
                {'vector', 'integer', 'positive', 'finite', ...
                'real', '<=', size(X,1)});
        end
    end
    
    if strcmpi(varargin{i}, 'NormalizeDistances')
        normDists = varargin{i+1};
        validateattributes(normDists, {'logical'}, {'scalar'});
    end

    if strcmpi(varargin{i}, 'Normalization')
        normType = lower(varargin{i+1});
        validateattributes(normType, {'char'}, {'vector'});
        assert(ismember(normType, allNormTypes), ...
            'Invalid normalization type');
    end

    if strcmpi(varargin{i}, 'CorrectFiniteSize')
        correctFiniteSize = varargin{i+1};
        validateattributes(correctFiniteSize, {'logical'}, {'scalar'});
    end

    if strcmpi(varargin{i}, 'EffectiveVolume')
        useEffectiveVolume = varargin{i+1};
        validateattributes(useEffectiveVolume, {'logical'}, {'scalar'});
    end
    
    if strcmpi(varargin{i}, 'PlotResults')
        plotResults = varargin{i+1};
        validateattributes(plotResults, {'logical'}, {'scalar'});
    end

    if strcmpi(varargin{i}, 'EffectiveVolumeMethod')
        effectiveVolumeMethod = lower(varargin{i+1});
        validateattributes(effectiveVolumeMethod, {'char'}, {'vector'});
        assert(ismember(effectiveVolumeMethod, allEVMethods), ...
            'Invalid effective volume method');
    end

    if strcmpi(varargin{i}, 'EffectiveVolumeGridSize')
        gridSize = varargin{i+1};
        validateattributes(gridSize, {'numeric'}, ...
                {'vector', 'integer', 'positive', 'finite', ...
                'real', 'numel', 2});
    end

end

% Override other settings if effective volume methods are selected
if useEffectiveVolume
    correctFiniteSize = false;
    normType = 'none';
end

% Generate Optional Input Parameters --------------------------------------

% Construct a Delaunay triangulation of the input point set
delTri = delaunayTriangulation(X);
[v, c] = voronoiDiagram(delTri);

% The edge list of the triangulation
% bondIDx = delTri.edges;

% Determine which cells lie on the boundary of the diagram
bdyCells = cellfun( @(x) ismember(1, x), c );

% Determine which subset of nodes to include in the calculation
if isempty(incNodeIDx), incNodeIDx = find(~bdyCells); end
N = numel(incNodeIDx);

% Calculate the area of each Voronoi polygon
cellAreas = zeros(numel(c), 1);
for i = 1:numel(c)
    if ~bdyCells(i)
        cellAreas(i) = polyarea(v(c{i},1), v(c{i},2));
    end
end

% Calculate an 'average cell size'
if isempty(cellSize)
    
    % Calculate an 'average cell size' based on pairwise distances
    % edgeLengths = X(bondIDx(:,2), :) - X(bondIDx(:,1), :);
    % edgeLengths = sqrt(sum(edgeLengths.^2, 2)); 
    % incEdges = all( ismember(bondIDx, incNodeIDx), 2 );
    % cellSize = mean(edgeLengths(incEdges));

    % Calculate an 'average cell size' based on cell areas
    cellSize = sqrt(mean(cellAreas(incNodeIDx)));
    
end

if normDists, cellAreas = cellAreas ./ cellSize.^2; end

% The uniform number density of the included region
numDensity = N / sum(cellAreas(incNodeIDx));

if useEffectiveVolume
    bdyPoly = extractBoundaryPolygon(delTri, v, c, incNodeIDx, false);
    if normDists
        bdyPoly.Vertices = bdyPoly.Vertices / cellSize;
    end
else
    bdyPoly = [];
end

%--------------------------------------------------------------------------
% Construct Pairwise Separation Vectors
%--------------------------------------------------------------------------

% Calculate the 2D separation vectors for each directed pair of nodes
delXij = zeros(N, N);
delYij = zeros(N, N);

for i = 1:N
    
    delXij(:,i) = X(incNodeIDx, 1) - X(incNodeIDx(i), 1);
    delYij(:,i) = X(incNodeIDx, 2) - X(incNodeIDx(i), 2);
    
end

% Compile the IDs of cells in each pair
[~, J] = ndgrid(1:numel(incNodeIDx), 1:numel(incNodeIDx));
J = J(:);

% Assemble into a single matrix
delRij = [ delXij(:), delYij(:) ];

% Remove self-reference entries
rmIDx = (sqrt(sum(delRij.^2, 2)) < 1e-13);
delRij(rmIDx, :) = [];
J(rmIDx, :) = [];

% Normalize distances by average cell size
if normDists, delRij = delRij ./ cellSize; end

% The magnitude of the separation vectors
R = sqrt(sum(delRij.^2, 2));

% Determine the maximum pairwise distance
if isempty(RMax), RMax = max(R); end

% Determine the number of bins
if isempty(numRBins)
    if normDists
        numRBins = ceil(4 * RMax);
    else
        numRBins = ceil(RMax / 2);
    end
end

%--------------------------------------------------------------------------
% Find the Axis of Symmetry of the Distribution Functions
%--------------------------------------------------------------------------

% Construct a set of sampling points over the polar angle for a fixed
% radius
theta = linspace(-pi, pi, numThetaBins+1).'; theta(end) = [];

% Construct a grid of bin edges
XBinEdges = linspace(-RMax, RMax, numGridBins);
YBinEdges = linspace(-RMax, RMax, numGridBins);

% Create a set of points on a ring several time larger than the average
% cell size
symR = symR * ones(numThetaBins, 1);
if ~normDists, symR = cellSize * symR; end
[x, y] = pol2cart(theta, symR);

% Determine the bin containing each sample point
[~, XInBin] = max(XBinEdges > x, [], 2); XInBin = XInBin-1;
[~, YInBin] = max(YBinEdges > y, [], 2); YInBin = YInBin-1;

% Convert these bin labels to linear indices
sampleBin = sub2ind( [length(YBinEdges)-1, length(XBinEdges)-1], ...
    XInBin, YInBin );

% Convert separation vectors into histogram
histRij = histcounts2( delRij(:,1), delRij(:,2), ...
    XBinEdges, YBinEdges );

% Extract angular histogram
histTheta = histRij(sampleBin);

% Determine the angle associated to the maximum value
[~, maxID] = max(histTheta);
thetaMax = theta(maxID);
XAngle = 0;
YAngle = pi/2;

%--------------------------------------------------------------------------
% Construct the Cut of the Radial Distribution Function
%--------------------------------------------------------------------------

% The width of the cut
if isempty(cutWidth)
    cutWidth = 0.25;
    if ~normDists, cutWidth = cutWidth * cellSize; end
end

% Construct radial sampling bins
r = linspace(0, RMax, numRBins);
dr = r(2)-r(1);

% Determine which pairs lie within the correct angular range
hasTheta_Max = ( (tan(thetaMax) .* delRij(:,1) - cutWidth / abs(cos(thetaMax))) <= delRij(:,2) ) ...
    & ( delRij(:,2) <= (tan(thetaMax) .* delRij(:,1) + cutWidth / abs(cos(thetaMax))) ) ...
    & ( (cos(thetaMax) .* delRij(:,1) + sin(thetaMax) .* delRij(:,2)) > 0 );

hasTheta_X = ( (tan(XAngle) .* delRij(:,1) - cutWidth / abs(cos(XAngle))) <= delRij(:,2) ) ...
    & ( delRij(:,2) <= (tan(XAngle) .* delRij(:,1) + cutWidth / abs(cos(XAngle))) ) ...
    & ( (cos(XAngle) .* delRij(:,1) + sin(XAngle) .* delRij(:,2)) > 0 );

hasTheta_Y = ( (tan(YAngle) .* delRij(:,1) - cutWidth / abs(cos(YAngle))) <= delRij(:,2) ) ...
    & ( delRij(:,2) <= (tan(YAngle) .* delRij(:,1) + cutWidth / abs(cos(YAngle))) ) ...
    & ( (cos(YAngle) .* delRij(:,1) + sin(YAngle) .* delRij(:,2)) > 0 );

% Accumulate Cell Counts --------------------------------------------------

% The correlation functions
GR_Max = zeros(1, numRBins);
GR_X = zeros(1, numRBins);
GR_Y = zeros(1, numRBins);
GR_All = zeros(1, numRBins);

% Normalize cell center coordinates
if normDists, X = X / cellSize; end

if useEffectiveVolume

    if strcmpi(effectiveVolumeMethod, 'grid')

        allVij = computeShellVolumeIntersectionGrid( ...
            bdyPoly, X(incNodeIDx,:), r, gridSize );

        allVij_Max = computeLinearShellVolumeIntersectionGrid( ...
            bdyPoly, X(incNodeIDx,:), r, thetaMax, cutWidth, gridSize );
        allVij_X = computeLinearShellVolumeIntersectionGrid( ...
            bdyPoly, X(incNodeIDx,:), r, XAngle, cutWidth, gridSize );
        allVij_Y = computeLinearShellVolumeIntersectionGrid( ...
            bdyPoly, X(incNodeIDx,:), r, YAngle, cutWidth, gridSize );

    elseif strcmpi(effectiveVolumeMethod, 'polygon')

        allVij = computeShellVolumeIntersection( ...
            bdyPoly, X(incNodeIDx,:), r );

        allVij_Max = computeLinearShellVolumeIntersection( ...
            bdyPoly, X(incNodeIDx,:), r, thetaMax, cutWidth );
        allVij_X = computeLinearShellVolumeIntersection( ...
            bdyPoly, X(incNodeIDx,:), r, XAngle, cutWidth );
        allVij_Y = computeLinearShellVolumeIntersection( ...
            bdyPoly, X(incNodeIDx,:), r, YAngle, cutWidth );

    else

        error('EV: This should not be reachable');

    end

    % Determine which bin each bond length belongs in
    [~, ~, rID] = histcounts(R, [r, (r(end)+dr)]);

    bondVij = sub2ind(size(allVij), J, rID);
    bondVij = 1./allVij(bondVij);
    bondVij(isinf(bondVij) | isnan(bondVij)) = 0;
    GR_All = full(sparse(1, rID, bondVij, 1, numRBins));

    rID_Max = hasTheta_Max .* rID;
    J_Max = hasTheta_Max .* J;
    J_Max(rID_Max == 0) = [];
    rID_Max(rID_Max == 0) = [];
    bondVij_Max = sub2ind(size(allVij_Max), J_Max, rID_Max);
    bondVij_Max = 1./ allVij_Max(bondVij_Max);
    bondVij_Max(isinf(bondVij_Max) | isnan(bondVij_Max)) = 0;
    GR_Max = full(sparse(1, rID_Max, bondVij_Max, 1, numRBins));

    rID_X = hasTheta_X .* rID;
    J_X = hasTheta_X .* J;
    J_X(rID_X == 0) = [];
    rID_X(rID_X == 0) = [];
    bondVij_X = sub2ind(size(allVij_X), J_X, rID_X);
    bondVij_X = 1./ allVij_X(bondVij_X);
    bondVij_X(isinf(bondVij_X) | isnan(bondVij_X)) = 0;
    GR_X = full(sparse(1, rID_X, bondVij_X, 1, numRBins));

    rID_Y = hasTheta_Y .* rID;
    J_Y = hasTheta_Y .* J;
    J_Y(rID_Y == 0) = [];
    rID_Y(rID_Y == 0) = [];
    bondVij_Y = sub2ind(size(allVij_Y), J_Y, rID_Y);
    bondVij_Y = 1./ allVij_Y(bondVij_Y);
    bondVij_Y(isinf(bondVij_Y) | isnan(bondVij_Y)) = 0;
    GR_Y = full(sparse(1, rID_Y, bondVij_Y, 1, numRBins));

    % Normalize results 
    GR_All = GR_All / (numDensity * (N-1));
    GR_Max = GR_Max / (numDensity * (N-1));
    GR_X = GR_X / (numDensity * (N-1));
    GR_Y = GR_Y / (numDensity * (N-1));

else

    allHasR = false(numel(R), numRBins);
    for i = 1:numRBins

        hasR = (r(i) <= R) & (R < (r(i)+dr));
        allHasR(:,i) = hasR;

        GR_Max(i) = sum(hasR .* hasTheta_Max);
        GR_X(i) = sum(hasR .* hasTheta_X);
        GR_Y(i) = sum(hasR .* hasTheta_Y);
        GR_All(i) = sum(hasR);

    end

end

%--------------------------------------------------------------------------
% Normalize Correlations
%--------------------------------------------------------------------------

cellsFromPair = @(x) (1+sqrt(1+4*x))/2;

switch normType

    case 'none'

        % This is a simple histogram count in each bin
        maxNorm = 1;
        XNorm = 1;
        YNorm = 1;
        allNorm = 1;

    case 'probability'

        % Bins are normalized by the total number of pairs of cells
        % included in the distribution function
        maxNorm = sum(GR_Max);
        XNorm = sum(GR_X);
        YNorm = sum(GR_Y);
        allNorm = sum(GR_All);

    case 'pairdist'

        % Bins are normalized so that N * g(r) * dr gives the number of
        % pairs of cells in the bin from [r, r+dr)
        maxNorm = cellsFromPair(sum(GR_Max)) * dr;
        XNorm = cellsFromPair(sum(GR_X)) * dr;
        YNorm = cellsFromPair(sum(GR_Y)) * dr;
        allNorm = cellsFromPair(sum(GR_All)) * dr;

    case 'radialdist'

        % Bins are normalized so that N * g(r) * rho * dV gives the number
        % of pairs of cells in the bin from [r, r+dr). Here rho is the cell
        % number density and dV is a volume element associated to the bin

        % Why is there a factor of 1/2 for these?
        maxNorm = N * numDensity * (dr/2);
        XNorm = N * numDensity * (dr/2);
        YNorm = N * numDensity * (dr/2);

        allNorm = (N * numDensity * (2 * pi * r * dr));

    otherwise

        error('Norm: This should not be reachable');

end

GR_Max = GR_Max ./ maxNorm;
GR_X = GR_X ./ XNorm;
GR_Y = GR_Y ./ YNorm;
GR_All = GR_All ./ allNorm;

%--------------------------------------------------------------------------
% Correct for Finite Size Effects
%--------------------------------------------------------------------------
% NOTE: Assumes the points in X have already been normalized by the cell
% size if necessary

if correctFiniteSize

    % Cut Along X-Direction -----------------------------------------------
    ksi = range(X(incNodeIDx,1));

    r1IDx = r <= ksi; r1 = r(r1IDx);
    r2IDx = ksi <= r; % r2 = r(r2IDx);

    finiteSizeFactor = zeros(size(GR_X));
    finiteSizeFactor(r1IDx) = 2 * (1-abs(r1/ksi));
    finiteSizeFactor(r2IDx) = 0;

    finiteSizeFactor = 1 ./ finiteSizeFactor;

    GR_X = GR_X .* finiteSizeFactor;

    % Cut Along Y-Direction -----------------------------------------------
    ksi = range(X(incNodeIDx,2));

    r1IDx = r <= ksi; r1 = r(r1IDx);
    r2IDx = ksi <= r; % r2 = r(r2IDx);

    finiteSizeFactor = zeros(size(GR_Y));
    finiteSizeFactor(r1IDx) = 2 * (1-abs(r1/ksi));
    finiteSizeFactor(r2IDx) = 0;

    finiteSizeFactor = 1 ./ finiteSizeFactor;

    GR_Y = GR_Y .* finiteSizeFactor;

    % Full 2D Radial Distribtion ------------------------------------------

    ksi = min(range(X(incNodeIDx,1)), range(X(incNodeIDx,2)));
    eta = max(range(X(incNodeIDx,1)), range(X(incNodeIDx,2)));

    r1IDx = (r <= ksi) & (r <= eta);
    r2IDx = (ksi <= r) & (r <= eta);
    r3IDx = (ksi <= r) & (eta <= r) & (r <= sqrt(ksi^2 + eta^2));
    r4IDx = (sqrt(ksi^2 + eta^2) <= r);

    r1 = r(r1IDx); r2 = r(r2IDx); r3 = r(r3IDx); % r4 = r(r4IDx);

    finiteSizeFactor = zeros(size(GR_All));
    finiteSizeFactor(r1IDx) = ...
        4 * r1 .* (pi/2 - r1/ksi - r1/eta + r1.^2/(2*ksi*eta));
    finiteSizeFactor(r2IDx) = ...
        4 * r2 .* (pi/2 - r2/ksi - acos(ksi./r2) - ...
        ksi/(2*eta) + sqrt(r2.^2/ksi^2-1));
    finiteSizeFactor(r3IDx) = ...
        4 * r3 .* (pi/2 - acos(eta./r3) - acos(ksi./r3) - eta/(2*ksi) - ...
        ksi/(2*eta) + sqrt(r3.^2/ksi^2-1) + ...
        sqrt(r3.^2/eta^2-1)-r3.^2/(2*ksi*eta));
    finiteSizeFactor(r4IDx) = 0;

    finiteSizeFactor = (2*pi*r) ./ finiteSizeFactor;
    finiteSizeFactor(isinf(finiteSizeFactor) | isnan(finiteSizeFactor)) = 0;

    GR_All = finiteSizeFactor .* GR_All;

end


%--------------------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------------------
if plotResults
    
    figure
    
    subplot(2,3,1)
    
    plot(r, GR_Max, '-k', 'LineWidth', 2);
    
    xlim([0.5 RMax]);
    
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    
    grid on
    
    if normDists
        ylabel('g({\itr/L}, \theta_{Max})');
        xlabel('{\itr/L} (Average cell lengths)');
    else
        ylabel('g({\itr}, \theta_{Max})');
        xlabel('{\itr} (Raw distance)');
    end
    
    subplot(2,3,2)
    
    plot(r, GR_X, '-k', 'LineWidth', 2);
    
    xlim([0.5 RMax]);
    
    % set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    
    grid on
    
    if normDists
        ylabel('g({\itr/L}, \theta_{X})');
        xlabel('{\itr/L} (Average cell lengths)');
    else
        ylabel('g({\itr}, \theta_{X})');
        xlabel('{\itr} (Raw distance)');
    end
    
    subplot(2,3,3)
    
    plot(r, GR_Y, '-k', 'LineWidth', 2);
    
    xlim([0.5 RMax]);
    
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    
    grid on
    
    if normDists
        ylabel('g({\itr/L}, \theta_{Y})');
        xlabel('{\itr/L} (Average cell lengths)');
    else
        ylabel('g({\itr}, \theta_{Y})');
        xlabel('{\itr} (Raw distance)');
    end
    
    subplot(2,3,4)
    
    scatter( delRij(hasTheta_Max, 1), delRij(hasTheta_Max, 2), ...
        'filled', 'r' );
    
    hold on
    scatter( delRij(~hasTheta_Max, 1), delRij(~hasTheta_Max, 2), ...
        'filled', 'k' );
    hold off
    
    axis equal
    xlim([-RMax, RMax]);
    ylim([-RMax, RMax]);
    
    title('Cut Along Axis of Symmetry');
    
    subplot(2,3,5)
    
    scatter( delRij(hasTheta_X, 1), delRij(hasTheta_X, 2), ...
        'filled', 'r' );
    
    hold on
    scatter( delRij(~hasTheta_X, 1), delRij(~hasTheta_X, 2), ...
        'filled', 'k' );
    hold off
    
    axis equal
    xlim([-RMax, RMax]);
    ylim([-RMax, RMax]);
    
    title('Cut Along X-Axis');
    
    subplot(2,3,6)
    
    scatter( delRij(hasTheta_Y, 1), delRij(hasTheta_Y, 2), ...
        'filled', 'r' );
    
    hold on
    scatter( delRij(~hasTheta_Y, 1), delRij(~hasTheta_Y, 2), ...
        'filled', 'k' );
    hold off
    
    axis equal
    xlim([-RMax, RMax]);
    ylim([-RMax, RMax]);
    
    title('Cut Along Y-Axis');

end

end

function bdyPoly = ...
    extractBoundaryPolygon(delTri, v, c, incNodeIDx, plotResults)
%EXTRACTBOUNDARYPOLYGON Extract the polyshape definining the 
%boundary of a subset of cells from a Voronoi tesselation
%
%   INPUT PARAMETERS:
%
%       - delTri:       The Delaunay triangulation dual to the Voronoi
%                       tesselation
%       - v:            #Vx2 Voronoi vertex coordinates
%       - c:            #Cx1 cell array of vertex IDs defining the
%                       ordered polygon of each Voronoi cell
%       - incNodeIDx:   #ICx1 list of included cells
%       - plotResults:  If true, the final results will be plotted
%
%   OUTPUT PARAMETERS:
%
%       - bdyPoly:      MATLAB-style polyshape defining the boundary
%
%   by Dillon Cislo 05/18/2022

if (nargin < 5), plotResults = false; end

% The cell IDs for each triangulation edge
bondIDx = delTri.edges;

% Determine which cells lie on the boundary of the diagram
bdyCells = cellfun( @(x) ismember(1, x), c );
includedBdyIDx = ismember(incNodeIDx, find(bdyCells));
if any(includedBdyIDx)
    warning(['True boundary cells cannot be included in ' ...
        'bondary polygon computation']);
    incNodeIDx(includedBdyIDx) = [];
end

% The vertices of the Voronoi domain intersected by the corresponding edge
% connecting cell centers
bondVIDx = zeros(size(bondIDx));
for i = 1:size(bondIDx, 1)

    cellIntersect = intersect( c{bondIDx(i,1)}, c{bondIDx(i,2)} );
    if ( numel(cellIntersect) > 2 )
        cellIntersect( cellIntersect == 1 ) = [];
    end

    bondVIDx(i,:) = cellIntersect;

end

% Boundary edges are defined to be those where one member of the edge is
% included in the subregion and the other is not
bdyEdgeIDx = sum(ismember(bondIDx, incNodeIDx), 2) == 1;
bdyEdges = bondVIDx(bdyEdgeIDx, :);

% The vertices definiing the boundary edges
bdyVertIDx = unique(bdyEdges(:));
bdyVerts = v(bdyVertIDx, :);

% Local IDs for the boundary vertices
locVertIDx = (1:numel(bdyVertIDx)).';

% Update bonds to be defined in terms of local IDs
bdyEdges = changem(bdyEdges, locVertIDx, bdyVertIDx);

% Construct a graph representation of the boundary polygon
A = sparse( bdyEdges(:), [bdyEdges(:,2); bdyEdges(:,1)], ...
    1, numel(locVertIDx), numel(locVertIDx) );
if any(sum(A,2) ~= 2)
    error('Vertices belong to more than two edges')
end
G = graph(A);

% Sort the edge list
newBdyEdges = dfsearch(G, 1, {'edgetonew'});
newBdyEdges = [ newBdyEdges; newBdyEdges(end,2), newBdyEdges(1) ];

bdyPoly = polyshape(bdyVerts(newBdyEdges(:,1), 1), ...
    bdyVerts(newBdyEdges(:,1), 2));

if plotResults

    voronoi(delTri, 'k');
    hold on
    scatter(delTri.Points(incNodeIDx,1), delTri.Points(incNodeIDx,2), ...
        'filled', 'r');
    plot(bdyPoly);
    hold off
    axis equal

end


end

function allVij = computeShellVolumeIntersection(P, X, r)
%COMPUTESHELLVOLUMEINTERSECTION Compute the area of intersection of an
%annular ring [r, r+dr) centered on a set of cells for a given set of
%radii (all differing by dr) with a polygon
%
%   INPUT PARAMETERS:
%
%       - P:    MATLAB-style polyshape representation of polygonal region
%       - X:    #Xx2 set of cell center coordinates
%       - r:    1x#R list of radii. r = [0 dr 2*dr ... R]
%
%   OUTPUT PARAMETERS:
%
%       - allVij:   #Xx#R set of intersection volumes
%
%   by Dillon Cislo 05/18/2022

% The difference between subsequent radii (assumed to be uniform)
dr = r(2)-r(1);

% Compute the area of intersection for each circle
try

    % Re-format the base polygon to match the Boost style
    % NOTE: Boost-style polygons are CLOSED and oriented CW
    PC = P.Vertices;
    PC = [PC; PC(1,:)];
    [PC(:,1), PC(:,2)] = poly2cw(PC(:,1), PC(:,2));
    
    % Generate a polygonal approximation of the unit circle
    % NOTE: Boost-style polygons are CLOSED and oriented CW
    numCircPnts = 101;
    circPnts = linspace(0, 2*pi, numCircPnts).';
    circPnts = flipud([cos(circPnts), sin(circPnts)]);

    QX = repmat(circPnts(:,1), 1, size(X,1) * numel(r));
    QX = repmat(repelem((r+dr), size(X,1)), numCircPnts, 1) .* QX;
    QX = QX + repmat(X(:,1).', numCircPnts, numel(r));

    QY = repmat(circPnts(:,2), 1, size(X,1) * numel(r));
    QY = repmat(repelem((r+dr), size(X,1)), numCircPnts, 1) .* QY;
    QY = QY + repmat(X(:,2).', numCircPnts, numel(r));

    allVij = polygon_intersection_area(PC, QX, QY);
    allVij = reshape(allVij, size(X,1), numel(r));
    allVij = [zeros(size(X,1), 1), allVij];

catch

    warning('Boost method failed - using MATLAB built-ins');

    % Generate a polygonal approximation of the unit circle
    numCircPnts = 100;
    circPnts = linspace(0, 2*pi, numCircPnts+1).'; circPnts(end) = [];
    circPnts = [cos(circPnts), sin(circPnts)];

    allVij = zeros(size(X,1), numel(r)+1);
    for i = 1:size(X,1)
        % progressbar(i, size(X,1)); % Needs 'gptoolbox'
        for j = 2:(numel(r)+1)

            curCircPnts = (r(j-1)+dr) * circPnts + X(i,:);
            circPoly = polyshape(curCircPnts(:,1), curCircPnts(:,2));

            allVij(i,j) = area(intersect(circPoly, P));

        end
    end

end

% Subtract circle intersection areas to find annular ring intersection
% areas
allVij = diff(allVij, 1, 2);

end

function allVij = ...
    computeLinearShellVolumeIntersection(P, X, r, ang, h)
%COMPUTELINEARSHELLVOLUMEINTERSECTION Compute the area of intersection of a
%rectangular region a distance [r, r+dr) centered on a set of cells for a
%given set of radii (all differing by dr) with a polygon
%
%   INPUT PARAMETERS:
%
%       - P:    MATLAB-style polyshape representation of polygonal region
%       - X:    #Xx2 set of cell center coordinates
%       - r:    1x#R list of radii. r = [0 dr 2*dr ... R]
%       - ang:  The orientation angle of the rectangular region
%       - h:    The thickness of the rectangular region
%
%   OUTPUT PARAMETERS:
%
%       - allVij:   #Xx#R set of intersection volumes
%
%   by Dillon Cislo 05/19/2022

% The difference between subsequent radii (assumed to be uniform)
dr = r(2)-r(1);

% Generate a polygonal represenation of the rectangular region with
% width dr and height h centered at (0,0) with horizontal axis along
% the orientation specified by ang
rectPnts = [-1 -1; 1 -1; 1 1; -1 1] / 2;
rectPnts(:,1) = dr * rectPnts(:,1);
rectPnts(:,2) = h * rectPnts(:,2);
rectPnts = exp(1i * ang) * complex(rectPnts(:,1), rectPnts(:,2));
rectPnts = [real(rectPnts), imag(rectPnts)];

% The unit vector along the orientation direction
n = [cos(ang), sin(ang)];

% Compute the area of intersection for each cell and radius
try

    % Re-format the base polygon to match the Boost style
    % NOTE: Boost-style polygons are CLOSED and oriented CW
    PC = P.Vertices;
    PC = [PC; PC(1,:)];
    [PC(:,1), PC(:,2)] = poly2cw(PC(:,1), PC(:,2));

    % Re-format the rotated rectangles to match the Boost style
    rPnts = flipud([rectPnts; rectPnts(1,:)]);

    QX1 = repmat(rPnts(:,1), 1, size(X,1) * numel(r));
    QX1 = QX1 + repmat(X(:,1).', 5, numel(r));
    QX2 = QX1;

    rShiftX = repmat(repelem((r+dr/2) * n(1), size(X,1)), 5, 1);
    QX1 = QX1 + rShiftX;
    QX2 = QX2 - rShiftX;

    QY1 = repmat(rPnts(:,2), 1, size(X,1) * numel(r));
    QY1 = QY1 + repmat(X(:,2).', 5, numel(r));
    QY2 = QY1;

    rShiftY = repmat(repelem((r+dr/2) * n(2), size(X,1)), 5, 1);
    QY1 = QY1 + rShiftY;
    QY2 = QY2 - rShiftY;

    QX = [QX1, QX2];
    QY = [QY1, QY2];

    allVij = polygon_intersection_area(PC, QX, QY);
    allVij = reshape(allVij, size(X,1), numel(r), 2);
    allVij = sum(allVij, 3);

catch

    warning('Boost method failed - using MATLAB built-ins');

    allVij = zeros(size(X,1), numel(r));
    for i = 1:size(X,1)
        % progressbar(i, size(X,1)); % Needs 'gptoolbox'
        for j = 1:numel(r)

            curRectPnts1 = rectPnts + (X(i,:) + (r(j)+dr/2) * n);
            rectPoly1 = polyshape(curRectPnts1(:,1), curRectPnts1(:,2));

            curRectPnts2 = rectPnts + (X(i,:) - (r(j)+dr/2) * n);
            rectPoly2 = polyshape(curRectPnts2(:,1), curRectPnts2(:,2));

            allVij(i,j) = area(intersect(rectPoly1, P)) + ...
                area(intersect(rectPoly2, P));

        end
    end

end

end

function allVij = computeShellVolumeIntersectionGrid(P, X, r, gridSize)
%COMPUTESHELLVOLUMEINTERSECTION Compute the area of intersection of an
%annular ring [r, r+dr) centered on a set of cells for a given set of
%radii (all differing by dr) with a polygon
%
%   INPUT PARAMETERS:
%
%       - P:    MATLAB-style polyshape representation of polygonal region
%       - X:        #Xx2 set of cell center coordinates
%       - r:        1x#R list of radii. r = [0 dr 2*dr ... R]
%       - gridSize: 1x2 dimensions of the sampling grid
%
%   OUTPUT PARAMETERS:
%
%       - allVij:   #Xx#R set of intersection volumes
%
%   by Dillon Cislo 05/18/2022

if (nargin < 4), gridSize = [100 100]; end

% The difference between subsequent radii (assumed to be uniform)
dr = r(2)-r(1);

% Generate a grid in the tight bounding box of the polygon
minX = min(P.Vertices(:,1)); maxX = max(P.Vertices(:,1));
minY = min(P.Vertices(:,2)); maxY = max(P.Vertices(:,2));

numX = gridSize(1);
numY = gridSize(2);         
[XGrid, YGrid] = ...
    meshgrid(linspace(minX, maxX, numX), linspace(minY, maxY, numY));

% Determine the area of each grid pixel
dx = XGrid(1,2)-XGrid(1,1);
dy = YGrid(2,1)-YGrid(1,1);
areaPix = dx * dy;

% Determine which grid pixels lie within the polygon
inPoly = reshape(P.isinterior(XGrid(:), YGrid(:)), size(XGrid));

% Compute the area of intersection for each annular ring
allVij = zeros(size(X,1), numel(r));
for i = 1:size(X,1)
    % progressbar(i, size(X,1)); % Needs 'gptoolbox'
    shiftXY2 = (XGrid-X(i,1)).^2 + (YGrid-X(i,2)).^2;
    for j = 1:numel(r)
        inRing = (r(j)^2 <= shiftXY2) & (shiftXY2 < (r(j)+dr)^2);
        allVij(i,j) = areaPix * sum(inRing(:) & inPoly(:));
    end
end

end

function allVij = ...
    computeLinearShellVolumeIntersectionGrid(P, X, r, ang, h, gridSize)
%COMPUTELINEARSHELLVOLUMEINTERSECTION Compute the area of intersection of a
%rectangular region a distance [r, r+dr) centered on a set of cells for a
%given set of radii (all differing by dr) with a polygon
%
%   INPUT PARAMETERS:
%
%       - P:    MATLAB-style polyshape representation of polygonal region
%       - X:        #Xx2 set of cell center coordinates
%       - r:        1x#R list of radii. r = [0 dr 2*dr ... R]
%       - ang:      The orientation angle of the rectangular region
%       - h:        The thickness of the rectangular region
%       - gridSize: 1x2 dimensions of the sampling grid
%
%   OUTPUT PARAMETERS:
%
%       - allVij:   #Xx#R set of intersection volumes
%
%   by Dillon Cislo 05/19/2022

if (nargin < 6), gridSize = [100 100]; end

% The difference between subsequent radii (assumed to be uniform)
dr = r(2)-r(1);

% Generate a grid in the tight bounding box of the polygon
minX = min(P.Vertices(:,1)); maxX = max(P.Vertices(:,1));
minY = min(P.Vertices(:,2)); maxY = max(P.Vertices(:,2));

numX = gridSize(1);
numY = gridSize(2);         
[XGrid, YGrid] = ...
    meshgrid(linspace(minX, maxX, numX), linspace(minY, maxY, numY));

% Determine the area of each grid pixel
dx = XGrid(1,2)-XGrid(1,1);
dy = YGrid(2,1)-YGrid(1,1);
areaPix = dx * dy;

% Determine which grid pixels lie within the polygon
inPoly = reshape(P.isinterior(XGrid(:), YGrid(:)), size(XGrid));

% Compute the area of intersection for each cell and radius
allVij = zeros(size(X,1), numel(r));
for i = 1:size(X,1)
    % progressbar(i, size(X,1)); % Needs 'gptoolbox'

    % Shift and rotate the grid points
    shiftX = XGrid-X(i,1);
    shiftY = YGrid-X(i,2);
    shiftXY = exp(-1i * ang) * complex(shiftX, shiftY);
    absShiftX = abs(real(shiftXY));
    absShiftY = abs(imag(shiftXY));

    for j = 1:numel(r)

        inBox = (r(j) <= absShiftX) & (absShiftX < (r(j)+dr));
        inBox = inBox & (absShiftY <= (h/2));

        allVij(i,j) = areaPix * sum(inBox(:) & inPoly(:));

    end
end
    
end


