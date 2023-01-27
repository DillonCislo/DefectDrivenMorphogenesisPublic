function [psiCorr, pairCorrs, numCellsR, r] = ...
    calculateBondOrientationalCorrelations(X, psi, varargin)
%CALCULATEBONDORIENTATIONALCORRELATIONS Calculate the two-point correlation
%function of the bond orientational order parameter of a set of points. It
%is assumed that the correlations depend only on the separation between two
%points
%
%   INPUT PARAMETERS:
%
%       - X:            #Px2 list of scattered point set (x,y) coordinates
%
%       - psi:          #Px1 list of complex bond orientational order
%                       parameters. Entries corresponding to boundary
%                       points of the Delaunay triangulation may be NaN
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
%       - ('IncludedNodes', incNodeIDx = []): #IPx1 vector of point IDs
%       indicated which nodes to include in the calculation
%
%       - ('OrderType', n = 0): The type of orientational order parameter
%       used to calculate the correlation function
%
%       - ('NormalizeDistances', normDists = true): If true, distances
%       between points are normalized by the cellSize parameter
%
%       - ('PlotResults', plotResults = false): If true, radial
%       correlations are plotted
%
%   OUTPUT PARAMETERS:
%
%       - psiCorr:      1xnumSamples two point correlation functions
%                       evaluated at a discrete set of distances
%
%       - pairCorrs:    (#P*(#P-1)/2)x1 list of all pairwise
%                       correlations among the input point set
%
%       - numCellsR:    1xnumSamples count of pairwise distances included
%                       in each bin used to calculate 'psiCorr' 
%
%   by Dillon Cislo 03/01/2021

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if (nargin < 1), error('Please input point set coordinates'); end
if (nargin < 2), error('Please input point set order parameters'); end

validateattributes(X, {'numeric'}, {'2d', 'ncols', 2, 'finite', 'real'});
validateattributes(psi, {'numeric'}, {'vector', 'numel', size(X,1)});
if (size(psi, 2) ~= 1), psi = psi.'; end

% Set default values
cellSize = [];
RMax = [];
numRBins = [];
incNodeIDx = [];
orderType = 0;
normDists = true;
plotResults = false;

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
    
    if strcmpi(varargin{i}, 'IncludedNodes')
        incNodeIDx = varargin{i+1};
        if ~isempty(incNodeIDx)
            validateattributes(incNodeIDx, {'numeric'}, ...
                {'vector', 'integer', 'positive', 'finite', ...
                'real', '<=', size(X,1)});
        end
    end
    
    if strcmpi(varargin{i}, 'OrderType')
        orderType = varargin{i+1};
        validateattributes(orderType, {'numeric'}, {'scalar'});
    end
    
    if strcmpi(varargin{i}, 'NormalizeDistances')
        normDists = varargin{i+1};
        validateattributes(normDists, {'logical'}, {'scalar'});
    end
    
    if strcmpi(varargin{i}, 'PlotResults')
        plotResults = varargin{i+1};
        validateattributes(plotResults, {'logical'}, {'scalar'});
    end
    
end

%--------------------------------------------------------------------------
% Calculate Pairwise Order Correlations for All Bonds
%--------------------------------------------------------------------------

% Construct a Delaunay triangulation of the input point set
delTri = delaunayTriangulation(X);
[~, c] = voronoiDiagram(delTri);

% The edge list of the triangulation
bondIDx = delTri.edges;

% Determine which cells lie on the boundary of the diagram
bdyCells = cellfun( @(x) ismember(1, x), c );

% Construct a list of all pairs of cell IDs
pairIDx = nchoosek( (1:size(X,1)).', 2 );

% Calculate the pairwise correlations
pairCorrs = psi(pairIDx(:,1)) .* conj(psi(pairIDx(:,2)));

% Calculate pairwise distances between cells
pairDist = X(pairIDx(:,2), :) - X(pairIDx(:,1), :);
pairDist = sqrt(sum(pairDist.^2, 2));

%--------------------------------------------------------------------------
% Construct Full Correlation Functions
%--------------------------------------------------------------------------

% Determine which subset of nodes to include in the calculation
if isempty(incNodeIDx), incNodeIDx = find(~bdyCells); end
assert(~any(isnan(psi(incNodeIDx))), ...
    'Included points contain NaN values as order parameters');

% Determine which pairs consist only of included nodes
incPairIDx = all( ismember(pairIDx, incNodeIDx), 2 );

% Calculate an 'average cell size' based on pairwise distances
if isempty(cellSize)

    % Determine which pairs correspond to edges of the Delaunay
    % triangulation
    isEdge = ismember( sort(pairIDx, 2), sort(bondIDx,2), 'rows' );
    
    cellSize = mean(pairDist(isEdge & incPairIDx));

end

% cellSize

% Normalize distances by average cell size
if normDists, pairDist = pairDist ./ cellSize; end

% Determine the maximum pairwise distance
if isempty(RMax), RMax = max(pairDist(incPairIDx)); end

% RMax

% Determine the number of bins
if isempty(numRBins)
    if normDists
        numRBins = ceil(4 * RMax);
    else
        numRBins = ceil(RMax / 2);
    end
end

% The radial sample points
r = linspace(0, RMax, numRBins);
dr = r(2) - r(1);

psiCorr = zeros(1, numRBins);
numCellsR = zeros(1, numRBins);

for i = 1:numRBins
    
    % The IDs of pairwise distances in the current bin
    hasR = (r(i) <= pairDist) & (pairDist < (r(i)+dr));
    
    numCellsR(i) = sum(hasR & incPairIDx);
    psiCorr(i) = mean( real( pairCorrs(hasR & incPairIDx) ) );
    
end

psiCorr( isnan(psiCorr) | isinf(psiCorr) ) = 0;

%--------------------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------------------
if plotResults
    
    figure;
    
    plotCorr = psiCorr;
    badR = numCellsR < -1;
    plotCorr(badR) = 0;
    
    plot(r, plotCorr, 'LineWidth', 1.5);
    
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    
    xlim([0 RMax]);
    % ylim([min(plotCorr) max(plotCorr)]);
    
    grid on
    
    if normDists
        ylabel(sprintf('\\langle\\psi_{%d}({\\itr/L})\\psi_{%d}(0)\\rangle', ...
            orderType, orderType));
        xlabel('{\it r/L} (Average cell lengths)');
    else
        ylabel(sprintf('\\langle\\psi_{%d}({\\itr})\\psi_{%d}(0)\\rangle', ...
            orderType, orderType));
        xlabel('{\it r} (Raw distance)');
    end
    
    title('Bond Orientational Correlation Function');
    
end

end

