function psiN = ...
    calculateBondOrientationalOrder(X, n, m, smoothOutput, plotResults)
%CALCULATEBONDORIENTATIONALORDER Calculate the n-fold bond orientational
%order parameter for a scattered point set. Point set connectivity is
%inferred by Voronoi tessellation
%
%   INPUT PARAMETERS:
%
%       - X:            #Px2 list of scattered point set (x,y) coordinates
%
%       - n:            Specifies which n-fold order parameter is
%                       calculated
%
%       - m:            The edge weighting power
%
%       - smoothOutput: If true, output is smoothed by Laplacian mesh
%                       smoothing on the Delaunay triangulation of the
%                       input point set
%
%       - plotResults:  If true, calculated order parameters are plotted
%
%   OUTPUT PARAMETERS
%
%       - psiN:         #Px1 complex order parameters for each point
%                       Points on the boundary of the Delaunay
%                       triangulation have their order parameters set to
%                       NaN
%
%   by Dillon Cislo 02/27/2021

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if (nargin < 1), error('Please input point set coordinates'); end
if (nargin < 2), error('Please input n-fold order type'); end
if (nargin < 3), m = 2; end
if (nargin < 4), smoothOutput = false; end
if (nargin < 5), plotResults = false; end

validateattributes(X, {'numeric'}, ...
    {'2d', 'ncols', 2, 'finite', 'real'});
validateattributes(n, {'numeric'}, ...
    {'scalar', 'positive', 'integer', 'finite', 'real'});
validateattributes(m, {'numeric'}, ...
    {'scalar', 'nonnegative', 'integer', 'finite', 'real'});
validateattributes(smoothOutput, {'logical'}, {'scalar'});
validateattributes(plotResults, {'logical'}, {'scalar'});

%--------------------------------------------------------------------------
% Process Node Connectivity
%--------------------------------------------------------------------------

% Construct Voronoi tessellation
delTri = delaunayTriangulation( X );
[v, c] = voronoiDiagram( delTri );

% The node IDs for each bond
bondIDx = delTri.edges;

% Determine which cells lie on the boundary of the diagram
bdyCells = cellfun( @(x) ismember(1, x), c );

% Convert the Voronoi cell connectivity list to a (NaN-padded) matrix
% Useful for plotting purposes
maxFaceSize = max(cellfun(@(x) numel(x), c));
voronoiFace = nan(size(c,1), maxFaceSize);
for i = 1:size(c,1)
    voronoiFace(i, 1:numel(c{i})) = c{i};
end

% Determine the connectivity of the Voronoi cells -------------------------

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

% The length of edge Voronoi edge
edgeLength = v(bondVIDx(:,2), :) - v(bondVIDx(:,1), :);
edgeLength = sqrt( sum( edgeLength.^2, 2 ) );

% The actual bond vector
bond = X(bondIDx(:,2), :) - X(bondIDx(:,1), :);

% The bond angle relative to the first cell in the bond
bondAngle = atan2(bond(:,2), bond(:,1));
bondAngle = wrapTo2Pi([bondAngle, pi+bondAngle]);

%--------------------------------------------------------------------------
% Calculate the Order Parameter
%--------------------------------------------------------------------------

psiN = zeros(size(X,1), 2);

for i = 1:size(bondIDx, 1)
    
    c1_IDx = bondIDx(i,1); % The ID of the first cell in the bond
    c2_IDx = bondIDx(i,2); % The ID of the second cell in the bond
    
    psiN( c1_IDx, 1 ) = psiN( c1_IDx, 1 ) + ...
        edgeLength(i).^m .* exp( 1i .* n .* bondAngle(i,1) );
    psiN( c1_IDx, 2 ) = psiN( c1_IDx, 2 ) + edgeLength(i).^m;
    
    psiN( c2_IDx, 1 ) = psiN( c2_IDx, 1 ) + ...
        edgeLength(i).^m .* exp( 1i .* n .* bondAngle(i,2) );
    psiN( c2_IDx, 2 ) = psiN( c2_IDx, 2 ) + edgeLength(i).^m;
    
end

% Normalize the order parameter
psiN = psiN(:,1) ./ psiN(:,2);

if smoothOutput
    
    % Boundary cells can have NaN entries. Set these to 0
    psiN( isnan(psiN) | isinf(psiN) ) = 0;
    
    try
        psiN = laplacian_smooth( X, delTri.ConnectivityList, ...
            'uniform', [], 0.001, 'implicit', psiN, 30 );
    catch
        warning(['Laplacian smoothing function not found! ' ...
            'Processing unsmoothed order parameters... ']);
    end
    
end

% Set boundary cell order parametres to NaN
psiN(bdyCells) = NaN;

%--------------------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------------------
if plotResults
    
    dispCells = ~bdyCells;
    
    % Color/Transparency Handling for Magnitude ---------------------------
    
    cmap = parula(256);
    
    cvals = abs(psiN);
    cvals( cvals < 0 ) = 0;
    cvals( cvals > 1 ) = 1;
    
    cvalsN = round( ((cvals-0) ./ (1-0)) .* (size(cmap,1)-1) ) + 1;
    
    absColors = nan(size(dispCells,1), 3);
    absColors(dispCells, :) = cmap(cvalsN(dispCells), :);
    
    absAlpha = zeros(size(X,1), 1);
    absAlpha(dispCells) = 1;
    
    % Color Transparency Handling for Angle -------------------------------
    
    cmap = phasemap(256);
    
    cvals = wrapTo2Pi(angle(psiN));
    
    cvalsN = round( ((cvals-0)) ./ (2*pi-0) .* (size(cmap,1)-1) ) + 1;
    
    angColors = nan(size(dispCells,1), 3);
    angColors(dispCells, :) = cmap(cvalsN(dispCells), :);
    
    angAlpha = zeros(size(X,1), 1);
    angAlpha(dispCells) = 1;
    
    % Generate Plots ------------------------------------------------------
    figure;
    
    subplot(1,2,1)
    
    patch( 'Faces', voronoiFace(dispCells, :), 'Vertices', v, ...
        'FaceVertexCData', absColors(dispCells, :), 'FaceColor', 'flat', ...
        'FaceVertexAlphaData', absAlpha(dispCells), 'FaceAlpha', 'flat', ...
        'AlphaDataMapping', 'none', 'EdgeColor', 'k', 'LineWidth', 2 );
    
    hold on
    
    scatter(X(bdyCells, 1), X(bdyCells, 2), 'filled', 'k');
    
    hold off
    
    axis equal
    xlim([min(X(:,1)) max(X(:,1))]);
    ylim([min(X(:,2)) max(X(:,2))]);
    
    colorbar
    set(gca, 'Clim', [0 1]);
    
    title( sprintf('\\mid \\psi_{%d} \\mid', n) );
    
    subplot(1,2,2)
    
    patch( 'Faces', voronoiFace(dispCells, :), 'Vertices', v, ...
        'FaceVertexCData', angColors(dispCells, :), 'FaceColor', 'flat', ...
        'FaceVertexAlphaData', angAlpha(dispCells), 'FaceAlpha', 'flat', ...
        'AlphaDataMapping', 'none', 'EdgeColor', 'k', 'LineWidth', 2 );
    
    hold on
    
    scatter(X(bdyCells, 1), X(bdyCells, 2), 'filled', 'k');
    
    hold off
    
    axis equal
    xlim([min(X(:,1)) max(X(:,1))]);
    ylim([min(X(:,2)) max(X(:,2))]);
    
    colorbar
    set(gca, 'Clim', [0 2*pi]);
    set(gca, 'Colormap', cmap);
    
    title( sprintf('arg[\\psi_{%d}]', n) );

end

end

