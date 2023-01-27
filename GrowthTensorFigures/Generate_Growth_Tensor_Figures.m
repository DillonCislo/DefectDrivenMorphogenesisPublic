%% Generate Growth Tensor Figures =========================================
% This is a script to generate some pretty figures describing the
% decomposition of a growth tensor into its isotropic and deviatoric parts.
% It generates the panels of Fig. 4F.
%
% by Dillon Cislo 2021/01/18
%==========================================================================

% Add relevant file structure to MATLAB path
[DDMDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
DDMDir = [DDMDir '/..'];
addpath(genpath(DDMDir));
clear DDMDir

%% Set Up Pipeline ========================================================
clear; close all; clc;

% Generate a Point Set for Cell Centers -----------------------------------

% Bounds of the lattice region
xBnd = [-100, 100];
yBnd = [-100, 100];

% Create rectangular lattice sites
nrows = 20; dy = ( yBnd(2)-yBnd(1) ) / ( nrows - 1 );
ncols = 20; dx = ( xBnd(2)-xBnd(1) ) / ( ncols - 1 );

numPoints = nrows * ncols;

[x, y] = meshgrid( (xBnd(1):dx:(xBnd(2))), ...
    (yBnd(1):dy:(yBnd(2))) );

x = [ x(:), y(:) ];

clear Y dx dy

% Add noise to the lattice sites
rng(1, 'twister');
x = awgn(x, 5);

% Construct Connectivity of the Undeformed Lattice ------------------------

delTri = delaunayTriangulation(x);
% F = delTri.ConnectivityList;

[v, c] = voronoiDiagram(delTri);

% Extract boundary cells
bdyFace = cellfun( @(x) ismember(1, x), c );

% Extract boundary nodes
bdyIDx = unique(freeBoundary(delTri));
bdyNodes = ismember((1:size(x,1)).', bdyIDx);

% Convert the voronoi cell connectivity list to a (NaN-padded)
% matrix Useful for plotting purposes
maxFaceSize = max(cellfun(@(x) numel(x), c));
voronoiFace = nan(size(c,1), maxFaceSize);
for i = 1:size(c,1)
    voronoiFace(i, 1:numel(c{i})) = c{i};
end

clear maxFaceSize

% Set Division Parameters -------------------------------------------------

% Match a desired division location to a node location
x0 = [37 -5];
divID = knnsearch(x, x0);
x0 = x(divID, :);

% Find the radius of the inclusion
divArea = polyarea(v(c{divID}, 1), v(c{divID}, 2));
a = sqrt(divArea / pi);

M = 1; %0.25; % The isotropic growth parameter
q = 0.25; %0.0625; % The deviatoric growth parameter
nu = 1/3; % The Poisson ratio
n = [0 1]; % The division axis
thin3D = true; % True 2D body or thin 3D body

% Set Visualization Parameters --------------------------------------------

windowSize = 75;

% cellColor = [134 0 187] / 255;
% cellColor = [203 61 223]/255;
cellColor = [242 171 6] / 255;

% arrowColor = [242 171 6] / 255;
arrowColor = [203 61 223]/255;
% arrowColor = [1 0 1];

pntColor = [223 61 80]/255;

patchLineWidth = 0.75;
arrowLineWidth = 1; %0.5;
maxArrowHeadSize = 1;

% arrowScaling = 4; % Raw length rescale
% arrowScaling = 0.5; % Internal scaling in 'quiver'

arrowScaling_Total = 4;
arrowScaling_Iso = 6;
arrowScaling_Ani = 10;

% In cm
% maxAxDim = 2.5;
% maxFigDim = 1.2 * maxAxDim;
maxFigDim = 3.5;


%% Visualize the Undeformed Lattice =======================================
close all; clc;

% % Construct a set of circle test points
% ntheta = 100;
% theta = linspace(0, 2*pi, ntheta+1); theta(end) = [];
% 
% circX = zeros(ntheta, 2);
% [circX(:,1), circX(:,2)] = pol2cart(theta, a);
% circX = circX + x0;

faceColors = ones( numel(bdyFace), 3 );
faceColors(divID, :) = cellColor;

% pntColor = [223 131 61]/255;
% pntColor = [80 223 61]/255;
% pntColor = [223 61 80]/255;


fig = figure('Color', [1 1 1]);

patch( 'Faces', voronoiFace(~bdyFace, :), 'Vertices', v, ...
    'FaceVertexCData', faceColors(~bdyFace, :), 'FaceColor', 'flat', ...
    'EdgeColor', 'k', 'LineWidth', patchLineWidth );

hold on

% scatter(x(~bdyNodes,1), x(~bdyNodes,2), ...
%     [], 'filled', 'MarkerFaceColor', pntColor);

% scatter(circX(:,1), circX(:,2), 'filled', 'g');
% viscircles(x0, a);

hold off

axis equal tight

box off
grid off
axis off

% xlim(xBnd); ylim(yBnd)
xlim(x0(1) + windowSize * [-1 1] / 2);
ylim(x0(2) + windowSize * [-1 1] / 2);

xticks([]);
yticks([]);

set(gca,'LooseInset',get(gca,'TightInset'))

% Re-size figures for paper -----------------------------------------------
set(fig, 'Units', 'centimeters');

ax = gca;
ax.ActivePositionProperty = 'position';

ratio = ax.Position(4) ./ ax.Position(3);
fig.Position(3) = maxFigDim;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');


%% Visualize the Total Deformation ========================================
close all; clc;

% % Construct a set of circle test points
% ntheta = 100;
% theta = linspace(0, 2*pi, ntheta+1); theta(end) = [];
% 
% circX = zeros(ntheta, 2);
% [circX(:,1), circX(:,2)] = pol2cart(theta, a);
% circX = circX + x0;
% 
% circU = circularInclusionDisplacement(circX, x0, n, a, M, q, nu, thin3D);
%
% dispCircX = circX + circU;

% Calculate the displacement due to the inclusion
[u, F, majL, minL] = ...
    circularInclusionDisplacement(x, x0, n, a, M, q, nu, thin3D);

% Construct the motion of the daughter cells
xChild = repmat(x0, 2, 1);
uChild = F .* [ n; -n];

% Assemble of a list of displaced cell locations following division
allX = x; allX(divID, :) = []; allX = [allX; xChild];
allU = u; allU(divID, :) = []; allU = [allU; uChild];
allDispX = allX + allU;
dispU = allU; dispU(1:(end-2), :);
isChild = false(size(allU, 1), 1); isChild((end-1):end) = true;

% Generate a Voronoi tesselation of the displaced cells
dispDelTri = delaunayTriangulation(allDispX);
dispF = delTri.ConnectivityList;
[dispV, dispC] = voronoiDiagram(dispDelTri);

% Extract boundary cells
dispBdyFace = cellfun( @(x) ismember(1, x), dispC );

% Extract boundary nodes
dispBdyIDx = unique(freeBoundary(dispDelTri));
dispBdyNodes = ismember((1:size(allDispX,1)).', dispBdyIDx);

% Convert the voronoi cell connectivity list to a (NaN-padded)
% matrix Useful for plotting purposes
maxFaceSize = max(cellfun(@(x) numel(x), dispC));
dispVoronoiFace = nan(size(dispC,1), maxFaceSize);
for i = 1:size(dispC,1)
    dispVoronoiFace(i, 1:numel(dispC{i})) = dispC{i};
end

clear maxFaceSize

%--------------------------------------------------------------------------
% Generate the Visuzalization
%--------------------------------------------------------------------------

faceColors = ones( numel(dispBdyFace), 3 );
faceColors(isChild, :) = repmat(cellColor, 2, 1);

fig = figure('Color', [1 1 1]);

patch( 'Faces', dispVoronoiFace(~dispBdyFace, :), 'Vertices', dispV, ...
    'FaceVertexCData', faceColors(~dispBdyFace, :), ...
    'FaceColor', 'flat', 'EdgeColor', 'k', ...
    'LineWidth', patchLineWidth );

hold on

% scatter(allDispX(~dispBdyNodes,1), allDispX(~dispBdyNodes,2), ...
%     [], 'filled', 'MarkerFaceColor', pntColor);

% scatter(dispCircX(:,1), dispCircX(:,2), 'filled', 'g');
% plot_ellipse( x0(1), x0(2), majL/2, minL/2, atan2(n(2), n(1)), ...
%     'PlotMajor', true, 'PlotMinor', true, 'PlotEllipse', true );
% hold on

% quiver( allDispX(~dispBdyNodes & ~isChild,1), ...
%     allDispX(~dispBdyNodes & ~isChild,2), ...
%     arrowScaling_Total * allU(~dispBdyNodes & ~isChild,1), ...
%     arrowScaling_Total * allU(~dispBdyNodes & ~isChild,2), ...
%     0, 'MaxHeadSize', maxArrowHeadSize, ...
%     'LineWidth', arrowLineWidth, 'Color', arrowColor );

axis equal tight

% xlim(xBnd); ylim(yBnd)
xLim = x0(1) + windowSize * [-1 1] / 2;
yLim = x0(2) + windowSize * [-1 1] / 2;
xlim(xLim); ylim(yLim);
xticks([]); yticks([]);

% Generate the velocity arrows --------------------------------------------

% Find the start and stop points of the arrow
inWindow = (xLim(1) < allDispX(:,1)) & (allDispX(:,1) < xLim(2)) & ...
    (yLim(1) < allDispX(:,2)) & (allDispX(:,2) < yLim(2));
start = allDispX(~dispBdyNodes & ~isChild & inWindow, :);
stop = start + arrowScaling_Total * ...
    allU(~dispBdyNodes & ~isChild & inWindow, :);

% The length of the arrowheads in PIXELS
length = sqrt(sum(allU(~dispBdyNodes & ~isChild & inWindow, :).^2, 2));
length = 12 .* length / max(length);

baseAngle = 90 .* ones(size(length));
tipAngle = 30 .* ones(size(length));

% The width of the tail in PIXELS
width = sqrt(sum(allU(~dispBdyNodes & ~isChild & inWindow, :).^2, 2));
width = 4 .* width / max(width);

axis(axis)
arrow(start, stop, length, baseAngle, tipAngle, width, ...
    'FaceColor', arrowColor, 'EdgeColor', 'none', 'LineWidth', 0.5);

hold off

box off
grid off

set(gca,'LooseInset',get(gca,'TightInset'))
axis off

% Re-size figures for paper -----------------------------------------------
set(fig, 'Units', 'centimeters');

ax = gca;
ax.ActivePositionProperty = 'position';

ratio = ax.Position(4) ./ ax.Position(3);
fig.Position(3) = maxFigDim;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');


%% Visualize the Isotropic Deformation ====================================
close all; clc;

% % Construct a set of circle test points
% ntheta = 100;
% theta = linspace(0, 2*pi, ntheta+1); theta(end) = [];
% 
% circX = zeros(ntheta, 2);
% [circX(:,1), circX(:,2)] = pol2cart(theta, a);
% circX = circX + x0;
% 
% circU = circularInclusionDisplacement(circX, x0, n, a, M, 0, nu, thin3D);
% 
% dispCircX = circX + circU;

% Calculate the displacement due to the inclusion
[u, F, majL, minL] = ...
    circularInclusionDisplacement(x, x0, n, a, M, 0, nu, thin3D);

% Construct the motion of the daughter cells
xChild = repmat(x0, 2, 1);
uChild = F .* [ n; -n];

% Assemble of a list of displaced cell locations following division
allX = x; allX(divID, :) = []; allX = [allX; xChild];
allU = u; allU(divID, :) = []; allU = [allU; uChild];
allDispX = allX + allU;
dispU = allU; dispU(1:(end-2), :);
isChild = false(size(allU, 1), 1); isChild((end-1):end) = true;

% Generate a Voronoi tesselation of the displaced cells
dispDelTri = delaunayTriangulation(allDispX);
dispF = delTri.ConnectivityList;
[dispV, dispC] = voronoiDiagram(dispDelTri);

% Extract boundary cells
dispBdyFace = cellfun( @(x) ismember(1, x), dispC );

% Extract boundary nodes
dispBdyIDx = unique(freeBoundary(dispDelTri));
dispBdyNodes = ismember((1:size(allDispX,1)).', dispBdyIDx);

% Convert the voronoi cell connectivity list to a (NaN-padded)
% matrix Useful for plotting purposes
maxFaceSize = max(cellfun(@(x) numel(x), dispC));
dispVoronoiFace = nan(size(dispC,1), maxFaceSize);
for i = 1:size(dispC,1)
    dispVoronoiFace(i, 1:numel(dispC{i})) = dispC{i};
end

clear maxFaceSize

%--------------------------------------------------------------------------
% Generate the Visuzalization
%--------------------------------------------------------------------------

faceColors = ones( numel(dispBdyFace), 3 );
faceColors(isChild, :) = repmat(cellColor, 2, 1);

fig = figure('Color', [1 1 1]);

patch( 'Faces', dispVoronoiFace(~dispBdyFace, :), 'Vertices', dispV, ...
    'FaceVertexCData', faceColors(~dispBdyFace, :), ...
    'FaceColor', 'flat', 'EdgeColor', 'k', ...
    'LineWidth', patchLineWidth );

hold on

% scatter(allDispX(~dispBdyNodes,1), allDispX(~dispBdyNodes,2), ...
%     [], 'filled', 'MarkerFaceColor', pntColor);

% scatter(dispCircX(:,1), dispCircX(:,2), 'filled', 'g');
% plot_ellipse( x0(1), x0(2), majL/2, minL/2, atan2(n(2), n(1)), ...
%     'PlotMajor', true, 'PlotMinor', true, 'PlotEllipse', true );
% hold on

% quiver( allDispX(~dispBdyNodes & ~isChild,1), ...
%     allDispX(~dispBdyNodes & ~isChild,2), ...
%     arrowScaling_Iso * allU(~dispBdyNodes & ~isChild,1), ...
%     arrowScaling_Iso * allU(~dispBdyNodes & ~isChild,2), ...
%     0, 'LineWidth', arrowLineWidth, 'Color', arrowColor );

axis equal tight

% xlim(xBnd); ylim(yBnd)
xLim = x0(1) + windowSize * [-1 1] / 2;
yLim = x0(2) + windowSize * [-1 1] / 2;
xlim(xLim); ylim(yLim);
xticks([]); yticks([]);

% Generate the velocity arrows --------------------------------------------

% Find the start and stop points of the arrow
inWindow = (xLim(1) < allDispX(:,1)) & (allDispX(:,1) < xLim(2)) & ...
    (yLim(1) < allDispX(:,2)) & (allDispX(:,2) < yLim(2));
start = allDispX(~dispBdyNodes & ~isChild & inWindow, :);
stop = start + arrowScaling_Iso * ...
    allU(~dispBdyNodes & ~isChild & inWindow, :);

% The length of the arrowheads in PIXELS
length = sqrt(sum(allU(~dispBdyNodes & ~isChild & inWindow, :).^2, 2));
length = 12 .* length / max(length);

baseAngle = 90 .* ones(size(length));
tipAngle = 30 .* ones(size(length));

% The width of the tail in PIXELS
width = sqrt(sum(allU(~dispBdyNodes & ~isChild & inWindow, :).^2, 2));
width = 4 .* width / max(width);

axis(axis)
arrow(start, stop, length, baseAngle, tipAngle, width, ...
    'FaceColor', arrowColor, 'EdgeColor', 'none', 'LineWidth', 0.5);

hold off

box off
grid off
axis off

xlim(xLim); ylim(yLim);

set(gca,'LooseInset',get(gca,'TightInset'))

% Re-size figures for paper -----------------------------------------------
set(fig, 'Units', 'centimeters');

ax = gca;
ax.ActivePositionProperty = 'position';

ratio = ax.Position(4) ./ ax.Position(3);
fig.Position(3) = maxFigDim;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

%% Visualize the Deviatoric Deformation ===================================
close all; clc;

% % Construct a set of circle test points
% ntheta = 100;
% theta = linspace(0, 2*pi, ntheta+1); theta(end) = [];
% 
% circX = zeros(ntheta, 2);
% [circX(:,1), circX(:,2)] = pol2cart(theta, a);
% circX = circX + x0;
% 
% circU = circularInclusionDisplacement(circX, x0, n, a, 0, q, nu, thin3D);
% 
% dispCircX = circX + circU;

% Calculate the displacement due to the inclusion
[u, F, majL, minL] = ...
    circularInclusionDisplacement(x, x0, n, a, 0, q, nu, thin3D);

% Construct the motion of the daughter cells
xChild = repmat(x0, 2, 1);
uChild = F .* [ n; -n];

% Assemble of a list of displaced cell locations following division
allX = x; allX(divID, :) = []; allX = [allX; xChild];
allU = u; allU(divID, :) = []; allU = [allU; uChild];
allDispX = allX + allU;
dispU = allU; dispU(1:(end-2), :);
isChild = false(size(allU, 1), 1); isChild((end-1):end) = true;

% Generate a Voronoi tesselation of the displaced cells
dispDelTri = delaunayTriangulation(allDispX);
dispF = delTri.ConnectivityList;
[dispV, dispC] = voronoiDiagram(dispDelTri);

% Extract boundary cells
dispBdyFace = cellfun( @(x) ismember(1, x), dispC );

% Extract boundary nodes
dispBdyIDx = unique(freeBoundary(dispDelTri));
dispBdyNodes = ismember((1:size(allDispX,1)).', dispBdyIDx);

% Convert the voronoi cell connectivity list to a (NaN-padded)
% matrix Useful for plotting purposes
maxFaceSize = max(cellfun(@(x) numel(x), dispC));
dispVoronoiFace = nan(size(dispC,1), maxFaceSize);
for i = 1:size(dispC,1)
    dispVoronoiFace(i, 1:numel(dispC{i})) = dispC{i};
end

clear maxFaceSize

%--------------------------------------------------------------------------
% Generate the Visuzalization
%--------------------------------------------------------------------------

faceColors = ones( numel(dispBdyFace), 3 );
faceColors(isChild, :) = repmat(cellColor, 2, 1);

fig = figure('Color', [1 1 1]);

patch( 'Faces', dispVoronoiFace(~dispBdyFace, :), 'Vertices', dispV, ...
    'FaceVertexCData', faceColors(~dispBdyFace, :), ...
    'FaceColor', 'flat', 'EdgeColor', 'k', ...
    'LineWidth', patchLineWidth );

hold on

% scatter(allDispX(~dispBdyNodes,1), allDispX(~dispBdyNodes,2), ...
%     [], 'filled', 'MarkerFaceColor', pntColor);

% scatter(dispCircX(:,1), dispCircX(:,2), 'filled', 'g');
% plot_ellipse( x0(1), x0(2), majL/2, minL/2, atan2(n(2), n(1)), ...
%     'PlotMajor', true, 'PlotMinor', true, 'PlotEllipse', true );
% hold on

% quiver( allDispX(~dispBdyNodes & ~isChild,1), ...
%     allDispX(~dispBdyNodes & ~isChild,2), ...
%     arrowScaling_Ani * allU(~dispBdyNodes & ~isChild,1), ...
%     arrowScaling_Ani * allU(~dispBdyNodes & ~isChild,2), ...
%     0, 'LineWidth', arrowLineWidth, 'Color', arrowColor );

axis equal tight

% xlim(xBnd); ylim(yBnd)
xLim = x0(1) + windowSize * [-1 1] / 2;
yLim = x0(2) + windowSize * [-1 1] / 2;
xlim(xLim); ylim(yLim);

% Generate the velocity arrows --------------------------------------------

% Find the start and stop points of the arrow
inWindow = (xLim(1) < allDispX(:,1)) & (allDispX(:,1) < xLim(2)) & ...
    (yLim(1) < allDispX(:,2)) & (allDispX(:,2) < yLim(2));
start = allDispX(~dispBdyNodes & ~isChild & inWindow, :);
stop = start + arrowScaling_Ani * ...
    allU(~dispBdyNodes & ~isChild & inWindow, :);

% The length of the arrowheads in PIXELS
length = sqrt(sum(allU(~dispBdyNodes & ~isChild & inWindow, :).^2, 2));
length = 12 .* length / max(length);

baseAngle = 90 .* ones(size(length));
tipAngle = 30 .* ones(size(length));

% The width of the tail in PIXELS
width = sqrt(sum(allU(~dispBdyNodes & ~isChild & inWindow, :).^2, 2));
width = 4 .* width / max(width);

axis(axis)
arrow(start, stop, length, baseAngle, tipAngle, width, ...
    'FaceColor', arrowColor, 'EdgeColor', 'none', 'LineWidth', 0.5);


hold off

box off
grid off
axis off

set(gca,'LooseInset',get(gca,'TightInset'))

% Re-size figures for paper -----------------------------------------------
set(fig, 'Units', 'centimeters');

ax = gca;
ax.ActivePositionProperty = 'position';

ratio = ax.Position(4) ./ ax.Position(3);
fig.Position(3) = maxFigDim;
fig.Position(4) = ratio * fig.Position(3);

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');
