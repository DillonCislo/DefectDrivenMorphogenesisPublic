function [hEllipse, hMajor, hMinor] = ...
    plot_ellipse(xC, yC, a, b, phi, varargin)
%PLOT_ELLIPSE Plots a set of ellipses on a set of user specified axes
%
%   INPUT PARAMETERS:
%
%       - xC:       Nx1 list of x-coordinates of ellipse centers
%       - yC:       Nx1 list of y-coordinates of ellipse centers
%       - a:        Nx1 list of semi-major axis lengths
%       - b:        Nx1 list of semi-minor axis lengths
%       - phi:      Nx1 list of counter-clockwise rotation angles from the
%                   x-axis to the major axis
%
%       (Name, Value)-Pairs:
%
%           - 'Axis':         The axis on which to plot the ellipses {gca}
%           - 'NumPoints':    The number of points used to discretize
%                             each ellipse {200}
%           - 'PlotEllipse:   Plot the outline of the ellipses {true};
%           - 'PlotMajor':    Plot the major axis of each ellipse {false};
%           - 'PlotMinor':    Plot the minor axis of each ellipse {false};
%           - 'LineWidth':    The width of the plot lines {2};
%           - 'EllipseColor': The color of the ellipse outline
%                             {[0.8500, 0.3250, 0.0980]}
%           - 'MajorColor':   The color of the major axis {'m'}
%           - 'MinorColor':   The color of the minor axis {'c'}
%
%   OUTPUT PARAMETERS:
%
%       - hEllipse: The graphics handle for the ellipses
%
%   by Dillon Cislo 02/10/2020

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Validate mandatory inputs -----------------------------------------------
validateattributes( xC, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real'} );
validateattributes( yC, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real', 'numel', numel(xC)} );
validateattributes( a, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real', 'numel', numel(xC)} );
validateattributes( b, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real', 'numel', numel(xC)} );
validateattributes( phi, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real', 'numel', numel(xC)} );

% Make all inputs column vectors
if (size(xC, 2) ~= 1), xC = xC.'; end
if (size(yC, 2) ~= 1), yC = yC.'; end
if (size(a, 2) ~= 1), a = a.'; end
if (size(b, 2) ~= 1), b = b.'; end
if (size(phi, 2) ~= 1), phi = phi.'; end

% Process optional inputs -------------------------------------------------

plotAxis = gca;
numPoints = 200;
plotEllipse = true;
plotMajor = false;
plotMinor = false;
lwidth = 2;
ellipseColor = [0.8500, 0.3250, 0.0980];
majorColor = 'm';
minorColor = 'c';

for i = 1:length(varargin)
    if isa(varargin{i}, 'double') 
        continue;
    end
    if isa(varargin{i}, 'logical')
        continue;
    end
    if isa(varargin{i}, 'matlab.graphics.axis.Axes')
        continue;
    end
    if ~isempty(regexp(varargin{i}, '^[Aa]xis', 'match'))
        plotAxis = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Nn]um[Nn]eighbors', 'match'))
        numPoints = lower(varargin{i+1});
    end
    if ~isempty(regexp(varargin{i}, '^[Pp]lot[Mm]ajor', 'match'))
        plotMajor = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Pp]lot[Mm]inor', 'match'))
        plotMinor = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Pp]lot[Ee]llipse', 'match'))
        plotEllipse = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ll]ine[Ww]idth', 'match'))
        lwidth = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ee]llipse[Cc]olor', 'match'))
        ellipseColor = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]ajor[Cc]olor', 'match'))
        majorColor = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]inor[Cc]olor', 'match'))
        minorColor = varargin{i+1};
    end
end

assert( isa( plotAxis, 'matlab.graphics.axis.Axes' ), ...
    'Plot axes must be a handle to a MATLAB "Axes" object' );
validateattributes( numPoints, {'numeric'}, ...
    {'scalar', 'integer', 'finite', 'nonnan', 'real', 'positive'});

%--------------------------------------------------------------------------
% GENERATE DISCRETIZATION OF THE ELLIPSES
%--------------------------------------------------------------------------

% The angle of each point (1 x numPoints)
theta = 0:(2*pi)/(numPoints-1):(2*pi);

% The (x,y)-coordinates of each point presuming the major axis aligns with
% the x-axis

xx = repmat(a, 1, numPoints) .* repmat(cos(theta), numel(a), 1);
xx = xx.'; xx = xx(:);

yy = repmat(b, 1, numPoints) .* repmat(sin(theta), numel(b), 1);
yy = yy.'; yy = yy(:);

% Rotate the points to properly align the major axes
repPhi = repmat( phi.', numPoints, 1 ); repPhi = repPhi(:);

x = ( cos(repPhi) .* xx - sin(repPhi) .* yy);
y = ( sin(repPhi) .* xx + cos(repPhi) .* yy);

% Translate ellipse centers
repXC = repmat( xC.', numPoints, 1 ); repXC = repXC(:);
repYC = repmat( yC.', numPoints, 1 ); repYC = repYC(:);

x = x + repXC; y = y + repYC;

% Insert NaNs in between ellipses for plotting purposes
x = reshape( x, numPoints, numel(xC) );
x = [ x; nan( 1, numel(xC) ) ];
x = x(:);

y = reshape( y, numPoints, numel(yC) );
y = [ y; nan( 1, numel(yC) ) ];
y = y(:);

%--------------------------------------------------------------------------
% PLOT ELLIPSES
%--------------------------------------------------------------------------

if plotEllipse
    
    hEllipse = plot( x, y, 'LineWidth', lwidth, ...
        'Color', ellipseColor );
    
else
    
    hEllipse = [];
    
end

%--------------------------------------------------------------------------
% PLOT THE MAJOR AXIS
%--------------------------------------------------------------------------

if plotMajor
    
    majorX = [ a .* cos(phi), -a .* cos(phi) ] + [xC xC];
    majorX = [ majorX.'; nan(1, numel(a)) ];
    majorX = majorX(:);
    
    majorY = [ a .* sin(phi), -a .* sin(phi) ] + [yC yC];
    majorY = [ majorY.'; nan(1, numel(a)) ];
    majorY = majorY(:);
    
    hold on
    hMajor = plot( majorX, majorY, ...
        'LineWidth', lwidth, 'Color', majorColor );
    hold off
    
else
    
    hMajor = [];
    
end

%--------------------------------------------------------------------------
% PLOT THE MINOR AXIS
%--------------------------------------------------------------------------

if plotMinor
    
    minorX = [ -b .* sin(phi), b .* sin(phi) ] + [xC xC];
    minorX = [ minorX.'; nan(1, numel(b)) ];
    minorX = minorX(:);
    
    minorY = [ b .* cos(phi), -b .* cos(phi) ] + [yC yC];
    minorY = [ minorY.'; nan(1, numel(b)) ];
    minorY = minorY(:);
    
    hold on
    hMinor = plot( minorX, minorY, ...
        'LineWidth', lwidth, 'Color', minorColor );
    hold off
    
else
    
    hMinor = [];
    
end


end

