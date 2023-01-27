function divOrder = determine_division_order(XY, divOrderOptions)
%DETERMINE_DIVISION_ORDER Generates a queue for a set of input points.
%Placement in the queue corresponds to the order in which the cells will
%divide. Assumes that each cell divides once before any other cell divides
%twice
%
%   INPUT PARAMETERS:
%
%       - XY:               #Px2 list of input point coordinates
%
%       - divOrderOptions:	A struct whose fields control the output
%                           Default values in parenthesis
%
%       - orderType: The scheme used to determine the order in which the
%       cells divide
%           - ('wavefront'): Order is determined by the order in which the
%           cells collide with an elliptical wavefront
%           - 'random': Order is determined randomly
%
%       - wavefrontOrigin: The origin of the expanding elliptical wavefront
%       ([0 0])
%
%       - wavefrontVelocity: The (x,y)-components of the wavefront
%       velocity. Setting both components equal will yield a spherical
%       wavefront ([1 1])
%
%       - plotOrder:    Plot the division order (false)
%
%   OUTPUT PARAMETERS:
%
%       - divOrder:     #Px1 list of indices indicating the division order
%
%   by Dillon Cislo

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if (nargin < 1), error('Please supply input point list'); end
if (nargin < 2), divOrderOptions = struct(); end

validateattributes(XY, {'numeric'}, {'2d', 'ncols', 2, 'finite', 'real'});

% Check for invalid fields
fieldNames = {'orderType', 'wavefrontOrigin', 'wavefrontVelocity', ...
    'plotOrder'};
assert(all(ismember(fieldnames(divOrderOptions), fieldNames)), ...
    'Invalid division order options supplied');

if isfield(divOrderOptions, 'orderType')
    allTypes = {'wavefront', 'random'};
    assert(ismember(divOrderOptions.orderType, allTypes), ...
        'Invalid order type supplied');
else
    divOrderOptions.orderType = 'wavefront';
end

if isfield(divOrderOptions, 'wavefrontOrigin')
    validateattributes(divOrderOptions.wavefrontOrigin, {'numeric'}, ...
        {'vector', 'numel', 2, 'finite', 'real'});
else
    divOrderOptions.wavefrontOrigin = [0 0];
end

if isfield(divOrderOptions, 'wavefrontVelocity')
    validateattributes(divOrderOptions.wavefrontVelocity, {'numeric'}, ...
        {'vector', 'numel', 2, 'finite', 'real'});
else
    divOrderOptions.wavefrontVelocity = [1 1];
end

if isfield(divOrderOptions, 'plotOrder')
    validateattributes(divOrderOptions.plotOrder, {'logical'}, {'scalar'});
else
    divOrderOptions.plotOrder = false;
end

%--------------------------------------------------------------------------
% Generate Division Ordering
%--------------------------------------------------------------------------

switch divOrderOptions.orderType
    
    case 'wavefront'
        divOrder = generate_wavefront_ordering( XY, divOrderOptions );
        
    case 'random'
        divOrder = randperm(size(XY,1)).';
        
end

if divOrderOptions.plotOrder
    plot_div_order(XY, divOrder, divOrderOptions);
end

end

function divOrder = generate_wavefront_ordering( XY, divOrderOptions )

x0 = divOrderOptions.wavefrontOrigin;
vx = divOrderOptions.wavefrontVelocity(1);
vy = divOrderOptions.wavefrontVelocity(2);

r = XY - x0;

tCollide = sqrt( (r(:,1).^2 ./ vx^2) + (r(:,2).^2 / vy^2) );
[~, divOrder] = sort(tCollide);
[~, divOrder] = ismember((1:numel(divOrder)).', divOrder );

end

function plot_div_order( XY, divOrder, divOrderOptions )

divColors = parula(numel(divOrder));

scatter(XY(:,1), XY(:,2), [], divColors(divOrder, :), 'filled');
axis equal

end




