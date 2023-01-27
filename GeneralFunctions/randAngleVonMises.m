function rangles = randAngleVonMises(numSamples, varargin)
%RANDANGLEVONMISES Draw random samples from a circular von Mises
%distribution by rejection sampling
%
%   INPUT PARAMETERS:
%       
%       - numSamples:   The desired number of random samples from the
%                       distribution
%
%   (Name, Value)-Pairs:
%
%       - ('Mean' or 'Mu', mu == 0): The mean of the circular distribution
%
%       - ('Concentration' or 'Kappa', k == 1): A measure of the
%       concentration of the distribution
%
%       - ('NumPeaks', n == 1): The number of peaks in the distribution
%
%       - ('Domain', domainType == 'Pi'): Whether or not the basic domain
%       of the angles is [-pi, pi] or [0, 2*pi]
%
%       - ('TruncateDomain', truncateDomain == false): Whether or not to
%       truncate the basic domain. Only useful for n > 1. If this option is
%       set to true, then the samples will only be drawn from the single
%       peak centered around the mean within the basic domain
%
%       - ('PlotSamples', plotSamples == false): Produce a plot of the
%       samples drawn comparing them to the desired distribution. Mostly
%       for debugging
%
%   OUTPUT PARAMETERS:
%
%       - rangles:  The random samples from the distribution
%
%   by Dillon Cislo 2021/02/05

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------
if (nargin < 1), numSamples = 1; end

% Set default parameters
mu = 0;
k = 1;
n = 1;
domainType = 'pi';
truncateDomain = false;
plotSamples = false;

for i = 1:length(varargin)
    
    if isa(varargin{i}, 'double'), continue; end
    if isa(varargin{i}, 'logical'), continue; end
    
    if (strcmpi(varargin{i}, 'Mu') || strcmpi(varargin{i}, 'Mean'))
        mu = varargin{i+1};
    end
    if (strcmpi(varargin{i}, 'Kappa') || ...
            strcmpi(varargin{i}, 'Concentration'))
        k = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'NumPeaks')
        n = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'Domain')
        domainType = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'TruncateDomain')
        truncateDomain = varargin{i+1};
    end
    if strcmpi(varargin{i}, 'PlotSamples')
        plotSamples = varargin{i+1};
    end
    
end

validateattributes(numSamples, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite', 'real'});
validateattributes(mu, {'numeric'}, {'scalar', 'finite', 'real'});
validateattributes(k, {'numeric'}, ...
    {'scalar', 'finite', 'nonnegative', 'real'});
validateattributes(n, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite', 'real'});
validateattributes(truncateDomain, {'logical'}, {'scalar'});
validateattributes(plotSamples, {'logical'}, {'scalar'});


if ~(strcmpi(domainType, 'Pi') || strcmpi(domainType, '2Pi'))
    error('Invalid domain type specified')
end

% The normalized von Mises distribution over the full domain
vonMises = @(x, mu, k, n) exp(k .* cos(n*(x-mu))) ./ ...
        (2 * pi * besseli(0,k));

%--------------------------------------------------------------------------
% GENERATE RANDOM SAMPLES
%--------------------------------------------------------------------------

% Initial samples are generated over the [-pi, pi] regardless of final
% domain type or truncation
domain = [-pi, pi];
mu = wrapToPi(mu);

rangles = nan(numSamples, 1);
samplesGathered = 0;

while (samplesGathered < numSamples)
    
    % Draw a number uniformly from the domain of the distribtion
    curX = diff(domain) .* rand(1) + domain(1);
    
    % Draw a test sample from the proposed (uniform) distribution
    curY = rand(1);
    
    if (curY < vonMises(curX, mu, k, n))
        samplesGathered = samplesGathered + 1;
        rangles(samplesGathered) = curX;
    end
    
end

% Wrap the output to the appropriate domain
if truncateDomain
    
    % Shift the angles to be centered about the mean
    rangles = rangles - mu;
    
    % Divide the domain by the number of peaks
    % The domain should still be symmetric around zero
    domain = domain / n;
    
    % The total length of the truncated domain
    D = diff(domain);
    
    % Logical vector indexing angles outside the truncated domain
    outIDx = (rangles < domain(1)) | (rangles > domain(2));
    
    % Wrap the angles to [-D/2, D/2]
    ranglesTmp = rangles + (D/2);
    positiveInput = (ranglesTmp > 0);
    ranglesTmp = mod(ranglesTmp, D);
    ranglesTmp((ranglesTmp == 0) & positiveInput) = D;
    ranglesTmp = ranglesTmp - D/2;
    rangles(outIDx) = ranglesTmp(outIDx);
    
    % Shift the angles back
    rangles = rangles + mu;
    
end

if strcmpi(domainType, 'Pi')
    rangles = wrapToPi(rangles);
else
    rangles = wrapTo2Pi(rangles);
end

%--------------------------------------------------------------------------
% GENERATE A PLOT OF THE SAMPLES
%--------------------------------------------------------------------------

if plotSamples
    
    
    plotTheta = linspace(domain(1), domain(2), 200);
    
    if truncateDomain
        
        plotTheta = plotTheta + mu;
        
        vonMises = @(x, mu, k, n) n .* exp(k .* cos(n*(x-mu))) ./ ...
            (2 * pi * besseli(0,k));
        
    end
    
    if strcmpi(domainType, 'Pi')
        plotTheta = wrapToPi(plotTheta);
    else
        plotTheta = wrapTo2Pi(plotTheta);
    end
    
    plotR = vonMises(plotTheta, mu, k, n);
    
    insertNanIndex = [0 (diff(plotTheta) < 0)];
    insertValue = (1-insertNanIndex) ./ 0;
    
    plotTheta_tmp = [plotTheta, insertValue];
    plotTheta = plotTheta_tmp(:).';
    plotTheta(isinf(plotTheta)) = [];
    
    plotR_tmp = [plotR, insertValue];
    plotR = plotR_tmp(:).';
    plotR(isinf(plotR)) = [];
    
    polarhistogram(rangles, 'Normalization', 'pdf');
    
    hold on
    
    polarplot(plotTheta, plotR, '-r', 'LineWidth', 2);
    
    maxR = 1.2 .* max(plotR);
    polarplot( domain(1) * [1 1] + mu, [0, maxR], '-g', 'LineWidth', 2);
    polarplot( domain(2) * [1 1] + mu, [0, maxR], '-g', 'LineWidth', 2);
    
    hold off
    
    set(gca, 'RLim', [0, maxR]);
    
    
end


end

