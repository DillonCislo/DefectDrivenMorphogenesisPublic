function [u, focusDist, majorL, minorL] = ...
    circularInclusionDisplacement(x, x0, n, a, M, q, nu, thin3D)
%CIRCULARINCLUSIONDISPLACEMENT Calculate the displacement field of an
%elastic body due to a circular Eshelby inclusion with uniform, but
%otherwise general, eigenstrain. Note that this displacement maps the
%originally circular inclusion into an ellipse.
%
%   INPUT PARAMETERS:
%
%       - x:        #Px2 list of query point coordinates
%       - x0:       1x2 location of the center of the inclusion
%       - n:        1x2 unit vector along the axis of elongation
%                   (Only meaningful for plastic deformations with
%                   non-vanishing deviatoric parts)
%       - a:        The radius of the circular inclusion
%       - M:        The strength of the isotropic part of the eigenstrain
%       - q:        The strength of the deviatoric part of the eigenstrain
%       - nu:       Poisson's ratio
%                   -1 <= nu <= 1 for true 2D bodies
%                   -1 <= nu <= 1/2 for thin 3D bodies
%       - thin3D:   Whether to treat the medium as a true 2D body or a thin
%                   3D body
%
%   OUTPUT PARAMETERS:
%
%       - u:            #Px2 list of displacements for each query point
%       - focusDist:    The distance of the deformed ellipse's foci to the
%                       centroid of the ellipse
%       - majorL:       The length of the major axis of the ellipse
%       - minorL:       The length of the minor axis of the ellipse
%
%   by Dillon Cislo 02/02/2021

% Reformat input variables for vectorized calculations
nVec = repmat(n, size(x, 1), 1);
RVec = x - x0;
R2 = sum(RVec.^2, 2);

if thin3D
    
    % Calculate the displacement within the inclusion
    u_in = M .* RVec ./ (4 * (1-nu)) + ...
        q * (3-4*nu)/(4*(1-nu)) * ...
        ( 2 .* dot(nVec, RVec, 2) .* nVec - RVec );
    
    % Calculate the displacment outside the inclusion
    u_out = q/(4*(1-nu)) .* (a^2 ./ R2) .* ...
        (2*(1-2*nu) + (a^2 ./ R2)) .* ...
        (2 .* dot(nVec, RVec, 2) .* nVec - RVec);
    u_out = u_out + ...
        M/(4*(1-nu)) .* (a^2 ./ R2).^2 .* RVec;
    u_out = u_out + ...
        1 / (2*(1-nu)) * (a^2 ./ R2) .* ...
        (1-(a^2 ./ R2)) .* ...
        (M/2 + q*(2 .* dot(nVec, RVec, 2).^2 ./ R2 - 1)) ...
        .* RVec;
    
    % Calculate the distance from either focus to the ellipse center
    focusDist = sqrt((a^2/(4*(1-nu)^2)) * (4*(1-nu)+M) * (3-4*nu) * q);
    
    % Calculate the major axis length
    majorL = sqrt( (a^2/(4*(1-nu)^2)) * ( 4*(1-nu) + M + q*(3-4*nu) )^2 );
    
    % Calculate the minor axis length
    minorL = sqrt( (a^2/(4*(1-nu)^2)) * ( 4*(1-nu) + M - q*(3-4*nu) )^2 );
    
else
    
    % Calculate the displacement within the inclusion
    u_in = M .* (1+nu) .* RVec ./ 4 + ...
        q * ((3-nu)/4) * ...
        ( 2 .* dot(nVec, RVec, 2) .* nVec - RVec );
    
    % Calculate the displacment outside the inclusion
    u_out = q * ((1+nu)/4) .* (a^2 ./ R2) .* ...
        (2*(1-nu)/(1+nu) + (a^2 ./ R2)) .* ...
        (2 .* dot(nVec, RVec, 2) .* nVec - RVec);
    u_out = u_out + ...
        M * ((1+nu)/4) .* (a^2 ./ R2).^2 .* RVec;
    u_out = u_out + ...
        ((1+nu)/2) * (a^2 ./ R2) .* ...
        (1-(a^2 ./ R2)) .* ...
        (M/2 + q*(2 .* dot(nVec, RVec, 2).^2 ./ R2 - 1)) ...
        .* RVec;
    
    % Calculate the distance from either focus to the ellipse center
    focusDist = sqrt( (a^2/4) * (4+M*(1+nu)) * (3-nu) * q);
    
    % Calculate the major axis length
    majorL = sqrt( (a^2/4) * ( 4 + M*(1+nu) + q*(3-nu) )^2 );
    
    % Calculate the minor axis length
    minorL = sqrt( (a^2/4) * ( 4 + M*(1+nu) - q*(3-nu) )^2 );
    
end

% Associate displacements to the geometric subregions
u = u_out;
u(sqrt(R2) <= a, :) = u_in(sqrt(R2) <= a, :); 

end

