function [alpha_walls, rt60_true] = findAbsCoeffsFromRT(room, rt60_target, abs_wall_ratios)
%ROOM_STATS Return some room parameters from statistical room acoustics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ROOM_STATS.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    abs_wall_ratios = ones(1,6);
end

if max(abs_wall_ratios)~=1
    abs_wall_ratios = abs_wall_ratios/max(abs_wall_ratios);
end

nBands = length(rt60_target);
rt60_true = zeros(size(rt60_target));
alpha_walls = zeros(nBands,6);
for nb=1:nBands
    rt60 = rt60_target(nb);
    fmin = @(alpha) abs(rt60 - getRTsabine(alpha, room, abs_wall_ratios));
    [alpha, minval] = fminsearch(fmin, 0.0001);
    rt60_true(nb) = rt60 + minval;
    alpha_walls(nb,:) = min(alpha*abs_wall_ratios,1);
end

end

function rt60 = getRTsabine(alpha, room, abs_wall_ratios)

c = 343;
l = room(1);
w = room(2);
h = room(3);
V = l*w*h;
Stot = 2*(l*w + l*h + w*h);

alpha_walls = alpha*abs_wall_ratios;
a_x = alpha_walls(1:2);
a_y = alpha_walls(3:4);
a_z = alpha_walls(5:6);
% mean absorption
a_mean = sum(w*h*a_x+l*h*a_y+l*w*a_z)/Stot;
rt60 = (55.25*V)/(c*Stot*a_mean);

end
