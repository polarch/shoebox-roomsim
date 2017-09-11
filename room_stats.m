function [RT60_sabine,  d_critical, d_mfpath] = room_stats(room, abs_wall)
%ROOM_STATS Return some room parameters from statistical room acoustics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ROOM_STATS.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c = 343;

l = room(1);
w = room(2);
h = room(3);

% total volume
V = l*w*h;

% total room wall area
Stot = 2*(l*w + l*h + w*h);

% analyse in frequency bands
Nbands = size(abs_wall,1);
for m = 1:Nbands
    a_x = abs_wall(m,1:2);
    a_y = abs_wall(m,3:4);
    a_z = abs_wall(m,5:6);
    % mean absorption
    a_mean(m) = sum(w*h*a_x+l*h*a_y+l*w*a_z)/Stot;


    RT60_sabine(m) = (55.25*V)/(c*Stot*a_mean(m));
end

d_critical = 0.1*sqrt(V./(pi*RT60_sabine));

d_mfpath = 4*V/Stot;

disp(['Room dimensions (m)          ' num2str(l) 'x' num2str(w) 'x' num2str(h)])
disp(['Room volume (m^3)            ' num2str(V)])
disp(['Mean absorption coeff        ' num2str(a_mean,2)])
disp(['Sabine Rev. Time 60dB (sec)  ' num2str(RT60_sabine,3)])
disp(['Critical distance (m)        ' num2str(d_critical,2)])
disp(['Mean free path (m)           ' num2str(d_mfpath,2)])
end
