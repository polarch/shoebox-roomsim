%% SETUP
datapath = '/home/politis/Devel/spatial_acoustic_analysis/data/';
addpath(datapath);

% room definition
room = [10.2 7.1 3.2];

% desired RT per octave band, and time to truncate the responses
rt60 = [1.0 0.8 0.7 0.6 0.5 0.4];
nBands = length(rt60);
% lowest octave band
band_centerfreqs(1) = 125;
for (nb=2:nBands) band_centerfreqs(nb) = 2*band_centerfreqs(nb-1); end
% absorption for approximately achieving the RT60 above - row per band
abs_wall = findAbsCoeffsFromRT(room, rt60);
% critical distance for the room
[~,  d_critical] = room_stats(room, abs_wall);

% receiver position
rec = [ 4.5 3.4 1.5;
        2.0 3.1 1.4];
nRec = size(rec,1);

% source positions
src = [ 6.2 2.0 1.8;
        7.9 3.3 1.75;
        5.8 5.0 1.9];
nSrc = size(src,1);

% SH orders for receivers
rec_orders = [1 3]; % rec1: first order(4ch), rec2: 3rd order (16ch)

% % convert source directions from listener-centric to room-centric
% [src_coords(:,1), src_coords(:,2), src_coords(:,3)] = sph2cart(src_dirs(:,1)*pi/180, ...
%     src_dirs(:,2)*pi/180, src_r);
% src_coords(:,2) = -src_coords(:,2);
% src = ones(nSrc,1)*rec(1,:) + src_coords;
% % check sources
% for n=1:nSrc
%     if (src(n,1)>room(1))||(src(n,2)>room(2))||(src(n,3)>room(3))
%         error('Source coordinates out of room boundaries')
%     end
% end

%%% RUN SIMULATOR %%%
%% Echogram
tic

maxlim = 1.5; % just stop if the echogram goes beyond that time ( or just set it to max(rt60) )
for nb = 1:nBands
    if (rt60(nb)<maxlim) limits(nb) = rt60(nb);
    else limits(nb,1) = maxlim;
    end
end
% compute echograms
[abs_echograms, rec_echograms, echograms] = compute_echograms_sh(room, src, rec, abs_wall, limits, rec_orders);

%% Rendering
% In this case all the information (e.g. SH directivities) are already
% encoded in the echograms, hence they are rendered directly to discrete RIRs
fs = 48000;
sh_rirs = render_sh_rirs(abs_echograms, band_centerfreqs, fs);

toc

%% Generate sound scenes
% Each source is convolved with the respective mic IR, and summed with
% the rest of the sources to create the microphone mixed signals

sourcepath = '4src_samples_voice_handclaps_fountain_piano.wav';
[src_sigs, fs_src] = audioread(sourcepath);
if (fs_src~=fs), resample(src_sigs, fs, fs_src); end % resample if it doesn't match the project's fs

sh_sigs = apply_source_signals_sh(sh_rirs, src_sigs);