%% SETUP
datapath = '~/Downloads/';
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
rec = [ 3.4 2.1 1.7];
nRec = size(rec,1);

% source start
xs0 = [ 7.4 2.1 1.8];
% source velocity (m/s)
s_vel = 1;
% other point in line
xs1 = [ 3.4 3.1 1.8];
% duration of source signal (5sec)
tSig = 5;
% block length for RIR interpolation (50msec)
tBlock = 0.05;
% number of blocks
nBlocks = ceil(tSig/tBlock);
% get source positions for block RIRs
nSrc = nBlocks+1;
src = zeros(nSrc,3);
for ns=1:nSrc
    src(ns,:) = xs0 + (xs1-xs0)/sqrt(sum((xs1-xs0).^2))*(ns-1)*tBlock;
end

%% DEFINE DIRECTIONAL IRs FOR RECEIVERS
fs = 48000;

% Receiver 3 example: Measured HRTFs (MAKE SURE THAT THEY ARE THE SAME
% SAMPLERATE AS THE REST OR RESAMPLE)
disp('Loading measured HRTF responses')
load('ownsurround_sim2016.mat','hrtf_dirs','hrtf_mtx')
grids{1} = hrtf_dirs*pi/180;
array_irs{1} = hrtf_mtx;
array_irs{1} = array_irs{1}(250+(1:256),:,:);
clear hrtf_mtx hrtf_dirs

%% RUN SIMULATOR %%%
%%% Echogram
tic

% limit the RIR by reflection order or by time-limit
type = 'maxTime';
maxlim = 0.5; % just cut if it's longer than that ( or set to max(rt60) )
for nb = 1:nBands
    if (rt60(nb)<maxlim) limits(nb) = rt60(nb);
    else limits(nb,1) = maxlim;
    end
end
% compute echograms
[abs_echograms, echograms] = compute_echograms_arrays(room, src, rec, abs_wall, limits);

%% Rendering
% In the array case, for each receiver position, a grid of directions
% should be provided, along with the corresponding multichannel IRs for the
% elements of the array at the grid directions. 
array_rirs = render_array_rirs(abs_echograms, band_centerfreqs, fs, grids, array_irs);

toc

%% Generate sound scenes
% Each source is convolved with the respective array IRs, and summed with
% the rest of the sources to create the microphone mixed signals

sourcepath = '4src_samples_voice_handclaps_fountain_piano.wav';
[src_sigs, fs_src] = audioread(sourcepath);
if (fs_src~=fs), resample(src_sigs, fs, fs_src); end % resample if it doesn't match the project's fs

array_sigs = apply_source_signals(array_rirs, src_sigs);
