function rirs = render_mic_rirs(echograms, band_centerfreqs, fs, endtime)
%
% echograms nSrc nRec nBands echogram structure
nRec = size(echograms,2);
nSrc = size(echograms,1);
nBands = size(echograms,3);

% sample echogram to a specific sampling rate with fractional interpolation
FRACTIONAL = 1;
% randomize slightly position of image source to avoid sweeping echoes
RAND_IMS = 1;
% decide on number of samples for all RIRs
if nargin<4 || isempty(endtime)
    endtime = 0;
    for ns=1:nSrc
        for nr=1:nRec
            temptime = echograms(ns, nr).time(end);
            if (temptime>endtime), endtime = temptime; end
        end
    end
end
L_rir = ceil(endtime*fs);
if (nBands>1)
    L_fbank = 1000;
else
    L_fbank = 0;
end
L_tot = L_rir + L_fbank;



% render responses and apply filterbank to combine different decays at different bands
rirs = zeros(L_tot, nRec, nSrc);
for ns=1:nSrc
    for nr=1:nRec
        % if randomization/jitter of image source positions is on, precompute
        % random time delays uniformly distributed on -dx:dx, where dx is the
        % maximum displacement along the source to receiver ray, from the true
        % image position (here 10cm)
        if RAND_IMS
            dx_max = 0.1;
            dt_max = dx_max/343; % speed of sound
            dt = 2*dt_max*rand(length(echograms(ns,nr,1).time),1)-dt_max;
            dt(1) = 0; % preserve timing of direct sound
        else
            dt = 0;
        end
        disp('')
        disp(['Rendering echogram: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
        tempIR = zeros(L_rir, nBands);
        for nb=1:nBands
            tempIR(:,nb) = render_rir(echograms(ns,nr,nb), endtime, fs, FRACTIONAL, RAND_IMS, dt);
        end
        disp('     Filtering and combining bands')
        rirs(:,nr,ns) = filter_rir(tempIR, band_centerfreqs, fs);
    end
end
