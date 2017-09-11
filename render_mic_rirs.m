function rirs = render_mic_rirs(echograms, band_centerfreqs, fs)
%
% echograms nSrc nRec nBands echogram structure
nRec = size(echograms,2);
nSrc = size(echograms,1);
nBands = size(echograms,3);

% sample echogram to a specific sampling rate with fractional interpolation
FRACTIONAL = 1;
% decide on number of samples for all RIRs
endtime = 0;
for ns=1:nSrc
    for nr=1:nRec
        temptime = echograms(ns, nr).time(end);
        if (temptime>endtime), endtime = temptime; end
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
rirs = zeros(L_tot, nr, ns);
for ns=1:nSrc
    for nr=1:nRec
        disp('')
        disp(['Rendering echogram: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
        tempIR = zeros(L_rir, nBands);
        for nb=1:nBands
            tempIR(:,nb) = render_rir(echograms(ns,nr,nb), endtime, fs, FRACTIONAL);
        end
        disp('     Filtering and combining bands')
        rirs(:,nr,ns) = filter_rir(tempIR, band_centerfreqs, fs);
    end
end
