function array_rirs = render_array_rirs(echograms, band_centerfreqs, fs, grids, array_irs)

% echograms nSrc nRec nBands echogram structure
nRec = size(echograms,2);
nSrc = size(echograms,1);
nBands = size(echograms,3);
if (nRec~=length(grids))||(nRec~=length(array_irs))
    error('The number of grids and sets of array responses should match the receiver positions')
end

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
if (nBands>1),  L_fbank = 1000;
else            L_fbank = 0; end

array_rirs = cell(nr,1);
for nr=1:nRec
    grid_dirs_rad = grids{nr};
    nGrid= size(grid_dirs_rad,1);    
    mic_irs = array_irs{nr}; % should be nSamples x nMics x nGrid
    L_resp = size(mic_irs,1);
    nMics = size(mic_irs,2);
    array_rirs{nr} = zeros(L_rir+L_fbank+L_resp-1, nMics, nSrc);
    for ns=1:nSrc
disp(['Rendering echogram: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
disp('      Quantize echograms to receiver grid')
        echo2gridMap = get_echo2gridMap(echograms(ns,nr), grid_dirs_rad);
        tempRIR = zeros(L_rir, nGrid, nBands);
        for nb = 1:nBands
            % First step: reflections are quantized to the grid directions
            q_echograms = quantise_echogram(echograms(ns, nr, nb), nGrid, echo2gridMap);
            % Second step: render quantized echograms
disp(['     Rendering quantized echograms: Band ' num2str(nb)])
            tempRIR(:,:,nb) = render_quantized(q_echograms, endtime, fs, FRACTIONAL);
        end
        tempRIR2 = zeros(L_rir+L_fbank, nGrid);
disp('      Filtering and combining bands')
        for ng=1:nGrid        
            tempRIR2(:,ng) = filter_rir(squeeze(tempRIR(:,ng,:)), band_centerfreqs, fs);
        end
        clear tempRIR
        
        % Third step: convolve with directional IRs at grid directions
        idx_nonzero = find(sum(tempRIR2.^2)>10e-12); % neglect grid directions with almost no energy
disp('      Convolving with receiver responses')
        tempRIR2 = [tempRIR2(:, idx_nonzero); zeros(L_resp-1, length(idx_nonzero))];
        for nm = 1:nMics
            tempResp = squeeze(mic_irs(:, nm, idx_nonzero));
            array_rirs{nr}(:, nm, ns) = sum(fftfilt(tempResp, tempRIR2), 2);
        end
    end
end
