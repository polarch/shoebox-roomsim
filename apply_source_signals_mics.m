function mic_sigs = apply_source_signals_mics(mic_rirs, src_sigs)

nRec = size(mic_rirs,2); 
nSrc = size(mic_rirs,3);
L_rir = size(mic_rirs,1);

nSigs = size(src_sigs,2);
if (nSigs<nSrc), error('the number of source signals should be at least as many as the source number in the simulation'); end
L_sig = size(src_sigs,1);

mic_sigs = zeros(L_sig+L_rir-1,nRec);
for nr=1:nRec
    for ns=1:nSrc
disp(['Convolving with source signal: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
        mic_sigs(:,nr) = mic_sigs(:,nr) + fftfilt(mic_rirs(:,nr,ns), [src_sigs(:,ns); zeros(L_rir-1,1)]);
    end
end