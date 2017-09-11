function [array_sigs, src_array_sigs] = apply_source_signals_arrays(array_rirs, src_sigs)

nRec = length(array_rirs);
nSrc = size(array_rirs{1},3);
nSigs = size(src_sigs,2);
if (nSigs<nSrc), error('the number of source signals should be at least as many as the source number in the simulation'); end

array_sigs = cell(nRec,1);
src_array_sigs = cell(nRec,1);
L_sig = size(src_sigs,1);
for nr=1:nRec
    temprirs = array_rirs{nr};    
    L_rir = size(temprirs,1);
    nMics = size(temprirs,2);  
    tempsigs = zeros(L_sig+L_rir-1, nMics,nSrc);
    for ns=1:nSrc
disp(['Convolving with source signal: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
        tempsigs(:,:,ns) = fftfilt(temprirs(:,:,ns), [src_sigs(:,ns); zeros(L_rir-1,1)]);
    end
    src_array_sigs{nr} = tempsigs;
    array_sigs{nr} = sum(tempsigs,3);
end
