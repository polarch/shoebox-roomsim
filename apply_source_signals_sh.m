function [sh_sigs, src_sh_sigs] = apply_source_signals_sh(sh_rirs, src_sigs)

nRec = size(sh_rirs,3);
nSrc = size(sh_rirs,4);
nSH = size(sh_rirs,2);
L_rir = size(sh_rirs,1);

nSigs = size(src_sigs,2);
if (nSigs<nSrc), error('the number of source signals should be at least as many as the source number in the simulation'); end
L_sig = size(src_sigs,1);

sh_sigs = zeros(L_sig+L_rir-1,nSH,nRec);
src_sh_sigs = zeros(L_sig+L_rir-1,nSH,nRec,nSrc);
for nr=1:nRec
    for ns=1:nSrc
        disp(['Convolving with source signal: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
        idx_nonzero = find(sum(squeeze(sh_rirs(:,:,nr,ns)))~=0);
        src_sh_sigs(:,idx_nonzero,nr,ns) = fftfilt(sh_rirs(:,idx_nonzero,nr,ns), [src_sigs(:,ns); zeros(L_rir-1,1)]);
%        sh_sigs(:,idx_nonzero,nr) = sh_sigs(:,idx_nonzero,nr) + src_sh_sigs(:,idx_nonzero,nr,ns);
    end
end
sh_sigs = sum(src_sh_sigs,4);