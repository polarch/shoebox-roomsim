function h = lagrange(N, delays)
%LAGRANGE  h=lagrange(N,delay) returns order N FIR 
%          filter h which implements given delay 
%          (in samples).  For best results, 
%          delay should be near N/2 +/- 1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LAGRANGE.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = (0:N)';
h = ones(N+1, length(delays));
for l=1:length(delays)
    for k = 0:N
        index = find(n ~= k);
        h(index, l) = h(index, l) *  (delays(l)-k)./ (n(index)-k);
    end
end

end
