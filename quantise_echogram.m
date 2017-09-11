function q_echogram = quantise_echogram(echogram, nGrid, echo2gridMap)
%QUANTISE_ECHOGRAM Quantises the reflections to specific rendering directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% QUANTISE_ECHOGRAM.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise echogram structure with each quantised echogram
if length(echo2gridMap)>length(echogram.time)
    echo2gridMap = echo2gridMap(1:length(echogram.time));
end
q_echogram(nGrid).value = [];
q_echogram(nGrid).time = [];
for n=1:nGrid
    q_echogram(n).value = echogram.value(echo2gridMap == n);
    q_echogram(n).time = echogram.time(echo2gridMap == n);
    if isempty(q_echogram(n).time), q_echogram(n).isActive = 0;
    else                            q_echogram(n).isActive = 1; end
end

end
