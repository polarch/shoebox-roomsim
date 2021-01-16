function [rt60, line_params] = revTimeRIR(H, fs)
%REVTIMERIR Computes the broadband reverberation time from RIRs
%   H: matrix of RIRs, one per column
%   line_params: line fitting parameters [[a1 b1]' [a2 b2]' ...], one line
%                per column/RIR

nRIR = size(H,2);
lRIR = size(H,1);

% line fitting parameters for RT estimation
offset_db = 5;
rt_interval_db = 30;
t = (0:size(H,1))/fs;
for nr=1:nRIR
    h = H(:,nr)/max(abs(H(:,nr)));
    % cumulative tail energy
%     h_flip = h(end:-1:1);
%     edc(1) = h_flip(1).^2;
%     for l=2:lRIR
%         edc(l) = edc(l-1) + h_flip(l).^2;
%     end
%     edc = edc(end:-1:1);
    edc = flipud(cumsum(flipud(h.^2))); % faster
    edc = edc/edc(1);    % normalized EDC
    edc_db = 10*log10(edc);
    % find line params y+ax+b
    [~,idx] = min(abs([edc_db+offset_db edc_db+(offset_db+rt_interval_db)]));
    y = edc_db(idx);
    x = t(idx);
    a = (y(2)-y(1))/(x(2)-x(1));
    b = y(1) - a*x(1);
    
    rt60(nr) = (-60-b)/a;
    line_params(:,nr) = [a b]';
end