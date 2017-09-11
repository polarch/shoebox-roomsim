 function [decay] = backward_int(reflections, DRAW)
%BACKWARD_INT Summary of this function goes here
%   Detailed explanation goes here

% If DRAW is true, plot echograms
if nargin <= 1
    DRAW = 0;
end

% find number of frequency bands
Nbands = length(reflections);

for m = 1:Nbands
    
    % backward integrate IR
    t = reflections(m).time;
    x = reflections(m).value;
    
    if t(1) ~= 0
        x(2:end + 1) = x;
        x(1) = 0;
        t(2:end + 1) = t;
        t(1) = 0;
    end    
    
    S_back = zeros(1,length(x));
    
    % integration of squared pressure (power)
%    S_back(end:-1:1) = abs(cumtrapz(t(end:-1:1), x(end:-1:1).^2));
%    S_tot = trapz(t, x.^2);
    
    % integration of absolute pressure
    S_back(end:-1:1) = abs(cumtrapz(t(end:-1:1), x(end:-1:1)));
    S_tot = trapz(t, x);    
    
    EDC = 10*log10(S_back/S_tot);
    

    decay(m).EDC = EDC;
    decay(m).time = t;
    
    t0 = t(2);
    
    %linear interpolation for -5dB
    idx1 = find(decay(m).EDC>-5, 1, 'last' );
    idx2 = idx1 +1;
    
    t5 = t(idx1) + (-5 - decay(m).EDC(idx1))*(t(idx2) - t(idx1))/ ...
        (decay(m).EDC(idx2) - decay(m).EDC(idx1));
    
    %linear interpolation for -10dB
    idx1 = find(decay(m).EDC>-10, 1, 'last' );
    idx2 = idx1 +1;
    
    t10 = t(idx1) + (-10 - decay(m).EDC(idx1))*(t(idx2) - t(idx1))/ ...
        (decay(m).EDC(idx2) - decay(m).EDC(idx1));
    
    EDT = t10 - t0;
    
    %linear interpolation for -35dB
    idx1 = find(decay(m).EDC>-35, 1, 'last' );
    idx2 = idx1 +1;
    
    t35 = t(idx1) + (-35 - decay(m).EDC(idx1))*(t(idx2) - t(idx1))/ ...
        (decay(m).EDC(idx2) - decay(m).EDC(idx1));
    
    RT30 = t35 - t5;
    
    decay(m).RT60 = 2*RT30;
    decay(m).EDT = EDT;

    a_EDT = ((-10)/(t10-t0))/10;
    a_RT30 = ((-35+5)/(t35-t5))/10;
    decay(m).slope = [a_EDT a_RT30];
    decay(m).points = [t0 t5 t10 t35];
    
end

if DRAW
    figure
    for m=1:Nbands
        subplot(ceil(Nbands/2),2,m)
        t = decay(m).time;
        plot(t, decay(m).EDC)
        hold on
        plot(decay(m).points([1 3]), [0 -10],'g','LineWidth',1)
        plot(decay(m).points([2 4]), [-5 -35],'r','LineWidth',1)      
        xlabel('Time (sec)')
        ylabel('Relative Amplitude')
        title(['Backward Integrated IR ' num2str(125*2^(m-1)) 'Hz'])
        text(0.8, -10, ['RT60 = ' num2str(decay(m).RT60,2) 'sec'])
        text(0.8, -5, ['EDT = ' num2str(decay(m).EDT,2) 'sec'])
        axis([0 1 -80 0])
        grid
    end
end

end