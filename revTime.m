function rt60 = revTime(room, alpha, method)
%REVTIME Computes reverberation time of a rectangular room using various methods.

    L = room(1);
    W = room(2);
    H = room(3);

    switch method
        case 'sabine'
            V = L*W*H;
            Sx = W*H;
            Sy = L*H;
            Sz = L*W;
            A = Sx*(alpha(1) + alpha(2)) + Sy*(alpha(3) + alpha(4)) + ...
                Sz*(alpha(5) + alpha(6));
            rt60 = 0.161*(V/A);

        case 'norriseyring'
            V = L*W*H;
            Sx = W*H;
            Sy = L*H;
            Sz = L*W;
            Stot = 2*Sx + 2*Sy + 2*Sz;
            A = Sx*(alpha(1) + alpha(2)) + Sy*(alpha(3) + alpha(4)) + ...
                Sz*(alpha(5) + alpha(6));
            Am = A/Stot;
            rt60 = 0.161*V/(-Stot*log(1-Am));

        case 'millingtonsette'
            V = L*W*H;
            Sx = W*H;
            Sy = L*H;
            Sz = L*W;
            A  = -(Sx*(log(1-alpha(1)) + log(1-alpha(2))) + ...
                   Sy*(log(1-alpha(3)) + log(1-alpha(4))) + ...
                   Sz*(log(1-alpha(5)) + log(1-alpha(6))));
            rt60 = 0.161*(V/A);

        case 'fitzroy'
            V = L*W*H;
            Sx = W*H;
            Sy = L*H;
            Sz = L*W;
            Stot = 2*Sx + 2*Sy + 2*Sz;        
            tx = -2*Sx/log(1-mean(alpha(1:2)));
            ty = -2*Sy/log(1-mean(alpha(3:4)));
            tz = -2*Sz/log(1-mean(alpha(5:6)));
            rt60 = 0.161*V/Stot^2 * (tx+ty+tz);

        case 'arau'
            V = L*W*H;
            Sx = W*H;
            Sy = L*H;
            Sz = L*W;
            Stot = 2*Sx + 2*Sy + 2*Sz;
            Tx = (0.161*V/(-Stot * log(1-mean(alpha(1:2)))))^(2*Sx/Stot);
            Ty = (0.161*V/(-Stot * log(1-mean(alpha(3:4)))))^(2*Sy/Stot);
            Tz = (0.161*V/(-Stot * log(1-mean(alpha(5:6)))))^(2*Sz/Stot);
            rt60 = Tx*Ty*Tz;

        case 'neubauerkostek'
            V = L*W*H;
            Sx = W*H;
            Sy = L*H;
            Sz = L*W;
            Stot = 2*Sx + 2*Sy + 2*Sz;        
            R   = 1 - alpha;
            Rxy = mean(R(1:4));
            Rz = mean(R(5:6));
            Rm  = mean(R);
            Axy = log(1/Rm) + (R(1)*(R(1)-Rxy)*Sx^2 + ...
                               R(2)*(R(2)-Rxy)*Sx^2 + ...
                               R(3)*(R(3)-Rxy)*Sy^2 + ...
                               R(4)*(R(4)-Rxy)*Sy^2)/ ...
                               ((Rxy*(2*Sx+2*Sy))^2);
            Az = log(1/Rm) + (R(5)*(R(5)-Rz)*Sz^2 + ...
                              R(6)*(R(6)-Rz)*Sz^2)/ ...
                              ((Rz*2*Sz)^2);
            rt60 = 0.32*V/(Stot^2) * (H*(L+W)/Axy + L*W/Az);
    end

end
