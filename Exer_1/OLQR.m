function [PP,FF]=OLQR(A,B,Q,R,S,N)
    
    PP = {S}; FF = {};

    for k = N-1:-1:0

        lastP = PP{length(PP)};

        % Calcolo F(k)
        F = -(R + B' * lastP * B)^-1 * B' * lastP * A;
        FF = [FF F];

        if k == 0
            break;
        end

        % Calcolo P(k)
        P = A' * lastP * A + Q - A' * lastP * B * (R + B' * lastP * B)^-1 * B' * lastP * A;
        PP = [PP P];

    end

    PP = flip(PP);
    FF = flip(FF);

end
