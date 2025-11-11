classdef OC_util
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        MC; % a Struct containing Q, R, M and S; 
        plant; % a Struct containing plant matrices A, B; 
        FF; PP; P; F;
    end


    methods
        % Construct an instance of this class
        function obj = OC_util(Q,R,M,S,A,B)
           
            obj.MC.Q = Q;
            obj.MC.R = R;
            obj.MC.M = M;
            obj.MC.S = S;

            obj.plant.A = A;
            obj.plant.B = B;

        end

        function [PP, FF] = OC_finite_horizon(obj,N)
            [PP, FF] = OLQR(obj.plant.A, obj.plant.B, obj.MC.Q, obj.MC.R, obj.MC.S, N);
        end

        function [x, u] = simulation(obj, N, x_0, finite_hor)
            x = [x_0];
            u = [];
            for k = 1:N-1    
                u_k = obj.FF{k} * x(:,k);
                u = [u u_k];
                x_k_forc = obj.plant.A * x(:,k) + obj.plant.B * u_k;
                x = [x x_k_forc];
            end
            u = [u obj.FF{N} * x(:,N)];
        end

    end
end