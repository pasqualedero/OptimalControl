function [F,Qe] = constrainedRHC(A,B,C,x0,u_max,y_max)
% Compute the shaping matrix Qe of terminal set and gain F, given: a
% polytopic model, euclidean norm input and output constraints  

[n,m]=size(B{1});
p = size(C,1);
Q1=eye(n);
R=eye(m);
yalmip('clear')
gamma=sdpvar;
Q=sdpvar(n);
Y=sdpvar(m,n,'full');
con=[Q>=0, [1 x0';x0 Q]>=0];

% for each vertex of the polytope
for j=1:length(A)
    Aj=A{j};
    Bj=B{j};
    % LMI of robust unconstrained RHC + LMI on output constraints
    con=[con, [Q Q*Aj'+Y'*Bj' Q*sqrt(Q1) Y'*sqrt(R);
                Aj*Q+Bj*Y Q zeros(n) zeros(n,m);
                sqrt(Q1)*Q zeros(n) gamma*eye(n) zeros(n,m);
                sqrt(R)*Y zeros(m,n) zeros(m,n) gamma*eye(m)]>=0, ...
                [Q (Aj*Q+Bj*Y)'*C';
                C*(Aj*Q+Bj*Y) y_max^2*eye(p)]>=0];
end

% LMI of input constraints
con = [con, [u_max^2*eye(m) Y; Y' Q]>=0];

% Define solver options
options = sdpsettings('solver', 'sedumi', 'verbose', 0);
% Solve the optimization problem
sol = optimize(con, gamma, options);
% Check if the problem is solved
if sol.problem == 0
    Qe = inv(value(Q));
    F=value(Y)*Qe;
else
    disp(sol.info);
    error('Something went wrong:');
end
end

