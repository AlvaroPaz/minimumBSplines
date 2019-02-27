function Example_01
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
clc; clear all; close all;
display('Description: Minimum-effort control problem for RR robot')

%--------------------------------------------------------------------------
% Set variables for BSplines
global N n s B dB ddB S J L1 L2 h qi qf G0

N = 10;   n = 2;   
ti = 0;   h = 0.01;   tf = 1;   s = [ti:h:tf];

%--------------------------------------------------------------------------
% Build BSpline bases and its partials wrt s and c
[~,B,~,dB,~,ddB] = buildBasisFunctions(N,n,s);

%--------------------------------------------------------------------------
% Set kinematic and dynamic parameters

m1 = 1;   m2 = 5;   rx1 = 0.5;   rx2 = 0.5;   L1 = 1;   L2 = 0.6;

r1 = [rx1;0;0];
I1 = eye(3);
J1 = [m1*eye(3),-m1*skew(r1);m1*skew(r1),I1-m1*skew(r1)*skew(r1)];
S1 = [0,0,0,0,0,1]';

r2 = [rx2;0;0];
I2 = eye(3);
J2 = [m2*eye(3),-m2*skew(r2);m2*skew(r2),I2-m2*skew(r2)*skew(r2)];
S2 = [0,0,0,0,0,1]';

S(:,:,1) = S1;   S(:,:,2) = S2;   J(:,:,1) = J1;   J(:,:,2) = J2;

G0(:,:,1) = eye(4);   G0(1,4,1) = 0;
G0(:,:,2) = eye(4);   G0(1,4,2) = L1;

%--------------------------------------------------------------------------
% Set variables for optimization process
options = optimoptions('fmincon','display','iter','Algorithm','SQP','GradObj','on');

qi = [-pi/2;0];   qf = [pi*0.75;0.4];

c0 = randn(N*n,1)*0; % Decision variable initialization

% Initial and final boundaries
Aeq = [B(:,:,1);B(:,:,end)];
beq = [qi;qf];

% Initial and final velocity boundaries
dqi = zeros(2,1);   dqf = zeros(2,1);
Aeq = [Aeq;dB(:,:,1);dB(:,:,end)];
beq = [beq;dqi;dqf];

A = [];      b = [];

% Joint boundaries
lb = [];   ub = [];
q_low = ones(n,1)*-3;   q_up = ones(n,1)*3;
lb = [lb;kron(ones(N,1),q_low)];
ub = [ub;kron(ones(N,1),q_up)];

%--------------------------------------------------------------------------
% Call optimizer
tic;
c = fmincon(@evalCost,c0,A,b,Aeq,beq,lb,ub,[],options);
display(toc,'Computing time: ')

q = zeros(n,size(s,2));   dq = zeros(n,size(s,2));   ddq = zeros(n,size(s,2));
for i = 1:size(s,2)
    q(:,i) = B(:,:,i)*c;
    dq(:,i) = dB(:,:,i)*c;
    ddq(:,i) = ddB(:,:,i)*c;
end

[x] = computeForwardKinematics(q);

counter = 1;
[xf] = computeForwardKinematics(qf);
[xi] = computeForwardKinematics(qi);
[xi_op] = computeForwardKinematics(c(1:2));
figure();
for i = 1:size(x,2)
    aux = reshape(x(:,i),n,3);
    plot(aux(1,:),aux(2,:),'b','LineWidth',2.5); grid on; hold on;
    scatter(aux(1,:),aux(2,:),'o','MarkerFaceColor','r','MarkerEdgeColor','g');
    to_p = [x(5,1:i);x(6,1:i)];
    plot(to_p(1,:),to_p(2,:),':k','LineWidth',2); plot(xf(5),xf(6),'*k'); plot(xi(5),xi(6),'*k'); plot(xi_op(5),xi_op(6),'*r');
	axis equal; grid on;
    axis equal; grid on;
    xlabel('[m]', 'Interpreter', 'latex'); % x-y-axis label
    ylabel('[m]', 'Interpreter', 'latex');
    pause(0);
    counter = counter + 1;
    hold off;
end


function [x] = computeForwardKinematics(q)
global n L2 G0 S
x = zeros(3*n,size(q,2));   G = zeros(4,4,n);
G3 = eye(4);   G3(1,4) = L2;
for j = 1:size(q,2)
    for i = 1:n
        sk = skew(S(4:6,:,i));
        sk_2 = sk^2;
        e_wq = eye(3) + sin(q(i,j))*sk + (1-cos(q(i,j)))*sk_2;
        T_wq = q(i,j)*eye(3) + (1-cos(q(i,j)))*sk + (q(i,j)-sin(q(i,j)))*sk_2;
        G(:,:,i) = G0(:,:,i)*[e_wq,T_wq*S(1:3,:,i);zeros(1,3),1];
    end
    A0 = G(:,:,1);    A1 = G(:,:,1)*G(:,:,2);    A2 = G(:,:,1)*G(:,:,2)*G3;
    x(:,j) = [A0(1:2,4);A1(1:2,4);A2(1:2,4)];
end


function [f,g] = evalCost(c)
global N n s B dB ddB S J G0 h

%--------------------------------------------------------------------------
% Compute objective function
f = 0;   g = zeros(size(c));
for k = 1:size(s,2)
    %----------------------------------------------------------------------
    % Eval k-th q, dq and ddq
    q = B(:,:,k)*c;
    dq = dB(:,:,k)*c;
    ddq = ddB(:,:,k)*c;
    
    %----------------------------------------------------------------------
    % Inverse Dynamics with Geoemtric Newton-Euler
    Z = zeros(6,6,n);   V = zeros(6,1,n+1);   dV = zeros(6,1,n+1);   G = zeros(4,4,n);
    dV(:,:,1) = [0;9.81;0;0;0;0]; % Acceleration due gravity
    for i = 1:1:n
        sk = skew(S(4:6,:,i));
        sk_2 = sk^2;
        e_wq = eye(3) + sin(q(i))*sk + (1-cos(q(i)))*sk_2;
        T_wq = q(i)*eye(3) + (1-cos(q(i)))*sk + (q(i)-sin(q(i)))*sk_2;
        G(:,:,i) = G0(:,:,i)*[e_wq,T_wq*S(1:3,:,i);zeros(1,3),1];
        
        Z(:,:,i) = Ad(G(:,:,i));
        V(:,:,i+1) = Z(:,:,i)*V(:,:,i) + S(:,:,i)*dq(i);
        dV(:,:,i+1) = Z(:,:,i)*dV(:,:,i) + S(:,:,i)*ddq(i) + omega(V(:,:,i+1))*S(:,:,i)*dq(i);
    end

    Z_dual = zeros(6,6,n);   F = zeros(6,1,n+1);   Tau = zeros(n,1);   ad_dual_V = zeros(6,6,n);
    for i = n:-1:1
        ad_dual_V(:,:,i) = omega_dual(V(:,:,i+1));
        if i == n
            F(:,:,i) = J(:,:,i)*dV(:,:,i+1) - ad_dual_V(:,:,i)*J(:,:,i)*V(:,:,i+1);
        else
            Z_dual(:,:,i) = Ad_dual(G(:,:,i+1));
            F(:,:,i) = Z_dual(:,:,i)*F(:,:,i+1) + J(:,:,i)*dV(:,:,i+1) - ad_dual_V(:,:,i)*J(:,:,i)*V(:,:,i+1);
        end
        Tau(i,1) = S(:,:,i)'*F(:,:,i);
    end
    f = f + 0.5*h*(Tau'*Tau);

    %----------------------------------------------------------------------
    % Compute required gradient
    if nargout >= 2
        %----------------------------------------------------------------
        % First Derivative of Inverse Dynamics wrt c

        D_q = B(:,:,k);    D_dq = dB(:,:,k);    D_ddq = ddB(:,:,k);

        % D_V, D_dV and D_F have n+1 elements due the initialization
        D_V = zeros(6,N*n,n+1);   D_dV = zeros(6,N*n,n+1);   adS = zeros(6,6,n);
        for i = 1:1:n
            adS(:,:,i) = omega(S(:,:,i));
            D_V(:,:,i+1) = Z(:,:,i)*D_V(:,:,i) - adS(:,:,i)*V(:,:,i+1)*D_q(i,:) + S(:,:,i)*D_dq(i,:);
            D_dV(:,:,i+1) = Z(:,:,i)*D_dV(:,:,i) + S(:,:,i)*D_ddq(i,:) - adS(:,:,i)*( Z(:,:,i)*dV(:,:,i)*D_q(i,:) + D_V(:,:,i+1)*dq(i) + V(:,:,i+1)*D_dq(i,:) );
        end

        D_F = zeros(6,N*n,n+1);   D_Tau = zeros(n,N*n);
        for i = n:-1:1
            aux_I = zeros(6*N*n,6);
            for ii = 1:1:N*n
                aux_I(1+6*(ii-1):6*ii,:) = omega_dual(D_V(:,ii,i+1));
            end
            aux_I = aux_I*J(:,:,i)*V(:,:,i+1);
            aux_I = reshape(aux_I,6,N*n);
            if i == n
                D_F(:,:,i) = J(:,:,i)*D_dV(:,:,i+1) - aux_I - ad_dual_V(:,:,i)*J(:,:,i)*D_V(:,:,i+1);
            else
                D_F(:,:,i) = Z_dual(:,:,i)*( omega_dual(-S(:,:,i+1))*F(:,:,i+1)*D_q(i+1,:) + D_F(:,:,i+1) ) + J(:,:,i)*D_dV(:,:,i+1) - aux_I - ad_dual_V(:,:,i)*J(:,:,i)*D_V(:,:,i+1);
            end
            D_Tau(i,:) = S(:,:,i)'*D_F(:,:,i);
        end
        g = g + h*(D_Tau'*Tau);
    end
end
