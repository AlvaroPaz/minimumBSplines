function Example_03
clc; clear all; close all;
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
display('Description: Minimum-effort control problem for biped robot walker')

%--------------------------------------------------------------------------
% Set variables for BSplines
global N n s B dB ddB S J L1 L2 L3 L4 L5 h qi qf son parent G0 LL

N = 7;    n = 5;
ti = 0;   h = 0.01;   tf = 1;   s = linspace(ti,tf,10);

%--------------------------------------------------------------------------
% Build BSpline bases and its partials wrt s and c
[~,B,~,dB,~,ddB] = buildBasisFunctions(N,n,s);

%--------------------------------------------------------------------------
% Set kinematic and dynamic parameters
parent = [0 1 2 2 4];
son = {2,[3 4],-1,5,-1};

link_len = 0.5;

L1 = link_len;      L2 = link_len;      L3 = link_len;      L4 = link_len;      L5 = link_len;
rx1 = link_len/2;   rx2 = link_len/2;   rx3 = link_len/2;   rx4 = link_len/2;   rx5 = link_len/2;
m1 = 1;             m2 = 1;             m3 = 1;             m4 = 1;             m5 = 1;

LL = [L1,L2,L3,L4,L5];

r1 = [rx1;0;0];
I1 = eye(3);
J1 = [m1*eye(3),-m1*skew(r1);m1*skew(r1),I1-m1*skew(r1)*skew(r1)];
S1 = [0,0,0,0,0,1]';

r2 = [rx2;0;0];
I2 = eye(3);
J2 = [m2*eye(3),-m2*skew(r2);m2*skew(r2),I2-m2*skew(r2)*skew(r2)];
S2 = [0,0,0,0,0,1]';

r3 = [rx3;0;0];
I3 = eye(3);
J3 = [m3*eye(3),-m3*skew(r3);m3*skew(r3),I3-m3*skew(r3)*skew(r3)];
S3 = [0,0,0,0,0,1]';

r4 = [rx4;0;0];
I4 = eye(3);
J4 = [m4*eye(3),-m4*skew(r4);m4*skew(r4),I4-m4*skew(r4)*skew(r4)];
S4 = [0,0,0,0,0,1]';

r5 = [rx5;0;0];
I5 = eye(3);
J5 = [m5*eye(3),-m5*skew(r5);m5*skew(r5),I5-m5*skew(r5)*skew(r5)];
S5 = [0,0,0,0,0,1]';

S(:,:,1) = S1;   S(:,:,2) = S2;   S(:,:,3) = S3;   S(:,:,4) = S4;   S(:,:,5) = S5;
J(:,:,1) = J1;   J(:,:,2) = J2;   J(:,:,3) = J3;   J(:,:,4) = J4;   J(:,:,5) = J5;

G0(:,:,1) = eye(4);   G0(1,4,1) = 0;
G0(:,:,2) = eye(4);   G0(1,4,2) = L1;
G0(:,:,3) = eye(4);   G0(1,4,3) = L2;
G0(:,:,4) = eye(4);   G0(1,4,4) = L2;
G0(:,:,5) = eye(4);   G0(1,4,5) = L4;

%--------------------------------------------------------------------------
% Set variables for optimization process
options = optimoptions('fmincon','Display','iter','Algorithm','SQP','GradObj','on');

%Flat terrain
qi = [(pi/2)-0.1;0.7;-0.8;1+(pi/2);-0.6285];
qf = [0.9715;0.6285;-0.2292;2*pi-1-(pi/2);-0.7];

% % Upstair
% qi = [1.25;1.27;-1.13;2.26;-0.43];
% qf = [1.2;0.43;-0.25;4.03;-1.27];

c0 = randn(N*n,1)*0; % Decision variable initialization

% Initial and final boundaries
Aeq = [B(:,:,1);B(:,:,end)];
beq = [qi;qf];

% Initial and final velocity boundaries
% dqi = zeros(n,1);   dqf = zeros(n,1);
% Aeq = [Aeq;dB(:,:,1);dB(:,:,end)];
% beq = [beq;dqi;dqf];

A = [];      b = [];

% Joint boundaries
lb = [];   ub = [];
% q_low = ones(n,1)*-3;   q_up = ones(n,1)*3;
% lb = [lb;kron(ones(N,1),q_low)];
% ub = [ub;kron(ones(N,1),q_up)];

%--------------------------------------------------------------------------
% Call optimizer
tic;
c = fmincon(@evalCost,c0,A,b,Aeq,beq,lb,ub,[],options);
display(toc,'Computing time')

q = zeros(n,size(s,2));   dq = zeros(n,size(s,2));   ddq = zeros(n,size(s,2));
for i = 1:1:size(s,2)
    q(:,i) = B(:,:,i)*c;
    dq(:,i) = dB(:,:,i)*c;
    ddq(:,i) = ddB(:,:,i)*c;
end

%--------------------------------------------------------------------------
% Plotting
[Rr,Pr] = computeForwardKinematics(q);

counter = 1;
figure();
for i = 1:size(Pr,3)
    
    r_ = 0.025;   r__ = 0.01;   np = 30;
    
    k = 3; counter2 = 1;
    for j = [1,2,3,3,5]
        l_ = LL(counter2);
        t_ = linspace(pi/2,3*(pi/2),np);
        aux_1 = [r_*cos(t_);r_*sin(t_)]; % left arc
        aux_2 = [linspace(0,l_,np);-r_*ones(1,np)]; % lower line
        aux_3 = [-r_*cos(t_)+l_;-r_*sin(t_)]; % right arc
        aux_4 = [linspace(l_,0,np);r_*ones(1,np)]; % upper line  
    
        R = Rr(:,k:k+1,1);   P = Pr(:,j,1);
        
        aux_5 = [aux_1,aux_2,aux_3,aux_4];
        aux_5 = R*aux_5;   aux_5(1,:) = aux_5(1,:)+P(1);   aux_5(2,:) = aux_5(2,:)+P(2);
        
        fill(aux_5(1,:),aux_5(2,:),'c'); hold on;
        plot(aux_5(1,:),aux_5(2,:),'-m','LineWidth',1); axis equal;
        
        t_ = linspace(0,2*pi,2*np);
        aux_6 = [r__*cos(t_);r__*sin(t_)]; % left joint
        aux_6 = R*aux_6;   aux_6(1,:) = aux_6(1,:)+P(1);   aux_6(2,:) = aux_6(2,:)+P(2);
        aux_7 = [r__*cos(t_)+l_;r__*sin(t_)];
        aux_7 = R*aux_7;   aux_7(1,:) = aux_7(1,:)+P(1);   aux_7(2,:) = aux_7(2,:)+P(2);
        fill(aux_6(1,:),aux_6(2,:),'m');
        fill(aux_7(1,:),aux_7(2,:),'m');
        k = k + 2;
        counter2 = counter2 + 1;
    end
    
    k = 3;   counter2 = 1;
    for j = [1,2,3,3,5]
        l_ = LL(counter2);
        t_ = linspace(pi/2,3*(pi/2),np);
        aux_1 = [r_*cos(t_);r_*sin(t_)]; % left arc
        aux_2 = [linspace(0,l_,np);-r_*ones(1,np)]; % lower line
        aux_3 = [-r_*cos(t_)+l_;-r_*sin(t_)]; % right arc
        aux_4 = [linspace(l_,0,np);r_*ones(1,np)]; % upper line  
        R = Rr(:,k:k+1,i);   P = Pr(:,j,i);
        
        aux_5 = [aux_1,aux_2,aux_3,aux_4];
        aux_5 = R*aux_5;   aux_5(1,:) = aux_5(1,:)+P(1);   aux_5(2,:) = aux_5(2,:)+P(2);
        
        fill(aux_5(1,:),aux_5(2,:),'g'); hold on;
        plot(aux_5(1,:),aux_5(2,:),'-r','LineWidth',1); axis equal;
        
        t_ = linspace(0,2*pi,2*np);
        aux_6 = [r__*cos(t_);r__*sin(t_)]; % left joint
        aux_6 = R*aux_6;   aux_6(1,:) = aux_6(1,:)+P(1);   aux_6(2,:) = aux_6(2,:)+P(2);
        aux_7 = [r__*cos(t_)+l_;r__*sin(t_)];
        aux_7 = R*aux_7;   aux_7(1,:) = aux_7(1,:)+P(1);   aux_7(2,:) = aux_7(2,:)+P(2);
        fill(aux_6(1,:),aux_6(2,:),'r');
        fill(aux_7(1,:),aux_7(2,:),'r');
        k = k + 2;
        
        for j_ = 1:n+1
            aux_8 = Pr(:,j_,1:i);   aux_8 = reshape(aux_8,2,[],1);
            plot(aux_8(1,:),aux_8(2,:),':k','LineWidth',1.5);
        end
        counter2 = counter2 + 1;
    end
    
    str = ['s = ',num2str(s(i)),'\rightarrow'];
    text(-0.7,1.2,str)

    axis equal; grid on; axis([-0.835 1 -0.235 1.6]);
    xlabel('\textbf{x} [m]', 'Interpreter', 'latex'); % x-y-axis label
    ylabel('\textbf{y} [m]', 'Interpreter', 'latex');
    set(gcf,'units','points','position',[100,100,500,400])
    pause(0);
    counter = counter + 1;
    hold off;
end


function [R,P] = computeForwardKinematics(q)
global n L1 L2 L3 L4 L5
R = zeros(2,2*(n+1),size(q,2));
P = zeros(2,n+1,size(q,2));

for i = 1:size(q,2)
    A0  = eye(4);              Mt0 = A0;
    R(:,1:2,i) = Mt0(1:2,1:2);   P(:,1,i) = Mt0(1:2,4);
    A1  = DH(q(1,i),0,L1,0);   Mt1 = Mt0*A1;
    R(:,3:4,i) = Mt1(1:2,1:2);   P(:,2,i) = Mt1(1:2,4);
    A2  = DH(q(2,i),0,L2,0);   Mt2 = Mt1*A2;
    R(:,5:6,i) = Mt2(1:2,1:2);   P(:,3,i) = Mt2(1:2,4);
    A3  = DH(q(3,i),0,L3,0);   Mt3 = Mt2*A3;
    R(:,7:8,i) = Mt3(1:2,1:2);   P(:,4,i) = Mt3(1:2,4);
    
    A4  = DH(q(4,i),0,L4,0);   Mt4 = Mt2*A4;
    R(:,9:10,i) = Mt4(1:2,1:2);   P(:,5,i) = Mt4(1:2,4);
    A5  = DH(q(5,i),0,L5,0);   Mt5 = Mt4*A5;
    R(:,11:12,i) = Mt5(1:2,1:2);   P(:,6,i) = Mt5(1:2,4);
end


function [f,g] = evalCost(c)
global N n s B dB ddB S J h son parent G0
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
    % Inverse Dynamics with Recursive Newton-Euler
    % Comment: V, dV and F have n+1 elements due the initialization
    Z = zeros(6,6,n);   V = zeros(6,1,n+1);   dV = zeros(6,1,n+1);   G = zeros(4,4,n);
    dV(:,:,1) = [0;9.81;0;0;0;0];
    for i = 1:1:n
        j = parent(i)+1;
        
        sk = skew(S(4:6,:,i));
        sk_2 = sk^2;
        e_wq = eye(3) + sin(q(i))*sk + (1-cos(q(i)))*sk_2;
        T_wq = q(i)*eye(3) + (1-cos(q(i)))*sk + (q(i)-sin(q(i)))*sk_2;
        G(:,:,i) = G0(:,:,i)*[e_wq,T_wq*S(1:3,:,i);zeros(1,3),1];
        
        Z(:,:,i) = Ad(G(:,:,i));
        V(:,:,i+1) = Z(:,:,i)*V(:,:,j) + S(:,:,i)*dq(i);
        dV(:,:,i+1) = Z(:,:,i)*dV(:,:,j) + S(:,:,i)*ddq(i) + omega(V(:,:,i+1))*S(:,:,i)*dq(i);
    end

    Z_dual = zeros(6,6,n);   F = zeros(6,1,n+1);   Tau = zeros(n,1);   ad_dual_V = zeros(6,6,n);
    for i = n:-1:1
        ad_dual_V(:,:,i) = omega_dual(V(:,:,i+1));
        F(:,:,i) = J(:,:,i)*dV(:,:,i+1) - ad_dual_V(:,:,i)*J(:,:,i)*V(:,:,i+1);
        %------------------------------------------------------------------
        if son{i} ~= -1        
            for j = son{i}
                Z_dual(:,:,j) = Ad_dual(G(:,:,j)); 
                F(:,:,i) = F(:,:,i) + Z_dual(:,:,j)*F(:,:,j);
            end
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
        for i = 1:n
            j = parent(i)+1;
            adS(:,:,i) = omega(S(:,:,i));
            D_V(:,:,i+1) = Z(:,:,i)*D_V(:,:,j) - adS(:,:,i)*V(:,:,i+1)*D_q(i,:) + S(:,:,i)*D_dq(i,:);
            D_dV(:,:,i+1) = Z(:,:,i)*D_dV(:,:,j) + S(:,:,i)*D_ddq(i,:) - adS(:,:,i)*( Z(:,:,i)*dV(:,:,j)*D_q(i,:) + D_V(:,:,i+1)*dq(i) + V(:,:,i+1)*D_dq(i,:) );
        end

        D_F = zeros(6,N*n,n+1);   D_Tau = zeros(n,N*n);
        for i = n:-1:1    
            %--------------------------------------------------------------
            aux_I = zeros(6*N*n,6);
            for i_ = 1:1:N*n
                aux_I(1+6*(i_-1):6*i_,:) = omega_dual(D_V(:,i_,i+1));
            end
            aux_I = aux_I*J(:,:,i)*V(:,:,i+1);
            aux_I = reshape(aux_I,6,N*n);
            %--------------------------------------------------------------
            D_F(:,:,i) = J(:,:,i)*D_dV(:,:,i+1) - aux_I - ad_dual_V(:,:,i)*J(:,:,i)*D_V(:,:,i+1);
            %--------------------------------------------------------------
            if son{i} ~= -1
                for j = son{i}
                    D_F(:,:,i) = D_F(:,:,i) + Z_dual(:,:,j)*( omega_dual(-S(:,:,j))*F(:,:,j)*D_q(j,:) + D_F(:,:,j) );
                end
            end
            D_Tau(i,:) = S(:,:,i)'*D_F(:,:,i);
        end
        g = g + h*(D_Tau'*Tau);
    end
end

