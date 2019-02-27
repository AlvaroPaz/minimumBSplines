%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
clc; clear all; close all;
startup % Be sure you have the Featherstone toolbox for comparing
display('Description: Validation of Newton-Euler geometric algorithms comparing with Featherstone and finite differences')

t = linspace(0,5,500);

ct = 1;
Q1 = sin(t);      Q2 = -ct*sin(t);    Q3 = ct*sin(t);      Q4 = -sin(t);    Q5 = sin(t);
dQ1 = cos(t);     dQ2 = -ct*cos(t);   dQ3 = ct*cos(t);     dQ4 = -cos(t);   dQ5 = cos(t);
ddQ1 = -sin(t);   ddQ2 = ct*sin(t);   ddQ3 = -ct*sin(t);   ddQ4 = sin(t);   ddQ5 = -sin(t);

n = 5;
parent = [0 1 2 1 4];

%--------------------------------------------------------------------------
% Featherstone toolbox
robot = planar(5);
robot.gravity = [0;-9.81];
robot.parent = parent;

for i = 1:n
    robot.I{i} = mcI(1,[0.5;0],1);
end
h = 2^-17;


m1 = 1;      m2 = 1;      m3 = 1;      m4 = 1;      m5 = 1;
rx1 = 0.5;   rx2 = 0.5;   rx3 = 0.5;   rx4 = 0.5;   rx5 = 0.5;
L1 = 1;      L2 = 1;      L3 = 1;      L4 = 1;      L5 = 1;

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
G0(:,:,4) = eye(4);   G0(1,4,4) = L1;
G0(:,:,5) = eye(4);   G0(1,4,5) = L4;

G = zeros(4,4,n);

TAU_t = zeros(n,1,size(t,2));       TAU_f = zeros(n,1,size(t,2));
D_TAU_t = zeros(n,2*n,size(t,2));   D_TAU_f = zeros(n,2*n,size(t,2));

time_1 = 0;   time_2 = 0;
for k = 1:size(t,2)
    q1 = Q1(k);       q2 = Q2(k);       q3 = Q3(k);          q4 = Q4(k);          q5 = Q5(k);
    dq1 = dQ1(k);     dq2 = dQ2(k);     dq3 = dQ3(k);        dq4 = dQ4(k);        dq5 = dQ5(k);
    ddq1 = ddQ1(k);   ddq2 = ddQ2(k);   ddq3 = ddQ3(k);      ddq4 = ddQ4(k);      ddq5 = ddQ5(k);

    s1 = sin(q1);   s2 = sin(q2);   s3 = sin(q3);   s4 = sin(q4);   s5 = sin(q5);
    c1 = cos(q1);   c2 = cos(q2);   c3 = cos(q3);   c4 = cos(q4);   c5 = cos(q5);
    q = [q1;q2;q3;q4;q5];   dq = [dq1;dq2;dq3;dq4;dq5];   ddq = [ddq1;ddq2;ddq3;ddq4;ddq5];
    
    %----------------------------------------------------------------------
    % Inverse Dynamics with Featherstone toolbox
    tic;
    TAU_f(:,:,k) = ID( robot, q, dq, ddq, {} );
    
    %----------------------------------------------------------------------
	% First Derivative of Featherstone using finite differences
    
    for i_ = 1:n
        q_h = q;    q_h(i_) = q_h(i_) + h;
        tau_h = ID( robot, q_h, dq, ddq, {} );
        D_TAU_f(:,i_,k) = (tau_h - TAU_f(:,:,k))/h;
        
        dq_h = dq;    dq_h(i_) = dq_h(i_) + h;
        tau_h = ID( robot, q, dq_h, ddq, {} );
        D_TAU_f(:,i_+n,k) = (tau_h - TAU_f(:,:,k))/h;
    end
    
    time_1 = time_1 + toc;
    
    %----------------------------------------------------------------------
    % Inverse Dynamics with Geometric Algorithms
    tic;
    % V, dV and F have n+1 elements due the initialization
    Z = zeros(6,6,n);   V = zeros(6,1,n+1);   dV = zeros(6,1,n+1);
    dV(:,:,1) = [0;9.81;0;0;0;0];
    for i = 1:n
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
        for j = find(parent==i)
            Z_dual(:,:,j) = Ad_dual(G(:,:,j)); 
            F(:,:,i) = F(:,:,i) + Z_dual(:,:,j)*F(:,:,j);
        end
        Tau(i,1) = S(:,:,i)'*F(:,:,i);
    end
    TAU_t(:,:,k) = Tau;    
    
    %----------------------------------------------------------------------
    % First Derivative of Inverse Dynamics w.r.t. the state x = [q;dq]
    
    D_X = eye(2*n);

    D_q = D_X(1:n,:);   D_dq = D_X(n+1:2*n,:);   D_ddq = zeros(n,2*n);

    % D_V, D_dV and D_F have n+1 elements due the initialization
    D_V = zeros(6,2*n,n+1);   D_dV = zeros(6,2*n,n+1);   adS = zeros(6,6,n);
    for i = 1:1:n
        j = parent(i)+1;
        adS(:,:,i) = omega(S(:,:,i));
        D_V(:,:,i+1) = Z(:,:,i)*D_V(:,:,j) - adS(:,:,i)*V(:,:,i+1)*D_q(i,:) + S(:,:,i)*D_dq(i,:);
        D_dV(:,:,i+1) = Z(:,:,i)*D_dV(:,:,j) - adS(:,:,i)*( Z(:,:,i)*dV(:,:,j)*D_q(i,:) + D_V(:,:,i+1)*dq(i) + V(:,:,i+1)*D_dq(i,:) );% + S(:,:,i)*D_ddq(i,:);
    end

    D_F = zeros(6,2*n,n+1);   D_Tau = zeros(n,2*n);
    for i = n:-1:1    
        
        aux_I = zeros(6*2*n,6);
        for i_ = 1:1:2*n
            aux_I(1+6*(i_-1):6*i_,:) = omega_dual(D_V(:,i_,i+1));
        end
        aux_I = aux_I*J(:,:,i)*V(:,:,i+1);
        aux_I = reshape(aux_I,6,2*n);
        
        D_F(:,:,i) = J(:,:,i)*D_dV(:,:,i+1) - aux_I - ad_dual_V(:,:,i)*J(:,:,i)*D_V(:,:,i+1);
        
        for j = find(parent==i)
            D_F(:,:,i) = D_F(:,:,i) + Z_dual(:,:,j)*( omega_dual(-S(:,:,j))*F(:,:,j)*D_q(j,:) + D_F(:,:,j) );
        end
        D_Tau(i,:) = S(:,:,i)'*D_F(:,:,i);
    end
    D_TAU_t(:,:,k) = D_Tau;
    
    time_2 = time_2 + toc;
end

display(time_1,'Computing time for numerical solution')
display(time_2,'Computing time for analytical solution')

display('In all graphics blue represents Geometric algorithms and green Featherstone algorithms')

%----------------------------------------------------------------------
% Plot Inverse Dynamics
figure(n+1);
for i = 1:n
    plot(t,reshape(TAU_f(i,1,:),1,[],1),'-g','LineWidth',1.5); grid on; hold on;
    plot(t,reshape(TAU_t(i,1,:),1,[],1),':b','LineWidth',1);
end
xlabel('Time [s]', 'Interpreter', 'latex'); % x-y-axis label
ylabel('Tau', 'Interpreter', 'latex');
title('Inverse Dynamics');

%----------------------------------------------------------------------
% Plot Differential Inverse Dynamics
for i = 1:n
    figure(i);
    for j = 1:2*n
        subplot(2,n,j);
        plot(t,reshape(D_TAU_f(i,j,:),1,[],1),'-g','LineWidth',1.5); grid on; hold on;
        plot(t,reshape(D_TAU_t(i,j,:),1,[],1),':b','LineWidth',1.5);
        xlabel('Time [s]', 'Interpreter', 'latex'); % x-y-axis label
        ylabel('$\nabla\tau$', 'Interpreter', 'latex');
        axis tight
    end
end

display(sum(sum(sum(abs(D_TAU_f-D_TAU_t)))),'Cumulative error between analytical al finite differentiation: ');
