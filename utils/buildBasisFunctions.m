function [b,B,db,dB,ddb,ddB] = buildBasisFunctions(N,n,s)
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
% Function to compute B(s) and its derivatives w.r.t. s

% N = number of control points
% s = vector of time parameter
% k = index of enabled control points
% m = grade --> default m = 3
% j = iterator through the knots
% n = DoF

m = 3;

kn = linspace(s(1),s(end),N-m+1);
t=[s(1)*ones(1,m),kn,s(end)*ones(1,m)];

b = zeros(1,N,size(s,2));     B = zeros(n,n*N,size(s,2));
db = zeros(1,N,size(s,2));    dB = zeros(n,n*N,size(s,2));
ddb = zeros(1,N,size(s,2));   ddB = zeros(n,n*N,size(s,2));
j = 4;   counter = 1;
for si = s
    if si > t(j+1)
        j = j + 1;
    end
    alp = zeros(6,1);   d_alp = zeros(6,1);
    aux1 = 1;
    for k = 1:m
        for i = j-m+k:j
            num = si-t(i);
            if num ~= 0
                alp(aux1) = num/(t(i+m+1-k)-t(i));            
            end
            d_alp(aux1) = 1/(t(i+m+1-k)-t(i)); 
            aux1 = aux1 + 1;
        end
    end
    alp_i = 1 - alp;   d_alp_i = -d_alp;
    
    bi = [alp_i(6)*alp_i(4)*alp_i(1);
          alp_i(6)*alp_i(4)*alp(1)+alp_i(6)*alp_i(2)*alp(4)+alp_i(5)*alp_i(2)*alp(6);
          alp_i(6)*alp(4)*alp(2)+alp_i(5)*alp(6)*alp(2)+alp_i(3)*alp(6)*alp(5);
          alp(6)*alp(5)*alp(3)];
    
    d_bi = [d_alp_i(6)*alp_i(4)*alp_i(1)+alp_i(6)*d_alp_i(4)*alp_i(1)+alp_i(6)*alp_i(4)*d_alp_i(1);
            d_alp_i(6)*alp_i(4)*alp(1)+alp_i(6)*d_alp_i(4)*alp(1)+alp_i(6)*alp_i(4)*d_alp(1) + d_alp_i(6)*alp_i(2)*alp(4)+alp_i(6)*d_alp_i(2)*alp(4)+alp_i(6)*alp_i(2)*d_alp(4) + d_alp_i(5)*alp_i(2)*alp(6)+alp_i(5)*d_alp_i(2)*alp(6)+alp_i(5)*alp_i(2)*d_alp(6);
            d_alp_i(6)*alp(4)*alp(2)+alp_i(6)*d_alp(4)*alp(2)+alp_i(6)*alp(4)*d_alp(2) + d_alp_i(5)*alp(6)*alp(2)+alp_i(5)*d_alp(6)*alp(2)+alp_i(5)*alp(6)*d_alp(2) + d_alp_i(3)*alp(6)*alp(5)+alp_i(3)*d_alp(6)*alp(5)+alp_i(3)*alp(6)*d_alp(5);
            d_alp(6)*alp(5)*alp(3)+alp(6)*d_alp(5)*alp(3)+alp(6)*alp(5)*d_alp(3)];
        
    dd_bi = [d_alp_i(6)*d_alp_i(4)*alp_i(1)+d_alp_i(6)*alp_i(4)*d_alp_i(1) + d_alp_i(6)*d_alp_i(4)*alp_i(1)+alp_i(6)*d_alp_i(4)*d_alp_i(1) + d_alp_i(6)*alp_i(4)*d_alp_i(1)+alp_i(6)*d_alp_i(4)*d_alp_i(1);
             d_alp_i(6)*d_alp_i(4)*alp(1)+d_alp_i(6)*alp_i(4)*d_alp(1) + d_alp_i(6)*d_alp_i(4)*alp(1)+alp_i(6)*d_alp_i(4)*d_alp(1) + d_alp_i(6)*alp_i(4)*d_alp(1)+alp_i(6)*d_alp_i(4)*d_alp(1)   +   d_alp_i(6)*d_alp_i(2)*alp(4)+d_alp_i(6)*alp_i(2)*d_alp(4) + d_alp_i(6)*d_alp_i(2)*alp(4)+alp_i(6)*d_alp_i(2)*d_alp(4) + d_alp_i(6)*alp_i(2)*d_alp(4)+alp_i(6)*d_alp_i(2)*d_alp(4)   +   d_alp_i(5)*d_alp_i(2)*alp(6)+d_alp_i(5)*alp_i(2)*d_alp(6) + d_alp_i(5)*d_alp_i(2)*alp(6)+alp_i(5)*d_alp_i(2)*d_alp(6) + d_alp_i(5)*alp_i(2)*d_alp(6)+alp_i(5)*d_alp_i(2)*d_alp(6);
             d_alp_i(6)*d_alp(4)*alp(2)+d_alp_i(6)*alp(4)*d_alp(2) + d_alp_i(6)*d_alp(4)*alp(2)+alp_i(6)*d_alp(4)*d_alp(2) + d_alp_i(6)*alp(4)*d_alp(2)+alp_i(6)*d_alp(4)*d_alp(2)   +   d_alp_i(5)*d_alp(6)*alp(2)+d_alp_i(5)*alp(6)*d_alp(2) + d_alp_i(5)*d_alp(6)*alp(2)+alp_i(5)*d_alp(6)*d_alp(2) + d_alp_i(5)*alp(6)*d_alp(2)+alp_i(5)*d_alp(6)*d_alp(2)   +   d_alp_i(3)*d_alp(6)*alp(5)+d_alp_i(3)*alp(6)*d_alp(5) + d_alp_i(3)*d_alp(6)*alp(5)+alp_i(3)*d_alp(6)*d_alp(5) + d_alp_i(3)*alp(6)*d_alp(5)+alp_i(3)*d_alp(6)*d_alp(5);
             d_alp(6)*d_alp(5)*alp(3)+d_alp(6)*alp(5)*d_alp(3) + d_alp(6)*d_alp(5)*alp(3)+alp(6)*d_alp(5)*d_alp(3) + d_alp(6)*alp(5)*d_alp(3)+alp(6)*d_alp(5)*d_alp(3)];   
    
    b(1,j-3:j,counter) = bi';
    db(1,j-3:j,counter) = d_bi';
    ddb(1,j-3:j,counter) = dd_bi';
    
    aux2 = repmat(bi,1,n);
    aux5 = repmat(d_bi,1,n);
    aux6 = repmat(dd_bi,1,n);
    
    for aux3 = 1:4
        aux4 = 1 + n*(j-4) + n*(aux3-1);
        B(:,aux4:aux4+n-1,counter) = diag(aux2(aux3,:));
        dB(:,aux4:aux4+n-1,counter) = diag(aux5(aux3,:));
        ddB(:,aux4:aux4+n-1,counter) = diag(aux6(aux3,:));
    end    
    
    counter = counter + 1;
end
display('Message:  BSpline built')
