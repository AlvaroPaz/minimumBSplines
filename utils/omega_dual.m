function [ Om ] = omega_dual( V )
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
%Computes ad_dual_{V}
    if numel(V) == 6
        v = skew(V(1:3))';   w = skew(V(4:6))';
        Om = [w,zeros(3);v,w];
    end
end