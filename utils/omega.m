function [ Om ] = omega( V )
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
%Computes ad_{V}
    if numel(V) == 6
        v = skew(V(1:3));   w = skew(V(4:6));
        Om = [w,v;zeros(3),w];
    end
end