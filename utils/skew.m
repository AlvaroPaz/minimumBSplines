function [ X ] = skew( x )
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
    if numel(x) == 3
        X = [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
    end
end