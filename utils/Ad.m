function [ Adjoint ] = Ad( G )
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
%Computes inv(Ad_{G_i-1_i})
    if numel(G) == 16
        R = G(1:3,1:3)';
        Adjoint = [R,-R*skew(G(1:3,4));zeros(3),R];
    end
end