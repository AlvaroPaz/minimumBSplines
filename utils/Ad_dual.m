function [ Adjoint ] = Ad_dual( G )
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
%Computes inv(Ad_dual{G_i_i+1})
    if numel(G) == 16
        R = G(1:3,1:3);
        Adjoint = [R,zeros(3);skew(G(1:3,4))*R,R];
    end
end