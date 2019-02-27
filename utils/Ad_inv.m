function [ Adjoint ] = Ad_inv( G )
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
% Computes Ad_{G_i-1_i}
    if numel(G) == 16
        R = G(1:3,1:3);
        Adjoint = [R,skew(G(1:3,4))*R;zeros(3),R];
    end
end