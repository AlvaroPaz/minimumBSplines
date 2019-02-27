function [ MT ] = DH( theta,d,a,alpha )
%%%--------------------------
% Alvaro Paz
% Cinvestav - Saltillo Campus
%%%--------------------------
% Rot(z)Tras(z)Tras(x)Rot(x)
MT=[cos(theta),     -sin(theta)*cos(alpha), sin(theta)*sin(alpha),  a*cos(theta);
    sin(theta),     cos(theta)*cos(alpha),  -cos(theta)*sin(alpha), a*sin(theta);
    0,              sin(alpha),             cos(alpha),             d;
    0,              0,                      0,                      1];
end