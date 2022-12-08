function [Lambda, Mu] = Lame(E0, Nu0)
% LAME.m compute lame coefficients from E (modulus) and Nu (poisson ratio).
% Usage:
%   E  = 1.0;  (MPa)
%   Nu = 0.35; (-)
%   [lam,mu] = Lame(E,Nu);

    Lambda = (Nu0 * E0) / ((1 + Nu0) * (1 - 2 * Nu0));
    Mu = E0 / (2 * (1 + Nu0));
end
