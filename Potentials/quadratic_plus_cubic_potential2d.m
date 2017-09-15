%% Computation of a quadratic plus quartic potential 
%% INPUTS:
%%          mass, omega_x,omega_y,alpha_x,alpha_y: Potential parameters (double)
%%          X,Y: Computational grids (matrix)
%% OUTPUT:
%%          Potential : values of the potential over the computational grid (matrix)

function [Potential] = quadratic_plus_cubic_potential2d(mass, omega_x, omega_y,alpha_x,alpha_y,X,Y)
Potential = mass*((omega_x*X).^2 + (omega_y*Y).^2)/2 + (alpha_x*X.^3 + alpha_y*Y.^3); 
% Computing the quadratic plus cubic potential