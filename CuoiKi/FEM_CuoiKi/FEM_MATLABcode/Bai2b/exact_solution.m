function [exactsol] = exact_solution(x,y) 
%---------------------------------
 exactsol = sin(2*pi*x) + cos(2*pi*y) + (x^2 + y^2)/4;
