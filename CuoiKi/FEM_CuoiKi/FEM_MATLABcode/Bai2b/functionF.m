function [valueF] = functionF(x,y)
% the exact solution =  sin(2*pi*x) + cos(2*pi*y) + (x^2 + y^2)/4;
valueF = (1 + 4*pi^2)*(sin(2*pi*x) + cos(2*pi*y)) + (x^2 + y^2)/4 - 1;
end
