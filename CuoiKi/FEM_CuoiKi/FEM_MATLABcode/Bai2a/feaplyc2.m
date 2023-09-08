 function [kk,ff]=feaplyc2(kk,ff,bcdof,bcval,Type,g)

%-------------------------------------------------------------------------%
%-- Purpose:                                                            --%
%--   Apply constraints to matrix equation [kk]{x}={ff}                 --%
%--                                                                     --%
%-- Synopsis:                                                           --%
%--    [kk,ff]=feaplybc(kk,ff,bcdof,bcval)                              --%
%--                                                                     --%
%-- Variable Description:                                               --%
%--   kk:  system matrix before applying constraints                    --%
%--   ff:  system vector before applying constraints                    --%
%--   bcdof:  a vector containging constrained d.o.f                    --%
%--   bcval:  a vector containing contained value                       --%
%--   Type :  type of boundary condition (Dirichlet,Neumann,...)        --%
%-------------------------------------------------------------------------%
if nargin < 5
    Type = 'Dirichlet';
    g = 0;
end
n = length(bcdof);
sdof = size(kk);
switch Type
    case 'Dirichlet'
        for i=1:n
            c = bcdof(i);
            for j=1:sdof
                kk(c,j) = 0;
            end
            kk(c,c) = 1;
            ff(c) = bcval(i);
        end
    case 'Neumann'
        for i=1:n
            c = bcdof(i);
            ff(c) = ff(c) + bcval(i);
        end
        kk = kk + g;
end