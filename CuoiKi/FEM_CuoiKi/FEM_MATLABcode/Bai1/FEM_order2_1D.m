%-------------------------------------------------------------------------%
% Finite Elements Method in 1 Dimensions - Order 4
clc
close all
clear all
format long
%-------------------------------------------------------------------------%
%Intervals
a=0;
b=1;
n=5;
h=(b-a)/n;
%-------------------------------------------------------------------------%
% Calculate nodes and elements
global node element;
for i=1:n+1
    node(i,1) = a+(i-1)*h;      % Position x(i)
end
for i=1:n
    element(i,1) = i;               % Index
    element(i,2) = i+1;
end
ele_type = 'L3';    % Consider order on 1D 
                    % number of node in 1 element 
                    % order 1: (2 node) L2, 
                    % order 2: (3 node) L3, 
                    % order 3: (4 node) L4.
%   Update node, element
[node,element] = Update_node_ele(node,element,ele_type);
numnode = size(node,1);         % Number of nodes
numelement = size(element,1);   % Number of elements
%-------------------------------------------------------------------------%
% Plot the mesh
for i = 1:numelement
    m = element(i,:);       % Consider interval [x(i),x(i+1)]
    x = node(m,1);          % Index x
    y = zeros(size(x));     % Height
    patch(x,y,'r');         % Draw vector (x(i),y(i)), (x(i+1),y(i+1))
    text(x(:),y(:),num2str(m(:)));
end
%-------------------------------------------------------------------------%
nne1U = size(element,2);    % Number nodes in 1 element
sdof = numnode;             % So gia tri Unknow
F = zeros(sdof,1);          % The left side of the siffness matrix
K = zeros(sdof,sdof);       % Stiffness matrix

%*************************************************************************%
%***                        P R O C E S S I N G                        ***%
%*************************************************************************%
% W: the weight of Gauss quadrature
% Q: the points of Gauss quadrature
% function quadrature calls: quadrature(the quadrature order, 'type of
% quadrature',the number of spacial dimentions of the problem)
[W,Q]=quadrature(3,'GAUSS' ,1);
% disp('COMPUTING STIFFNESS MATRIX');

for i=1:numelement                  %Consider on each element of the mesh
    sctrU = element(i,:);
    %         edofU(1:nne1U) = sctrU;
    A = zeros(nne1U,nne1U);
    B = zeros(1,nne1U);
    %         consider on [-1,1]
    for j=1:size(W,1)   %quadrature loop
        
        pt = Q(j);      %quadrature point
        wt = W(j);      %quadrature weight
        [Nu,dNdxiu] = lagrange_basis(ele_type,pt);     %element shape functions
        J0 = h/2;               %element Jacobian matrix
        dNdxu = dNdxiu/J0;      %dao ham lagrange tren [-1 1]
        %----------------------------%
        %---   Compute B matrix   ---%
        %----------------------------%
        B = zeros(1,nne1U);
        C = zeros(1,nne1U);
        B(1,:) = dNdxu(:,1);
        C(1,:) = Nu(:,1);
        A = A + (B'*B + C'*C)*wt*abs(det(J0));
        x = node(sctrU,:)'*Nu;
        f = F1(x);
        F(sctrU,1) = F(sctrU,1) + f*Nu*wt*abs(det(J0));
    end                 %end quarature loop
    K(sctrU,sctrU) = K(sctrU,sctrU) + A;
end                     %end element loop
%-------------------------------------------------------------------------%
%   APPLY ESSENTIAL BOUNDARY CONDITIONS
%   disp('   APPLYING BOUNDARY CONDITIONS   ')
%   bcdof: determine the order number of vertices on the boundary Omega
%   bcval: determine values of these boundary vertices bcdof
%-------------------------------------------------------------------------%
bcdof = [1 n+1];        % 2 boundary point
bcval(1) = -( 1 + 2*pi );
bcval(2) = ( 2 + 2*pi );
g = [0 1];
[K,F] = feaplyc2(K,F,bcdof,bcval,'Neumann',g);
u_de = K\F;
%   Fix back u,x dicrete
[x,u] = fixback(node,element,u_de,ele_type);
m = 50;
k = (b-a)/m;
y = a:k:b;
uex = zeros(m+1,1);
for i=1:m+1
    uex(i) = uexact(y(i));
end
figure
plot(x,u,y,uex)
