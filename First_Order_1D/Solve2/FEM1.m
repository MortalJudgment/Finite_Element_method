%-------------------------------------------------------------------------%
% Finite Elements Method in 1 Dimensions
clc
close all
clear all
format long
%-------------------------------------------------------------------------%
%Intervals
a=0;
b=1;
n=10;
h=(b-a)/n;
%-------------------------------------------------------------------------%
% Calculate nodes and elements 
global node element;
for i=1:n+1
    node(i,1) = a+(i-1)*h;      % Vi tri cua x(i)
end
for i=1:n
    element(i,1)=i;             % Chi so i chi toi x(i)
    element(i,2)=i+1;           % Chi so i+1 chi toi x(i+1)
end
numnode = size(node,1);         % Number of nodes
numelement = size(element,1);   % Number of elements
%-------------------------------------------------------------------------%
% Plot the mesh
for i = 1:numelement
    n = element(i,:);       % Consider interval [i,i+1]
    x = node(n,1);          % Value  x(i),x(i+1)
    y = zeros(size(x));     % Value y(i),y(i+1)
    patch(x,y,'r');         % Draw vector (x(i),y(i)), (x(i+1),y(i+1))
    text(x(:),y(:),num2str(n(:)));      % Add coffiecents 
end
%-------------------------------------------------------------------------%
nne1U = size(element,2);    % So dinh trong mot element
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
[W,Q]=quadrature(2,'GAUSS' ,1);

% disp('COMPUTING STIFFNESS MATRIX');

for i=1:numelement                  %Consider on each element of the mesh
    sctrU = element(i,:);
    edofU(1:nne1U) = sctrU;
    A = zeros(nne1U,nne1U);
    %consider on [-1,1] 
    for j=1:size(W,1)   %quadrature loop
        
        pt = Q(j,:);    %quadrature point
        wt = W(j);      %quadrature weight
        [Nu,dNdxiu] = lagrange_basis('L3',pt);     %element shape functions
        J0 = h/2;         %element Jacobian matrix
        dNdxu = dNdxiu/J0;    %dao ham lagrange tren [-1 1]
        %----------------------------%
        %---   Compute B matrix   ---%
        %----------------------------%
        B = zeros(1,nne1U);
        C = zeros(1,nne1U);
        B(1,1:nne1U) = dNdxu(:,1);
        C(1,1:nne1U) = Nu(:,1);
        A = A + (B'*B*wt + C'*C*wt)*abs(det(J0));
        x = node(sctrU,:)'*Nu;
        f = F1(x);
        F(sctrU,1) = F(sctrU,1) + f*Nu*W(j)*abs(det(J0));
    end                 %end quarature loop
    K(sctrU,sctrU) = K(sctrU,sctrU) + A;
end                     %end element loop
%-------------------------------------------------------------------------%
% APPLY ESSENTIAL BOUNDARY CONDITIONS
% disp('   APPLYING BOUNDARY CONDITIONS   ')
% bcdof: determine the order number of vertices on the boundary Omega
% bcval: determine values of these boundary vertices bcdof
bcdof = [1 sdof];
bcval(1,1) = uexact(node(bcdof(1),:));
bcval(1,2) = uexact(node(bcdof(2),:));
% bcval = [0 0];
% K1 = K;
[K,F] = feaplyc2(K,F,bcdof,bcval);
U = K\F;

% U exactly
uex = zeros(sdof,1);  
for i=1:sdof
    uex(i) = uexact(node(i));
end
figure
plot(node,U,'-.r')
hold on
plot(node,uex,'x b')
axis auto
hold off


