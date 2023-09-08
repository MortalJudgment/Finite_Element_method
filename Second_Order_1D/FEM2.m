%-------------------------------------------------------------------------%
% Finite Elements Method in 1 Dimensions-Second order
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
degree = 2;
%-------------------------------------------------------------------------%
% Calculate nodes and elements
global node element;
for i=1:n+1
    node(i,1) = a+(i-1)*h;      % Vi tri cua y(i)
end
for i=1:n
    element(i,1) = i;
    element(i,2) = i+1;
end
numnode = size(node,1);         % Number of nodes
numelement = size(element,1);   % Number of elements
switch degree
    case 1
        eletype = 'L2';
    case 2
        eletype = 'L3';
        % Update node, element
        for j=1:numelement
            strc = element(j,:);
            i = numnode + j;
            element(j,3) = i;
            node(i,:) = (node(strc(1),:) + node(strc(2),:))/2;
        end
        numnode = size(node,1);         % Number of nodes
end
%-------------------------------------------------------------------------%
% Plot the mesh
% for i = 1:numelement
%     m = element(i,:);       % Consider interval [i,i+1]
%     x = node(m,1);          % Value  y(2*i-1),y(2*i),y(2*i+1)
%     y = zeros(size(x));     % Cao do
%     patch(x,y,'r');         % Draw vector (x(i),y(i)), (x(i+1),y(i+1))
%     text(x(:),y(:),num2str(m(:)));
% end
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
    edofU(1:nne1U) = sctrU;
    A = zeros(nne1U,nne1U);
    %consider on [-1,1]
    for j=1:size(W,1)   %quadrature loop
        
        pt = Q(j);      %quadrature point
        wt = W(j);      %quadrature weight
        [Nu,dNdxiu] = lagrange_basis(eletype,pt);     %element shape functions
        J0 = h/2;               %element Jacobian matrix
        dNdxu = dNdxiu/J0;      %dao ham lagrange tren [-1 1]
        %----------------------------%
        %---   Compute B matrix   ---%
        %----------------------------%
        B = zeros(1,nne1U);
        B(1,:) = dNdxu(:,1);
        A = A + B'*B*W(j)*abs(det(J0));
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
bcdof = [1 n+1];
bcval(1,1) = uexact(node(bcdof(1),:));
bcval(1,2) = uexact(node(bcdof(2),:));
% bcval = [0 0];
% K1 = K;
[K,F] = feaplyc2(K,F,bcdof,bcval);
U = K\F;

% U exactly
uex = zeros(sdof,1);
x = zeros(sdof,1);
j=1;
uh = zeros(sdof,1);
for i=1:numelement
    x(j) = node(element(i,1));
    uh(j) = U(element(i,1));
    x(j+1) = node(element(i,3));
    uh(j+1) = U(element(i,3));
    x(j+2) = node(element(i,2));
    uh(j+2) = U(element(i,2));    
    j=j+2;
end
for i=1:numnode
    uex(i) = uexact(x(i));
end
% j=1;
% Uca = zeros(sdof,1);
% for i=1:numelement
%     Uca(j) = U(element(i,1));
%     Uca(j+1) = U(element(i,3));
%     Uca(j+2) = U(element(i,2));
%     j=j+2;
% end
figure
plot(x,uh,'- r')
hold on
plot(x,uex,'-. b')
axis([0 1 -1 1])
hold off


