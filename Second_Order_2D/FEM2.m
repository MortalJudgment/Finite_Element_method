%-------------------------------------------------------------------------%
% Finite Elements Method in 2 Dimensions-Second order
clc
close all
clear all
format long
%-------------------------------------------------------------------------%
%Intervals
a=0;
b=1;
c=0;
d=1;
n=2;
m=1;
h=(b-a)/n;
k=(d-c)/m;
degree = 2;
%-------------------------------------------------------------------------%
% Calculate nodes and elements
global node element;
r=1;
for j=1:m+1
    for i=1:n+1
        node(r,1) = a+(i-1)*h;
        node(r,2) = c+(j-1)*k;
        r=r+1;
    end
end
r=1;l=1;
for j = 1:m+1
    for i = 1:n
        element(r,1) = l;
        element(r,2) = l+1;
        r=r+1;
        l=l+1;
    end
    l=l+1;
end
l=1;
for i = 1:n+1
    for j = 1:m
        element(r,1) = l;
        element(r,2) = l+n+1;
        r=r+1;
        l=l+1;
    end
end
numnode = size(node,1);        % Number of nodes
numelement = size(element,1);   % Number of elements
switch degree
    case 1
        eletype = 'Q4';
    case 2
        eletype = 'Q9';
%         Update node, element
        for j=1:numelement
            strc = element(j,:);
            i = numnode + j;
            element(j,3) = i;
            node(i,:) = (node(strc(1),:) + node(strc(2),:))/2;
        end
        numnode = size(node,1);
        r = numelement+1;
        s = n+1;
        for i=1:n
            element(r,1) = element(i,3);
            element(r,2) = element(s,3);
            j = numnode + i;
            element(r,3) = j;
            node(j,:) = (node(element(r,1),:) + node(element(r,2),:))/2;
            r = r+1;
            s = s+1;
        end
        numnode = size(node,1);         % Number of nodes
        numelement = size(element,1);   % Number of elements
end
% %-------------------------------------------------------------------------%
% % Plot number node
% for i = 1:numnode
%     x = node(i,1);          
%     y = node(i,2);          
%     text(x(:),y(:),num2str(i(:)));
% end
% %Plot the mesh
% for i = 1:numelement
%     m = element(i,:);       % Consider interval 
%     x = node(m,1);          % 
%     y = node(m,2);          % 
%     patch(x,y,'r');         % Draw vector (x(i),y(i)), (x(i+1),y(i+1))
%     axis([a-1 b+1 c-1 d+1])
% end
% %-------------------------------------------------------------------------%


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
[W,Q]=quadrature(3,'GAUSS' ,2);
% disp('COMPUTING STIFFNESS MATRIX');

for i=1:numelement                  %Consider on each element of the mesh
    sctrU = element(i,:);
    edofU(1:nne1U) = sctrU;
    A = zeros(nne1U,nne1U);
    
    for j=1:size(W,1)   %quadrature loop
        pt = Q(j,:);      %quadrature point
        wt = W(j);      %quadrature weight
        [Nu,dNdxiu] = lagrange_basis(eletype,pt,2);     %element shape functions
        J0 = [4,0;0,2];                %element Jacobian matrix
        dNdxu = dNdxiu * J0^(-1);      %dao ham lagrange
        %----------------------------%
        %---   Compute B matrix   ---%
        %----------------------------%
%         B = zeros(2,nne1U);
        B = dNdxu(sctrU,:)';
        A = A + B'*B*W(j)^2*abs(det(J0));
        x = node(sctrU,1)'*Nu(sctrU);
        y = node(sctrU,2)'*Nu(sctrU);
        f = F1(x,y);
        F(sctrU,1) = F(sctrU,1) + f*Nu(sctrU,:)*W(j)*abs(det(J0));
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
% x = zeros(sdof,1);
% j=1;
% uh = zeros(sdof,1);
% for i=1:numelement
%     x(j) = node(element(i,1));
%     uh(j) = U(element(i,1));
%     x(j+1) = node(element(i,3));
%     uh(j+1) = U(element(i,3));
%     x(j+2) = node(element(i,2));
%     uh(j+2) = U(element(i,2));    
%     j=j+2;
% end
x = zeros(n+1,1);
y = zeros(m+1,1);
uh = zeros(sdof,1);
for i=1:n
    x(i) = node(element(i,1),1);
    x(i+1) = node(element(i,3),1);
    x(i+2) = node(element(i,2),1);
end
for j=1:m
    y(j) = node(element(i,1),2);
    y(j+1) = node(element(i,3),2);
    y(j+2) = node(element(i,2),2);
end
for i=1:numnode
    uex(i) = uexact(x(i),y(i));
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
plot(x,y,U,'x r')
hold on
plot3(x,y,uex,'-. b')
hold off


