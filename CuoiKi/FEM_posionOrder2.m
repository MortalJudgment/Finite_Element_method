clc
clear all
close all
format long

%-------------------------------------------------------------------------%
%----------------- POSSION                           ---------------------%
%-------------- - u_xx + u(x) = x + (1+4*pi^2)*sin(2*pi*x)----------------%
%--------------  u'(0)        = 1 + 2*pi                  ----------------%
%--------------  u'(1) + u(1) = 2 + 2*pi                  ----------------%
%-------------------------------------------------------------------------%

% colordef black
% node: coordinates of vertices
% element: the order number of all vertices of each element in the mesh
global node element;
% Mien Omega = [a,b]
a = 0;
b = 1;
N = 8;
% node: the coordinate of nodes
%the number of nodes = N+1
h = (b-a)/N;
degree = 2;

if degree == 1
    order = 'L2';
elseif degree == 2
    order = 'L3';
end

node(N,1) = b;
for i = 1:N+1
    node(i,1) = a + (i-1)*h;
end
% the number of elements of the primal mesh = N
for i = 1:N
    element(i,1) = i;
    element(i,2) = i+1;
end
node2(:,1) = node(:,1);
numelem = size(element,1);  % number of elements

j = 1;
node2(1,1) = node(1,1);
for iel=1:numelem
    n = element(iel,:);
    i = numelem + iel + 1;
    element(iel,3) = i;
    node(i,:) = (node(n(1),:)+node(n(2),:))/2;
    node2(j+1,:) = node(i,:);
    node2(j+2,:) = node(n(2),:);
    j = j+2;
end

numnode = size(node,1);     % number of nodes

% plot the mesh
% plot_mesh(node,element,elemType, 'b-');
for iel=1:numelem
    n = element(iel,:);
    x = node(n,1);
    %--------------------
    % 1D
    y = zeros(size(x));
    %--------------------
    % 2D
    % y = node(n,2);
    %--------------------
    patch(x,y,'w');
    text(x(:),y(:),num2str(n(:)));
end

%--------------------------------------------------------------------------
nnelU = size(element,2);  % so dinh trong mot element = 3
sdof = numnode; % the number of unknows
F = zeros(sdof,1); % the left hand side of the stiffness matrix
K = zeros(sdof,sdof); % stiffness matrix

%**************************************************************************
%***                      P R O C E S S I N G                          ***
%**************************************************************************
% W: the weight of Gauss quadrature
% Q: the points of Gauss quadrature on [-1,1]
% [W,Q] = quadrature(the number of Gauss points, 'GAUSS', 1D);
[W,Q] = quadrature(3, 'GAUSS', 1);
% W: trong so cau phuong
% Q: diem cau phuong xet tren doan [-1,1]
disp('COMPUTING STIFFNESS MATRIX')
for e = 1:numelem              % consider on each element of the mesh
    sctrU = element(e,:);      % element scatter vector
    edofU(1:nnelU) = sctrU;
    A = zeros(nnelU,nnelU);
    
    % consider on [-1,1]
    for q = 1:size(W,1)     %quadrature loop
        pt = Q(q,:);        %quadrature point
        wt = W(q);          %quadrature weight
        [Nu,dNdxiu] = lagrange_basis(order,pt);  %element shape functions
        
        J0 = h/2;           %element Jacobian matrix
        dNdxu = dNdxiu/J0;
        % dNdxu: dao ham ham Lagrange tren doan [-1,1]
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% COMPUTE B MATRIX
        B = zeros(1,nnelU);
        B(1,1:nnelU) = dNdxu(:,1);
        
        % COMPUTE C MATRIX
        C = zeros(1,nnelU);
        C(1,1:nnelU)   = Nu(:,1)';
        %----------------------------------------------------
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        A = A + (B'*B + C'*C)*W(q)*abs(det(J0));
        x = node(sctrU,:)'*Nu;
        f = source_term(x);
        F(sctrU,1) = F(sctrU,1) + f*Nu*W(q)*abs(det(J0));
    end   % of quadrature loop
    K(sctrU,sctrU) = K(sctrU,sctrU) + A;
end
% of element loop
%-----------------------------------------------------------------
% % APPLY ESSENTIAL BOUNDARY CONDITIONS
% disp(' APPLYING BOUNDARY CONDITIONS')
% % ENFORCE BOUNDRY CONDITION
% % bcdof: determine the order number of vertices on the boundary Omega
% % bcval: determine valuse of these boundary vertices bcdof
bcdof = [1 N+1];
bcval(1,1)= exact_sol(node(bcdof(1),:));
bcval(1,2)= exact_sol(node(bcdof(2),:));
%bcval = [0 0];
[K,F] = feaplyc2(K,F,bcdof,bcval);
% U = conjgrad(K,F,10e-6);
U = K\F;

% Edit node
x = zeros(sdof,1);
j = 1;
uh = zeros(sdof,1);
for i = 1:numelem
    x(j) = node(element(i,1));
    uh(j) = U(element(i,1));
    x(j+1) = node(element(i,3));
    uh(j+1) = U(element(i,3));
    x(j+2) = node(element(i,2));
    uh(j+2) = U(element(i,2));
    j = j+2;
end
% U exact
U_ex = zeros(sdof,1);
for i = 1:sdof
    U_ex(i) = exact_sol(x(i));
end

figure
plot(x, uh, '*r',x, U_ex ,'b');
% figure
% plot(node, U_ex ,'b');