% The Finite Element Method (FEM)
% for Poisson problem
% *************************************************************************
clc
close all
clear all
format long
global node element;
% *************************************************************************
% ***                            I N P U T                              ***
% *************************************************************************
% tic;
disp('************************************************')
disp('***          S T A R T I N G    R U N        ***')
disp('************************************************')
tic
% disp([num2str(toc),'   START'])

% PLATE PROPERTIES
L=1; D=1; % length and wide of the Omerga domain
%-------------------------------------------------------------------------%
% Number of elements in x-direction.
numx = 4; % mesh 16 (quadrange), mesh 32 (triangle)
% numx = 8 % mesh 64 (quadrange), mesh 128 (triangle)
% numx = 16% mesh 256 (quadrange), mesh 512 (triangle)
% numx = 32% mesh 1024 (quadrange), mesh 2048 (triangle)
% numx = 64
% Number of elements in y-direction.
numy=numx;
%-------------------------------------------------------------------------%

% MESH PROPERTIES
%  elemType = 'T3';
% the element type used in the FEM simulation; 'T3' is for a
% three node constant strain triangular element, 'Q4' is for
% a four node quadrilateral element, and 'Q9' is for a nine
% node quadrila  teral element.

elemType = 'T3';

%plotMesh = 1;    % A flag that if set to 1 plots the initial mesh (to make
%sure u that the mesh is correct)


% *************************************************************************
% ***                 P R E - P R O C E S S I N G                       ***
% *************************************************************************
% GENERATE FINITE ELEMENT MESH
%
% These connectivity matricies refer to the node numbers defined in the
% coordinate matrix node.

disp([num2str(toc),'   GENERATING MESH'])
switch elemType
    case 'Q4'           % here we generate the mesh of Q4 elements
        air=0.0;
        node = singularmesh_cavityflow(L,D,numx,numy,air);
        
        nnx=numx+1;
        nny=numy+1;
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        
        element=make_elem(node_pattern,numx,numy,inc_u,inc_v,1,elemType);
    case 'Q9'           % here we generate the mesh of Q9 elements
        air=0.0;
        node = singularmesh_cavityflow(L,D,numx,numy,air);
        
        nnx=numx+1;
        nny=numy+1;
        inc_u=[1 1 1 1 1 2 1 2 2];
        inc_v=nnx;
        numnode=size(node,1);           % number of nodes
        extrax = numnode + nnx;
        extray = numnode + nnx + 2*numx;
        node_pattern=[ 1 2 nnx+2 nnx+1 numnode+1 extrax+2 extray+1 extrax extrax+1];
        
        element = make_elem(node_pattern,numx,numy,inc_u,inc_v,2,elemType);
        % Update node
        node = [node;zeros(numx*(numy+1)+(2*numx+1)*numy,2)];
        for j=1:size(element,1)
            temp = element(j,:);
            node(temp(5),2) = node(temp(1),2);
            node(temp(5),1) = (node(temp(1),1)+ node(temp(2),1))/2;
            node(temp(6),1) = node(temp(2),1);
            node(temp(6),2) = (node(temp(2),2) + node(temp(3),2))/2;
            node(temp(7),2) = node(temp(3),2);
            node(temp(7),1) = (node(temp(3),1) + node(temp(4),1))/2;
            node(temp(8),1) = node(temp(4),1);
            node(temp(8),2) = (node(temp(4),2) + node(temp(1),2))/2;
            node(temp(9),:) = (node(temp(5),:) + node(temp(7),:))/2;
        end
    case'T3'    % and last but not least T3 elements
        air=0.45;
        node = singularmesh_cavityflow(L,D,numx,numy,air);
        
        nnx=numx+1;
        nny=numy+1;
        inc_u=1;
        inc_v=nnx;

        node_pattern1=[ 1 2 nnx+1 ];
        node_pattern2=[ 2 nnx+2 nnx+1 ];
        
        element=[make_elem(node_pattern1,numx,numy,inc_u,inc_v,1,elemType);
            make_elem(node_pattern2,numx,numy,inc_u,inc_v,1,elemType) ];
    case'T6'
        air=0.45;
        node = singularmesh_cavityflow(L,D,numx,numy,air);
        
        nnx=numx+1;
        nny=numy+1;
        inc_u=[1 1 1 1 2 2];
        inc_v=nnx;
        numnode=size(node,1);           % number of nodes
        extrax = numnode + nnx;
        extray = numnode + nnx + 2*numx;
        
        node_pattern1=[ 1 2 nnx+1 numnode+1 extrax+1 extrax];
        node_pattern2=[ nnx+2 nnx+1 2 extray+1 extrax+1 extrax+2];
        
        element=[make_elem(node_pattern1,numx,numy,inc_u,inc_v,2,elemType);
            make_elem(node_pattern2,numx,numy,inc_u,inc_v,2,elemType) ];
        %Update node
        node = [node;zeros(numx*(numy+1)+(2*numx+1)*numy,2)];
        for j=1:size(element,1)
            temp = element(j,:);
            node(temp(4),2) = node(temp(1),2);
            node(temp(4),1) = (node(temp(1),1)+ node(temp(2),1))/2;
            node(temp(5),:) = (node(temp(2),:) + node(temp(3),:))/2;
            node(temp(6),1) = node(temp(3),1);
            node(temp(6),2) = (node(temp(3),2) + node(temp(1),2))/2;
        end
    otherwise
        disp('Not supported yet!!!')
end
%-------------------------------------------------------------------------%
% Area of each element of the primal mesh
% Dien tich cua moi element
AreaEle = func_area(element,node);
%-------------------------------------------------------------------------%
% Plot meshes
% figure (1)
% plot_mesh(node,element,elemType,'b-');
%-------------------------------------------------------------------------%
numnode=size(node,1);    % number of nodes
numelement=size(element,1); % number of elements
%-------------------------------------------------------------------------%
figure (2)
draw_writedown_numenode(node,element,numelement,elemType)
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   INITIALIZING DATA STRUCTURES'])

nne1U = size(element,2);    % Nummer nodes in 1 element
sdof = numnode;     % So gia tri Unknow
F=zeros(sdof,1);    % external load vector
K=zeros(sdof,sdof); % stiffness matrix

% *************************************************************************
% ***                       P R O C E S S I N G                         ***
% *************************************************************************
%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FECC
disp([num2str(toc),' COMPUTING K MATRIX'])
[W,Q]=quadrature(5, 'GAUSS', 2); % Gauss quadrature order 2 in 2D
valueF = 0;

for e=1:numelement
    sctrU = element(e,:);                        % element scatter vector
    A=zeros(nne1U,nne1U);
    for q=1:size(W,1)                           % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        
        [Nu,dNdxiu]=lagrange_basis(elemType,pt);
        %-----------------------------------------------------------------%
        J0 = node(sctrU,:)'*dNdxiu;             % element Jacobian matrix
        dNdxu=dNdxiu/J0;
        %-----------------------------------------------------------------%
        % COMPUTE B MATRIX 
        B=zeros(2,nne1U);
        B(1,1:nne1U)   = dNdxu(:,1)';
        B(2,1:nne1U)   = dNdxu(:,2)';
        % COMPUTE C MATRIX
        C=zeros(1,nne1U);
        C(1,:)= Nu(:,1)';
        %-----------------------------------------------------------------%
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        A = A + (B'*B + C'*C)*W(q)*abs(det(J0));
    end  % of quadrature loop
    K(sctrU,sctrU)=K(sctrU,sctrU) + A;
    %=====================================================================%
    % Compute for the source term
    for i = 1 : size(sctrU,2)
        x = node(sctrU(i),1);
        y  = node(sctrU(i),2);
        valueF = functionF(x,y);
        Func = AreaEle(e)*valueF/size(sctrU,2);
        F(sctrU(i),1) = F(sctrU(i),1) + Func;
    end
    
end    % of element loop

% spy(K)
%-------------------------------------------------------------------------%
% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])
%Devide Domain into 4 Edge
BelowEdge = [];
AboveEdge = [];
RightEdge = [];
LeftEdge = [];
for i = 1 : numnode
    if (abs(node(i,2)) < 1e-12)
        BelowEdge = [BelowEdge,i];
    end
    if (abs(node(i,1)) < 1e-12)
        LeftEdge = [LeftEdge,i];
    end
    if (abs(node(i,2) - 1) < 1e-12)
        AboveEdge = [AboveEdge,i];
    end
    if (abs(node(i,1) - 1) < 1e-12)
        RightEdge = [RightEdge,i];
    end
end
% Adding boudary condition for each Edge
% Adding at Dirichlet conditon at the end

% bcdof - a vector containging constrained d.o.f
% bcval - a vector containing contained value
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Apply Neumann condition

% Left boundary condition
bcdof = LeftEdge;
nodeLeft = zeros(length(LeftEdge),1);
for i=1:length(LeftEdge)
    nodeLeft(i) = node(bcdof(i),2);     % Consider y-axis
end
bcval = IntergralOnLineF(nodeLeft,bcdof,D/numy,'Left');
g1 = IntergralOnStrLineK(nodeLeft,bcdof,D/numy,sdof);
[K,F]=feaplyc2(K,F,bcdof,bcval,'Neumann',-g1);
% Right boundary condition
bcdof = RightEdge;
nodeRight = zeros(length(RightEdge),1);
for i=1:length(RightEdge)
    nodeRight(i) = node(bcdof(i),2);     % Consider y-axis
end
bcval = IntergralOnLineF(nodeRight,bcdof,D/numy,'Right');
g2 = IntergralOnStrLineK(nodeRight,bcdof,D/numy,sdof);
[K,F]=feaplyc2(K,F,bcdof,bcval,'Neumann',-g2);
% Below boundary condition
bcdof = BelowEdge;
n = length(BelowEdge);
nodeBelow = zeros(n,1);
bcval = zeros(n,1);
for i=1:n
    nodeBelow(i) = node(bcdof(i),1);     % Consider x-axis
end
for i=1:n
    bcval(i) = EgdeCondition(nodeBelow(i),0,'Below');
end
[K,F]=feaplyc2(K,F,bcdof,bcval,'Dirichlet');
% Above boundary condition
bcdof = AboveEdge;
n = length(AboveEdge);
nodeAbove = zeros(n,1);
bcval = zeros(n,1);
for i=1:n
    nodeAbove(i) = node(bcdof(i),1);     % Consider x-axis
end
for i=1:n
    bcval(i) = EgdeCondition(nodeAbove(i),1,'Above');
end
[K,F]=feaplyc2(K,F,bcdof,bcval,'Dirichlet');
%-------------------------------------------------------------------------%
% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\F;
disp([num2str(toc),'   END PROGRAM'])       
for  i = 1 : size(node,1)
    x = node(i,1);
    y = node(i,2);
    Uexact(i,1) = exact_solution(x,y);
end

X=zeros(numx+1,1); Y=zeros(numx+1,1);
for i=1:numx+1
    X(i,1)=(i-1)/numx; Y(i,1)=(i-1)/numx;
end

%[XX,YY]=meshgrid(X,Y);
UU=zeros(numx+1, numx+1); UU_ex=zeros(numx+1, numx+1);
for i=1:(numx+1)*(numx+1)
    ix=mod((i-1),numx+1);
    iy=(i-ix-1)/(numx+1);
    UU(ix+1,iy+1)=U(i,1); UU_ex(ix+1,iy+1)=Uexact(i,1);
end
figure (3)
h1=surf(X,Y,UU);
xlim([0 1]);
ylim([0 1]);
zlim([-2.2 3.2]);
title('Dicrete-solution')

figure (4)
h2=surf(X,Y,UU_ex);
xlim([0 1]);
ylim([0 1]);
zlim([-2.2 3.2]);
title('Exact-solution')