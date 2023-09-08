% Approximatly intergral on line of u\phi dS
function G = IntergralOnStrLineK(node,ele,h,numbernode)
% Fixback ele positon
temp1 = 0; temp2 = 0;
for i = 1 : length(ele)-1
    for j = i+1:length(ele)
        if node(i) > node(j)
            %Swap value
            temp1 = node(i);
            node(i) = node(j);
            node(j) = temp1;
            temp2 = ele(i);
            ele(i) = ele(j);
            ele(j) = temp2;
        end
    end
end
%---------------------------------------%
% Build an element
element = []; k=1;
for i = 1 : length(ele)-1
    element = [element; k k+1];
    k = k+1;
end
numnode = size(node,1);
numelement = size(element,1);
nne1U = size(element,2);
sdof = numnode;
K = zeros(sdof,sdof);
[W,Q]=quadrature(3,'GAUSS' ,1);
for i=1:numelement                  %Consider on each element of the mesh
    sctrU = element(i,:);
    A = zeros(nne1U,nne1U);
    for j=1:size(W,1)   %quadrature loop
        pt = Q(j);      %quadrature point
        wt = W(j);      %quadrature weight
        [Nu,dNdxiu] = lagrange_basis('L2',pt);     %element shape functions
        J0 = h/2;               %element Jacobian matrix
        dNdxu = dNdxiu/J0;      %dao ham lagrange tren [-1 1]
        %----------------------------%
        %---   Compute C matrix   ---%
        %----------------------------%
        C = zeros(1,nne1U);
        C(1,:) = Nu(:,1);
        A = A + (C'*C)*wt*abs(det(J0));
    end                 %end quarature loop
    K(sctrU,sctrU) = K(sctrU,sctrU) + A;
end
% Tranform depend on position
G = zeros(numbernode,numbernode);
%Set up follow row
for i=1:length(ele)
    for j=1:length(ele)
        m = ele(i);
        n = ele(j);
        G(m,n) = K(i,j);
    end
end    
    
    
    
    