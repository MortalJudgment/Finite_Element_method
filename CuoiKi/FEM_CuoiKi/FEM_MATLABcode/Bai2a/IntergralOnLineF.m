% Approximatly intergral on line of h\phi dS
function F = IntergralOnLineF(node,ele,h,position)
% Fixback ele positon
temp = 0;
for i = 1 : length(ele)-1
    for j = i+1:length(ele)
        if node(i) > node(j)
            %Swap value
            temp = node(i);
            node(i) = node(j);
            node(j) = temp;
        end
    end
end
%---------------------------------------%
% Build an element
element = [];
k = 1;
for i = 1 : length(ele)-1
    element = [element; k k+1];
    k = k+1;
end
numnode = size(node,1);
numelement = size(element,1);
nne1U = size(element,2);
sdof = numnode;
F = zeros(sdof,1);
[W,Q]=quadrature(3,'GAUSS' ,1);
for i=1:numelement                  %Consider on each element of the mesh
    sctrU = element(i,:);
    for j=1:size(W,1)   %quadrature loop
        pt = Q(j);      %quadrature point
        wt = W(j);      %quadrature weight
        [Nu,dNdxiu] = lagrange_basis('L2',pt);     %element shape functions
        J0 = h/2;               %element Jacobian matrix
        dNdxu = dNdxiu/J0;      %dao ham lagrange tren [-1 1]
        x = node(sctrU,:)'*Nu;
        f = EgdeCondition(0,x,position);
        F(sctrU,1) = F(sctrU,1) + f*Nu*wt*abs(det(J0));
    end                 %end quarature loop
end
% Return
for i = 1 : length(ele)-1
    for j = i+1:length(ele)
        if ele(i) > ele(j)
            %Swap value
            temp = F(i);
            F(i) = node(j);
            F(j) = temp;
        end
    end
end  
end
