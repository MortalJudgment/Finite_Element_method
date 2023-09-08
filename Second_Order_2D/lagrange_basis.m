function [N,dNdxi]=lagrange_basis(type,coord,dim)
% type: dang ham Lagrange
% coord: toa do diem Gauss
% dim: so chieu
% Ham Lagrange duoc xac dinh tren doan [-1;1] hay tu giac/tam giac tham khao

% returns the lagrange interpolant basis and its gradients w.r.t the
% parent coordinate system.
%
%         [N(xi),dNdxi(xi)]=lagrange_basis(type-order,coord,dim)
%
%   type is the toplogical class of finite element it is in the general
%   form 'topology-#of nodes' ie a three node triangel is T3 a four
%   node quadralateral is Q4 a 4 node tetrahedra is H4 a 27 node brick
%   is B27 etc
%
%   coord is the parent coordinates at which the basis and its
%   gradients are to be evaluated at.
%
%   presently defined are L2, L3, T3, Q4,
%
%   If dim is set to 2 then the vector representation of the N
%   matrix is returned.
%
% written by Jack Chessa
%            j-chessa@northwestern.edu
% Department of Mechanical Engineering
% Northwestern University

if ( nargin == 2 )
    dim=1;
end

switch type
    case 'L2'
        %:::::::::::::::::::: L2 TWO NODE LINE ELEMENT :::::::::::::::::::%
        %
        %    1----------------2
        %
        %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
        
        if size(coord,2) < 1
            disp('Error coordinate needed for the L2 element')
        else
            xi=coord(1);
            N=([1-xi,1+xi]/2)'; % N=[N1,N2]
            dNdxi=[-1;1]/2;
        end
        
    case 'L3'
        %:::::::::::::::::: L3 THREE NODE LINE ELEMENT :::::::::::::::::::%
        %
        %    1--------2--------3
        %
        %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
        
        if size(coord,2) < 1
            disp('Error two coordinates needed for the L3 element')
        else
            xi=coord(1);
            N=[(1-xi)*xi/(-2);(1+xi)*xi/2;1-xi^2]; %[N1 N3 N2]
            dNdxi=[xi-.5;xi+.5;-2*xi];
        end
        
    case 'L4'
        %:::::::::::::::::: L4 FOUR NODE LINE ELEMENT ::::::::::::::::::::%
        %
        %    1------2------3------4
        %
        %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
        if size(coord,2) < 1
            disp('Error two coordinates needed for the L3 element')
        else
            xi=coord(1);
            N=[(-9*xi^3 + 9*xi^2 + xi - 1)/16 ;16*(xi+1)*(xi-1/3)*(xi-1)/27;1-xi^2]; %[N1 N3 N2]
            dNdxi=[(-27*xi^2 + 18*xi +1)/16 ;(16*(xi-1)*(xi-1/3))/27+((16*xi+16)*(xi-1))/27+((16*xi+16)*(xi-1/3))/27;-2*xi];
        end
        
    case 'T3'
        %%%%%%%%%%%%%%%% T3 THREE NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%
        %
        %               3
        %             /  \
        %            /    \
        %           /      \
        %          /        \
        %         /          \
        %        /            \
        %       /              \
        %      /                \
        %     /                  \
        %    1--------------------2
        %
        
        if size(coord,2) < 2
            disp('Error two coordinates needed for the T3 element')
        else
            xi=coord(1); eta=coord(2);
            N=[1-xi-eta;xi;eta];
            dNdxi=[-1,-1;1,0;0,1];
        end
        
        
    case 'Q4'
        %%%%%%%%%%%%%%% Q4 FOUR NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
        %
        %    4--------------------3
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    1--------------------2
        %
        if size(coord,2) < 2
            disp('Error two coordinates needed for the Q4 element')
        else
            xi=coord(1); eta=coord(2);
            N=1/4*[ (1-xi)*(1-eta);
                (1+xi)*(1-eta);
                (1+xi)*(1+eta);
                (1-xi)*(1+eta)];
            dNdxi=1/4*[-(1-eta), -(1-xi);
                1-eta,    -(1+xi);
                1+eta,      1+xi;
                -(1+eta),   1-xi];
        end
    case 'Q9'
        %%%%%%%%%%%%%%% Q9 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
        %
        %    4----------7---------3
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    8          9         6
        %    |                    |
        %    |                    |
        %    |                    |
        %    |                    |
        %    1----------5---------2
        %
        if size(coord,2) < 2
            disp('Error two coordinates needed for the Q9 element')
        else
            xi=coord(1); eta=coord(2);
            N=1/4*[ (xi^2*eta^2 - xi^2*eta - xi*eta^2 +xi*eta);
                    (xi^2*eta^2 - xi^2*eta + xi*eta^2 -xi*eta);
                    (xi^2*eta^2 + xi^2*eta + xi*eta^2 + xi*eta);
                    (xi^2*eta^2 + xi^2*eta - xi*eta^2 - xi*eta);
                    (2*(eta^2 - xi^2*eta^2 + xi^2*eta - eta));
                    (2*(xi^2 - xi^2*eta^2 - xi*eta^2 + xi));
                    (2*(eta^2 - xi^2*eta^2 - xi^2*eta + eta));
                    (2*(xi^2 - xi^2*eta^2 + xi*eta^2 - xi));
                    (4*(- xi^2 - eta^2 + xi^2*eta^2 + 1)) ];
            dNdxi=1/4*[ (eta - 2*eta*xi + 2*eta^2*xi - eta^2), (xi - 2*eta*xi + 2*eta*xi^2 - xi^2);
                    (2*eta*xi - xi + 2*eta*xi^2 - xi^2), (2*eta^2*xi - 2*eta*xi - eta + eta^2);
                    (xi + 2*eta*xi + 2*eta*xi^2 + xi^2), (eta + 2*eta*xi + 2*eta^2*xi + eta^2);
                    (2*eta*xi^2 - 2*eta*xi - xi + xi^2), (2*eta*xi - eta + 2*eta^2*xi - eta^2);
                    (4*eta - 4*eta*xi^2 + 2*xi^2 - 2), (- 4*xi*eta^2 + 4*xi*eta);
                    (- 4*eta*xi^2 - 4*eta*xi), (4*xi - 4*eta^2*xi - 2*eta^2 + 2);
                    (4*eta - 4*eta*xi^2 - 2*xi^2 + 2), (- 4*xi*eta^2 - 4*xi*eta);
                    (- 4*eta*xi^2 + 4*eta*xi), (4*xi - 4*eta^2*xi + 2*eta^2 - 2);
                    (8*eta*xi^2 - 8*eta), (8*xi*eta^2 - 8*xi)  ];
        end
        
    otherwise
        disp(['Element ',type,' not yet supported'])
        N=[]; dNdxi=[];
end

I=eye(dim);
Nv=[];
for i=1:size(N,1)
    Nv=[Nv;I*N(i)];
end

