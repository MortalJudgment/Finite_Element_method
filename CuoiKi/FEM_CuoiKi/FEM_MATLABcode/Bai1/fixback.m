% Fix back x and u dicrete
function [x,u] = fixback(node,element,u_before,ele_type)
if nargin == 3
    ele_type = 'L2'; % Default value
end
numnode = size(node,1);
numele = size(element,1);
x = zeros(numnode,1);
u = zeros(numnode,1);
switch ele_type
    case 'L2'
        x = node;
        u = u_before;
    case 'L3'
        j=1;
        for i = 1:numele
            x(j) = node(element(i,1));
            u(j) = u_before(element(i,1));
            x(j+1) = node(element(i,3));
            u(j+1) = u_before(element(i,3));
            x(j+2) = node(element(i,2));
            u(j+2) = u_before(element(i,2));
            j = j+2;
        end
    case 'L4'
        j=1;
        for i = 1:numele
            x(j) = node(element(i,1));
            u(j) = u_before(element(i,1));
            x(j+1) = node(element(i,3));
            u(j+1) = u_before(element(i,3));
            x(j+2) = node(element(i,4));
            u(j+2) = u_before(element(i,4));
            x(j+3) = node(element(i,2));
            u(j+3) = u_before(element(i,2));
            j = j+3;
        end
end
    