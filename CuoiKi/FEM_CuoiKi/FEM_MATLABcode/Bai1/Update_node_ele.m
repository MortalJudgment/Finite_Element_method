% Update node and element depend on order that we consider
function [newnode,newelement] = Update_node_ele(node,element,ele_type)
if nargin < 3
    ele_type = 'L2';
end
numnode = size(node,1);         % Number of nodes
numelement = size(element,1);   % Number of elements
switch ele_type
    case 'L2'               % 2 node in one element
        % Not change
        newnode = node;
        newelement = element;
    case 'L3'
        i = numnode + 1;
        for k=1:numelement
            strc = element(k,:);
            element(k,3) = i;
            node(i,:) = (node(strc(1),:) + node(strc(2),:))/2;
            i = i+1;
        end
        newnode = node;
        newelement = element;
    case 'L4'
        i = numnode + 1;
        for k=1:numelement
            strc = element(k,:);
            element(k,3) = i;
            j = i+1;
            element(k,4) = j;
            node(i,:) = (2*node(strc(1),:) + node(strc(2),:))/3;
            node(j,:) = (node(strc(1),:) + 2*node(strc(2),:))/3;
            i = i+2;
        end
        newnode = node;
        newelement = element;
    otherwise
        disp('Not supported yet!!')
end
end