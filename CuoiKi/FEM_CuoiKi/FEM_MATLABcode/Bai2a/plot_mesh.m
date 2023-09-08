function plot_mesh(X,connect,elem_type,se)
% plots a nodal mesh and an associated connectivity.
% X: is the nodal coordinates, (gia tri diem)
% connect: is the connectivity, (cau noi)
% elem_type: is either 'T3', 'T6', 'Q4', 'Q9',....
% depending on the element topology.
% se: color and line edit
if ( nargin < 4 )
    se='p-';
end

holdState=ishold;       % Check that if "hold on" is on
hold on

% fill X if needed
if (size(X,2) < 3)
    for c=size(X,2)+1:3
        X(:,c)=zeros(size(X,1),1);
    end
end

for e=1:size(connect,1)
    if ( strcmp(elem_type,'Q4') )
        ord=[1,2,3,4,1];
        %%%%%%%%%%%%%%% Q4 FOUR NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%
        %
        %    4---------------3
        %    |               |
        %    |               |
        %    |               |
        %    |               |
        %    |               |
        %    1---------------2
        %
    elseif ( strcmp(elem_type,'Q8') )
        ord=[1,5,2,6,3,7,4,8,1];
        %%%%%%%%%%%%%% Q8 EIGHT NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%
        %
        %    4-------7-------3
        %    |               |
        %    |               |
        %    8               6
        %    |               |
        %    |               |
        %    1-------5-------2
        %
    elseif ( strcmp(elem_type,'Q9') )
        ord=[1,5,2,6,3,7,4,8,1,5,9,7,4,8,9,6];
        %%%%%%%%%%%%%%% Q9 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%
        %
        %    4-------7-------3
        %    |       |       |
        %    |       |       |
        %    8-------9-------6
        %    |       |       |
        %    |       |       |
        %    1-------5-------2
        %
    elseif ( strcmp(elem_type,'Q16') )  % 3-node triangle element
        %       ord=[1,5,6,2,7,8,3,9,10,4,11,12,1];
        %       ord=[1,2,3,4,8,12,16,15,14,13,9,5,1];
        ord=[1,2,3,4,5,6,7,8,9,10,11,12,1];
        %%%%%%%%%%%% Q16 SIXTEEN NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%
        %
        %   10----9----8----7
        %    |              |
        %   11              6
        %    |              |
        %   12              5
        %    |              |
        %    1----2----3----4
        %
    elseif ( strcmp(elem_type,'T3') )  % 3-node triangle element
        ord=[1,2,3,1];
        %%%%%%%%%%%%%%%% T3 THREE NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%% 
    %   
    %           3
    %          / \
    %         /   \
    %        /     \
    %       /       \
    %      /         \
    %     /           \
    %    1-------------2
    %
    elseif ( strcmp(elem_type,'T6') )  % 6-node triangle element
        ord=[1,4,2,5,3,6,1];
        %%%%%%%%%%%%%%%% T6 SIX NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%% 
    %   
    %           3
    %          / \
    %         /   \
    %        /     \
    %       6       5
    %      /         \
    %     /           \
    %    1------4------2
    %
    elseif ( strcmp(elem_type,'H4') )  % 4-node tet element
        ord=[1,2,4,1,3,4,2,3];
    %%%%%%%%%%%%%%%% H4 FOUR NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%%%%%
    %
    %             4
    %           / | \
    %          /  |  \
    %         /   |   \ 
    %        /    |    \ 
    %       /     |     \
    %      1 -----|-----3
    %         \   |   /
    %           \ | /
    %             2  
    elseif ( strcmp(elem_type,'B8') )  % 8-node brick element
        ord=[1,5,6,2,3,7,8,4,1,2,3,4,8,5,6,7];
        %%%%%%%%%%%%%%%%%%% B8 EIGHT NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%%%%
    % 
    %                  8 
    %               /    \    
    %            /          \
    %         /                \
    %      5                     \
    %      |\                     7
    %      |   \                / |
    %      |     \     4    /     |
    %      |        \    /        |
    %      |           6          |
    %      1           |          |
    %       \          |          3
    %          \       |        /
    %            \     |     /
    %               \  |  /
    %                  2
    %                
    end
    for j=1:size(ord,2)
        xpt(j)=X(connect(e,ord(j)),1);
        ypt(j)=X(connect(e,ord(j)),2);
        zpt(j)=X(connect(e,ord(j)),3);
    end
    plot3(xpt,ypt,zpt,se)
end

%rotate3d on
axis auto

if ( ~holdState )
    hold off
end
