function gcoord = singularmesh_cavityflow(lengthx, lengthy, lx,ly, air)

% lengthx:           length of x-axis side of problem
% lengthy:           length of y-axis side of problem
% lx     :           number of element in x-axis
% ly     :           number of element in y-axis
% air    :           singularity scale coefficient   
%--------------------------------------------
%input data for nodal coordinate values
%gcoord(i,j) where i->node no.  and j->x or y
%--------------------------------------------
% nel=lx*ly;
dx=lengthx/lx;          % the length of side of element in x-axis
dy=lengthy/ly;          % the length of side of element in y-axis
gcoord=[];
for i=1:ly+1
    for j=1:lx+1
        gcoord=[gcoord; (j-1)*dx (i-1)*dy;]; 
    end
end
    
nn=0;
for ip=1:lx+1                % sampling node coordiantes for discretisation
    for iq=1:ly+1 
        nn=nn+1;     

%        r=random('beta',1,1);        % r=[0 1];
%        r=air*(2*r-1);               % project r=[-air air];
         r = 0;

        if ip==1 || iq==ly/2+1|| ip==lx+1 || iq==1 || iq==ly+1 %| (iq-1)*bb==5 %|(ip-1)*aa==L/2
             r=0;
        end
        gcoord(nn,1)=gcoord(nn,1)+dx*r;
        gcoord(nn,2)=gcoord(nn,2)+dy*r;  
    end
end
%-------------------------------------------------------------------------%
%---------- Draw mesh and write down the number of node ------------------%
%-------------------------------------------------------------------------%
% clear nn;
% nodes=[];
% for i=1:lx
%     for j=1:ly
%         nodes=[nodes; (ly+1)*(i-1)+j (ly+1)*i+j (ly+1)*i+j+1 (ly+1)*(i-1)+j+1;]
%     end
% end
% 
% figure;
% h=gcf;
% %set(h,'name','plate problem 2D');
% %set(h,'NumberTitle','off');
% axis equal;
% %title('Mesh configuration');
% for iel=1:nel
%    n1=nodes(iel,1);
%    n2=nodes(iel,2);
%    n3=nodes(iel,3);
%    n4=nodes(iel,4);
%    x(1)=gcoord(n1,1);
%    x(2)=gcoord(n2,1);    
%    x(3)=gcoord(n3,1);
%    x(4)=gcoord(n4,1);
%    y(1)=gcoord(n1,2);    
%    y(2)=gcoord(n2,2);    
%    y(3)=gcoord(n3,2);
%    y(4)=gcoord(n4,2);
%    
%    patch(x(1:4),y(1:4),'w');
%    text(x(1),y(1),num2str(n1));
%    text(x(2),y(2),num2str(n2));
%    text(x(3),y(3),num2str(n3));
%    text(x(4),y(4),num2str(n4));
%    text(xx,yy,num2str(iel)); 
% end