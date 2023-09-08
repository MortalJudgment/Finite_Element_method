function draw_writedown_numenode(nodes,elements,numelement,type)
% node   :           position of every single point
% elements:          element
% numelement:        number of element
if( strcmp(type,'Q9') == 1 )
    n = zeros(9,1);
    x = zeros(9,1);
    y = zeros(9,1);
    for iel=1:numelement
        for i=1:9
            n(i) = elements(iel,i);
            x(i) = nodes(n(i),1);
            y(i) = nodes(n(i),2);
        end
        patch(x(1:4),y(1:4),'w');
        patch([x(5),x(7)],[y(5),y(7)],'w');
        patch([x(6),x(8)],[y(6),y(8)],'w');
        for i=1:9
            text(x(i),y(i),num2str(n(i)));
        end
    end
elseif (strcmp(type,'T6') == 1 )
    n = zeros(6,1);
    x = zeros(6,1);
    y = zeros(6,1);
    for iel=1:numelement
        for i=1:6
            n(i) = elements(iel,i);
            x(i) = nodes(n(i),1);
            y(i) = nodes(n(i),2);
        end
        patch(x(1:3),y(1:3),'w');
        for i=1:6
            text(x(i),y(i),num2str(n(i)));
        end
    end
elseif (strcmp(type,'T3') == 1 )||(strcmp(type,'Q4') == 1 )
    
    for iel=1:numelement
        n=elements(iel,:);
        x=nodes(n,1);
        y=nodes(n,2);
        patch(x,y,'w');
        text(x(:),y(:),num2str(n(:)));
    end
else
    disp('Not supported yet')
end

