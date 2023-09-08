function element=make_elem(node_pattern,num_u,num_v,inc_u,inc_v,order,elem_type)
% creates a connectivity list
% order of the function
if ( nargin < 6 )
    order = 1;
end

if ( nargin < 5 )
    disp(['Not enough parameters specified for make_elem function'])
end
if order == 1
    inc=zeros(1,size(node_pattern,2));
    e=1;
    element=zeros(num_u*num_v,size(node_pattern,2));
    
    for row=1:num_v
        for col=1:num_u
            element(e,:)=node_pattern+inc;
            inc=inc+inc_u;
            e=e+1;
        end
        inc=row*inc_v;
    end
elseif order == 2
    if ( strcmp(elem_type,'Q9') )
        inc=zeros(1,size(node_pattern,2));
        e=1;
        element=zeros(num_u*num_v,size(node_pattern,2));
        
        for row=1:num_v
            for col=1:num_u
                element(e,:)=node_pattern+inc;
                inc=inc+inc_u;
                e=e+1;
            end
            inc=row*(inc_v*ones(1,size(inc,2)) + [0 0 0 0 2*(inc_v-1) 2*(inc_v-1) 2*(inc_v-1) 2*(inc_v-1) 2*(inc_v-1)]);
        end
    elseif ( strcmp(elem_type,'T6') )
        inc=zeros(1,size(node_pattern,2));
        e=1;
        element=zeros(num_u*num_v,size(node_pattern,2));
        
        for row=1:num_v
            for col=1:num_u
                element(e,:)=node_pattern+inc;
                inc=inc+inc_u;
                e=e+1;
            end
            inc=row*(inc_v*ones(1,size(inc,2)) + [0 0 0 2*(inc_v-1) 2*(inc_v-1) 2*(inc_v-1)]);
        end
    end
else
    disp('Not supported yet')
end
