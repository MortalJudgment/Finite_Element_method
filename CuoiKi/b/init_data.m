function [gcoord,nodes]=init_data(elemType,sizemesh)

switch elemType
    case 'Q4'
        if sizemesh == 1
            fid1=fopen('mesh/node1_quad_size1.txt','r');
            fid2=fopen('mesh/element1_quad_size1.txt','r');
        end
        if sizemesh == 2
            fid1=fopen('mesh/node1_quad_size2.txt','r');
            fid2=fopen('mesh/element1_quad_size2.txt','r');
        end
        if sizemesh == 4
            fid1=fopen('mesh/node1_quad_size4.txt','r');
            fid2=fopen('mesh/element1_quad_size4.txt','r');
        end
        
    case 'T3'
        if sizemesh == 1
            fid1=fopen('mesh/node1_tri_size1.txt','r');
            fid2=fopen('mesh/element1_tri_size1.txt','r');
        end
        if sizemesh == 2
            fid1=fopen('mesh/node1_tri_size2.txt','r');
            fid2=fopen('mesh/element1_tri_size2.txt','r');
        end
        if sizemesh == 4
            fid1=fopen('mesh/node1_tri_size4.txt','r');
            fid2=fopen('mesh/element1_tri_size4.txt','r');
        end
        
    otherwise
        disp('This primal mesh is not considered.')
end
        


% INPUT THE COORDINATE INFORMATION
count=0;
while 1
    tline = fgetl(fid1);
    if isnumeric(tline)%~ischar(tline)
        break
    else
        [C1]= sscanf(tline,'%f %f');
        if size(C1,1)>0
            count=count+1;
            gcoord(count,1:2)=[C1(1) C1(2)];
        end
    end
end
nnode=size(gcoord,1);
fclose(fid1);

count=0;
while 1
    tline = fgetl(fid2);
    if isnumeric(tline)%~ischar(tline)
        break
    else
        if elemType == 'Q4'
            [C2]= sscanf(tline,'%d %f %f %f %f');
            if size(C2,1)>0
                count=count+1;
                nodes(count,1:4)=[C2(1) C2(2) C2(3) C2(4)];
            end
        end
        if elemType == 'T3'
            [C2]= sscanf(tline,'%d %f %f %f');
            if size(C2,1)>0
                count=count+1;
                nodes(count,1:3)=[C2(1) C2(2) C2(3)];
            end
        end

    end
end
fclose(fid2);

    