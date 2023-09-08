function [AreaEle]= func_area(element,node)
for i = 1 : size(element,1)
    for j = 1 : size(element,2)
        X(j) = node(element(i,j),1);
        Y(j) = node(element(i,j),2);
    end
    AreaEle(i) = polyarea(X,Y);
end
clear X Y