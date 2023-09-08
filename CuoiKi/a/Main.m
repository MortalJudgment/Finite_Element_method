clc
close all
clear all
format long
disp('************************************************')
disp('***          S T A R T I N G    R U N        ***')
disp('************************************************')
tic
L=1; D=1; % length and wide of the Omerga domain
numx = 4;               %Number of elements in x-direction.
numy = numx;            %Number of elements in y-direction.
elemType = 'T3';
[X,Y,U] = FEM_ForPoissonEq(L,D,numx,numy,elemType);
k=1;
for  i = 1 : numy
    for j = i : numx
        Uexact(k) = exact_solution(X(j),Y(i));
        k=k+1;
    end
end
UU=zeros(numx+1, numx+1); UU_ex=zeros(numx+1, numx+1);
for i=1:(numx+1)*(numx+1)
    ix=mod((i-1),numx+1);
    iy=(i-ix-1)/(numx+1);
    UU(ix+1,iy+1)=U(i,1); UU_ex(ix+1,iy+1)=Uexact(i);
end
figure
h1=surf(X,Y,UU);
title('Dicrete-solution')
%xlim([0 1]);
%ylim([0 1]);
%zlim([0.6 2.2]);
figure 
h2=surf(X,Y,UU_ex);
title('Exact-solution')
%xlim([0 1]);
%ylim([0 1]);
%zlim([0.6 2.2]);
% disp([num2str(toc),'   END PROGRAM'])
