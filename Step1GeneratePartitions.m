%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Generate the relevant partitions needed. 
%          Recall that not all of the partitions are needed
%          (Only need the last graph in the partition weakly 
%          equivalent to T19).
%          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vector that lists all edges in K_8
% this list enumerates the edges using the rows of the matrix
AlphS2D4 = [1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 
 3 4; 3 5; 3 6; 3 7; 3 8; 4 5; 4 6; 4 7; 4 8; 5 6; 5 7; 5 8; 6 7; 6 8; 7 8]; 
[rows_S2D4, columns_S2D4] = size(AlphS2D4);


% Up to a group action, uu1 corresponds to T19
uu1=[1,3,5,7,8,10,21]; % partition corresponding to Gamma_1
uu2=[2,0,0,0,0,0,0];   % Gamma_2
uu3=[4,0,0,0,0,0,0];   % Gamma_3
uu4=[6,0,0,0,0,0,0];   % Gamma_4

% remaining edges to assign
v=[9,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28];

AA2=nchoosek(v,6); % last 6 entries of Gamma_2 (including those with cycles) 
[rows_AA2, columns_AA2] = size(AA2)

% initialize starting array
ZZT2=int16.empty(0,7);
count=0;

% start to fill in the remaining possible positions
for i=1:rows_AA2
% start by filling in Gamma_2
uu2(1,1)=2;
uu2(1,2:7)=AA2(i,1:6);
% use PartToVect function to get edges
ww=PartToVect(uu2,AlphS2D4);
% split edges
ww1(1,1:7)=ww(1,1:7);
ww2(1,1:7)=ww(2,1:7);
% make into a graph and search for cycles
% exclude if there are cycles
G = graph(ww1,ww2);
if ~(hascycles(G))
ZZT2=[ZZT2;uu2];
count=count+1;
end
% counter to be sure nothing is delaying the program
if (i==1000)|(i==2000)|(i==4000)|(i==5000)|(i==10000)
i 
count
end
end

% put in the result if it worked out to not have a cycle
AA2CF=ZZT2; % last 6 entries of Gamma_2 (only those that are cycle-free) 

[rows_AA2CF, columns_AA2CF] = size(AA2CF)

% start to continue the process with the other graphs
% similar process, but makes the possible edges left smaller
ZZT4=int16.empty(0,28);
count=0;
for i=1:rows_AA2CF
xx2(1,1:6)=AA2CF(i,2:7);
uu2(1,1:7)=AA2CF(i,1:7);
rest=setdiff(v,xx2);
CC3=nchoosek(rest,6);
for jj=1:924
 xx3=CC3(jj,1:6);
uu3(1,1)=4;
uu3(1,2:7)=xx3(1,1:6);
% whatever is not in Gamma_3 had to be in Gamma_4
xx4=setdiff(rest,xx3);
uu4(1,1)=6;
uu4(1,2:7)=xx4(1,1:6);
ww3=PartToVect(uu3,AlphS2D4);
ww31(1,1:7)=ww3(1,1:7);
ww32(1,1:7)=ww3(2,1:7);
G3 = graph(ww31,ww32);
ww4=PartToVect(uu4,AlphS2D4);
ww41(1,1:7)=ww4(1,1:7);
ww42(1,1:7)=ww4(2,1:7);
G4 = graph(ww41,ww42);
if (~(hascycles(G3)))&(~(hascycles(G4)))
    wwwf=[uu1,uu2,uu3,uu4];
ZZT4=[ZZT4;wwwf];
count=count+1;
end
end
if (i==100)|(i==1000)|(i==2000)|(i==3000)|(i==4000)|(i==5000)
i 
count
% this will be 5880
end
end


count
writematrix(ZZT4,'CycleFreePartitionsGamma1X19','Delimiter','space');
% the file CycleFreePartitionsGamma1X19 has a list of all cycle-free partitions 
% with Gamma_1=X19, (1,3)\in Gamma_2, (1,5)\in Gamma_3 and (1,7)\in Gamma_4



%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Converts partition into n x 2 vector 
%          Each row represents an edge
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xxx=PartToVect(uu,AlphS2D4)
% associates to a 7x1 vector uu the corespondig matrix of edges 
xxx=zeros(2,7);
for i=1:7
    xxx(1,i)=AlphS2D4(uu(i),1);
    xxx(2,i)=AlphS2D4(uu(i),2);
end
end
