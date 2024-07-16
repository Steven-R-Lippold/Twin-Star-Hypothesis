AlphS2D4 = [1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 
 3 4; 3 5; 3 6; 3 7; 3 8; 4 5; 4 6; 4 7; 4 8; 5 6; 5 7; 5 8; 6 7; 6 8; 7 8];
% vector that lists all edges in K_8
[rows_S2D4, columns_S2D4] = size(AlphS2D4);

AlphS2ABC = [1 2 3; 1 2 4; 1 2 5; 1 2 6; 1 2 7; 1 2 8; 1 3 4; 1 3 5; 1 3 6; 
1 3 7; 1 3 8; 1 4 5; 1 4 6; 1 4 7; 1 4 8; 1 5 6; 1 5 7; 1 5 8; 1 6 7; 1 6 8; 
1 7 8; 2 3 4; 2 3 5; 2 3 6; 2 3 7; 2 3 8; 2 4 5; 2 4 6; 2 4 7; 2 4 8; 2 5 6; 
2 5 7; 2 5 8; 2 6 7; 2 6 8; 2 7 8; 3 4 5; 3 4 6; 3 4 7; 3 4 8; 3 5 6; 3 5 7; 
3 5 8; 3 6 7; 3 6 8; 3 7 8; 4 5 6; 4 5 7; 4 5 8; 4 6 7; 4 6 8; 4 7 8; 5 6 7; 
5 6 8; 5 7 8; 6 7 8];% all triples (a,b,c) with 1\leq a<b<c\leq 8
[rows_S2abc, columns_S2abc] = size(AlphS2ABC);


AA2CF = readmatrix('CycleFreePartitionsKlein.txt'); % the file 
% CycleFreePartitionsKlein has representatives under Klein group action on
% the set of cycle-free 4-partitions with Gamma_1=X19, 
% (1,3)\in Gamma_2, (1,5)\in Gamma_3 and (1,7)\in Gamma_4

[rows_AA2CF, columns_AA2CF] = size(AA2CF)


vHH=[4,4,1,1,1,1,1,1];% incidence vector for the graph TS4

count=0;
ZZT4=int16.empty(0,28);



for ii=1:rows_AA2CF
    wwtt(1,1:28)=AA2CF(ii,1:28);
if (~(Step3DS(wwtt,vHH,AlphS2D4,AlphS2ABC,rows_S2abc)))
    ZZT4=[ZZT4;wwtt];
    count=count+1;
end
if (ii==10)|(ii==500)|(ii==1000)|(ii==5000)|(ii==10000)|(ii==25000)|(ii==50000)
|(ii==75000)|(ii==100000)|(ii==125000)|(ii==150000)
ii % ii will go up to 154272 
count
end
end

count % to get our claim then count should be zero 
writematrix(ZZT4,'NoPartitions','Delimiter','space');
%  and this file should be empty



%%%%%%%%%%%



function ttt=Step3DS(vxx,vHH,AlphS2D4,AlphS2ABC,rows_S2abc)
scount1=0;
ttt=false; 
while (~(ttt))&(scount1 <rows_S2abc)
    scount1=scount1+1;
    aa1=AlphS2ABC(scount1,1);
    bb1=AlphS2ABC(scount1,2);
    cc1=AlphS2ABC(scount1,3);
    vabc1=ActABC(aa1,bb1,cc1,vxx,AlphS2D4);
 scount2=0;   
while (~(ttt))&(scount2 <rows_S2abc)
    scount2=scount2+1;
    aa2=AlphS2ABC(scount2,1);
    bb2=AlphS2ABC(scount2,2);
    cc2=AlphS2ABC(scount2,3);
    vabc2=ActABC(aa2,bb2,cc2,vabc1,AlphS2D4);
    scount3=0;
  while (~(ttt))&(scount3 <rows_S2abc)
    scount3=scount3+1;
    aa3=AlphS2ABC(scount3,1);
    bb3=AlphS2ABC(scount3,2);
    cc3=AlphS2ABC(scount3,3);
    vabc3=ActABC(aa3,bb3,cc3,vabc2,AlphS2D4);
   ttt=hasDS(vabc3,vHH,AlphS2D4);
  end
end
end
end



function tt=hasDS(txx,vHH,AlphS2D4)
tt=false;
    tx1(1,1:7)=txx(1,1:7);
    tx2(1,1:7)=txx(1,8:14);
    tx3(1,1:7)=txx(1,15:21);
    tx4(1,1:7)=txx(1,22:28);
    mm1=PartToVect(tx1,AlphS2D4);
    cm1=CountVal(mm1);
    mm2=PartToVect(tx2,AlphS2D4);
    cm2=CountVal(mm2);
    mm3=PartToVect(tx3,AlphS2D4);
    cm3=CountVal(mm3);
    mm4=PartToVect(tx4,AlphS2D4);
    cm4=CountVal(mm4);
    if (cm1==vHH)|(cm2==vHH)|(cm3==vHH)|(cm4==vHH)
        tt=true;
    end
end



 
function ccx=CountVal(uu)% associates to a 2x7 vector uu the corespondig 
% vertex multiplicity matrix 1x8 sorted by the largest value
ccx=zeros(1,8);
for i=1:8
    ccx(1,i)=sum(uu(:) == i);
end
ccx=sort(ccx,'descend');
end

 
function xxx=PartToVect(uu,AlphS2D4)
% associates to a 1x7 vector uu the corespondig 2x7 matrix of edges 
xxx=zeros(2,7);
for i=1:7
    xxx(1,i)=AlphS2D4(uu(1,i),1);
    xxx(2,i)=AlphS2D4(uu(1,i),2);
end
end




function swrABC=ActABC(a,b,c,w28,AlphS2D4) 
% gives the action of  the triple [a,b,c]  on the cycle-free partition vector w28 
ABC=zeros(5,28);
ab1=Findy([a,b],AlphS2D4);
ab2=Findy([a,c],AlphS2D4);
ab3=Findy([b,c],AlphS2D4);

for ii=1:28
    if (~(w28(1,ii)==ab1))&(~(w28(1,ii)==ab2))&(~(w28(1,ii)==ab3))
        for jj=1:5
            ABC(jj,ii)=w28(1,ii);
        end
    elseif (w28(1,ii)==ab1)
        ABC(1,ii)=ab3;
        ABC(2,ii)=ab1;
        ABC(3,ii)=ab2;
        ABC(4,ii)=ab2;
        ABC(5,ii)=ab3;
      
    elseif (w28(1,ii)==ab2)
        ABC(1,ii)=ab2;  
        ABC(2,ii)=ab3;
        ABC(3,ii)=ab1;
        ABC(4,ii)=ab3;
        ABC(5,ii)=ab1;
     
    elseif (w28(1,ii)==ab3)
        ABC(1,ii)=ab1;
        ABC(2,ii)=ab2;
        ABC(3,ii)=ab3;
        ABC(4,ii)=ab1;
        ABC(5,ii)=ab2;
    
    end
end 
for kk=1:5
uu1=sort(ABC(kk,1:7));
uu2=sort(ABC(kk,8:14));
uu3=sort(ABC(kk,15:21));
uu4=sort(ABC(kk,22:28));

ww1=PartToVect(uu1,AlphS2D4);
ww2=PartToVect(uu2,AlphS2D4);
ww3=PartToVect(uu3,AlphS2D4);
ww4=PartToVect(uu4,AlphS2D4);
ww11(1,1:7)=ww1(1,1:7);
ww12(1,1:7)=ww1(2,1:7);
ww21(1,1:7)=ww2(1,1:7);
ww22(1,1:7)=ww2(2,1:7);
ww31(1,1:7)=ww3(1,1:7);
ww32(1,1:7)=ww3(2,1:7);
ww41(1,1:7)=ww4(1,1:7);
ww42(1,1:7)=ww4(2,1:7);
G1 = graph(ww11,ww12);
G2 = graph(ww21,ww22);
G3 = graph(ww31,ww32);
G4 = graph(ww41,ww42);
uuu=[uu1,uu2,uu3,uu4];
if (~(hascycles(G1)))&(~(hascycles(G2)))&(~(hascycles(G3)))&(~(hascycles(G4)))
&(~(isequal(uuu,w28)))
sw1=[uu1,uu2,uu3,uu4];
swrABC=rearange(sw1);
end 
end
end





function yyy=Findy(vv,AlphS2D4) %find the row of the pair vv in AlphS2D4
yyy=0;
tty=true; 
while (tty)&(yyy < 28)
    yyy=yyy+1;
    if (vv(1,1)==AlphS2D4(yyy,1))&(vv(1,2)==AlphS2D4(yyy,2))
        tty=false; 
    end 
end 
end 




function re28=rearange(ww)
re28=zeros(1,28);
ww20a=zeros(1,7);
ww20b=zeros(1,7);
ww20c=zeros(1,7);
ww20d=zeros(1,7);
for qqq=1:7
    ww20a(1,qqq)=ww(1,qqq);
    ww20b(1,qqq)=ww(1,qqq+7);
    ww20c(1,qqq)=ww(1,qqq+14);
    ww20d(1,qqq)=ww(1,qqq+21);
end 
re20a=sort(ww20a);
re20b=sort(ww20b);
re20c=sort(ww20c);
re20d=sort(ww20d);

for ppp=1:7
    re28(1,ppp)=re20a(1,ppp);
    re28(1,ppp+7)=re20b(1,ppp);
    re28(1,ppp+14)=re20c(1,ppp);
    re28(1,ppp+21)=re20d(1,ppp);
end 
end
