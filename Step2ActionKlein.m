AlphS2D4 = [1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 
 3 4; 3 5; 3 6; 3 7; 3 8; 4 5; 4 6; 4 7; 4 8; 5 6; 5 7; 5 8; 6 7; 6 8; 7 8];
[rows_S2D4, columns_S2D4] = size(AlphS2D4);



AA2CF = readmatrix('CycleFreePartitionsGamma1X19.txt');
[rows_AA2CF, columns_AA2CF] = size(AA2CF)




sigma8=[1,2,3,4,5,8,7,6];
sigma35=[1,2,5,4,3,6,7,8];

count=0;
ZZT4=int16.empty(0,28);
YY=zeros(4,28);
XX=AA2CF;
[rows_xx, columns_xx] = size(XX); 
while rows_xx>0 
   wwt(1,1:28)=XX(1,1:28);
 wws1=Perm23(sxw(sigma35,wwt,AlphS2D4));
 wws2=sxw(sigma8,wws1,AlphS2D4);
 wws3=sxw(sigma8,wwt,AlphS2D4);
YY(1,1:28)=wwt;
YY(2,1:28)=wws1;
YY(3,1:28)=wws2;
YY(4,1:28)=wws3;
  
   ZZT4=[ZZT4;wwt];
   count=count+1;
XX=setdiff(XX,YY,'rows','legacy');
 if (count==100)|(count==1000)|(count==10000)(count==25000)|(count==50000)
|(count==75000)|(count==100000)|(count==125000)|(count==150000)
     count %the count will go up to 154272 
     rows_xx
 end
 [rows_xx, columns_xx] = size(XX);
end

count
writematrix(ZZT4,'CycleFreePartitionsKlein','Delimiter','space'); 






function per23=Perm23(w28) %gives the action of (2,3) on the partition vector
% w28=(G1,G2,G3,G4) (from V^{28}). 
pp1(1,1:7)=w28(1,1:7);
pp2(1,1:7)=w28(1,8:14);
pp3(1,1:7)=w28(1,15:21);
pp4(1,1:7)=w28(1,22:28);

per23=[pp1,pp3,pp2,pp4];
end



function swr=sxw(sigma8,w28,AlphS2D4) %gives the action of sigma8 (from S_8) on 
%the vector w28 (from V^{28}). 
sw=zeros(1,28);
for ii=1:28
    su1=sigma8(1,AlphS2D4(w28(1,ii),1));
    su2=sigma8(1,AlphS2D4(w28(1,ii),2));
    vv=sort([su1,su2]);
    sw(1,ii)=Findy(vv,AlphS2D4);
end 
swr=rearange(sw);
end





function yyy=Findy(vv,AlphS2D4) %finds the row of the pair vv in AlphS2D4
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
