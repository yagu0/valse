%clear all
load donneesSelec.mat
donnees=BB;

Res=mean(donnees,2);
%Ind=mean(donnees(:,1:100),2);
for i=1:340
    X2(:,i)=Res(48*(i-1)+1:48*i);
end
for i=1:339
    signal(:,i)=[X2(:,i);X2(:,i+1)];
end

[p1,n1]=size(X2);
for i=1:n1
    [C,L]=wavedec(X2(:,i)',4,'haar');
    xProj(i,:)=C(4:12);
end

X=xProj(1:end-2,:);
Y=xProj(2:end-1,:);
for i=1:n1
    Xrecon(i,:)=waverec([zeros(1,3),xProj(i,:),zeros(1,36)],L,'haar');
end

for i=1:339
    signal2(:,i)=[Xrecon(i,:),Xrecon(i+1,:)];
end
