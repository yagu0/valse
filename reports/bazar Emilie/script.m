%clear all
load donneesSelec.mat
%donnees=[ind100,res100];
donnees = BB;
Res=mean(donnees,2);
for i=1:349
    X2(:,i)=Res(48*(i-1)+1:48*i);
end
for i=1:348
    signal(:,i)=[X2(:,i);X2(:,i+1)];
end

[p1,n1]=size(X2);
for i=1:n1
    for j=1:p1
        XC(i,j)=X2(j,i)-mean(X2(:,i));
    end
    
    [C,L]=wavedec(XC(i,:)',4,'haar');
    xProj(i,:)=C(1:12);
end

X=xProj(1:end-2,:);
Y=xProj(2:end-1,:);
for i=1:n1
    Xrecon(i,:)=waverec([xProj(i,:),zeros(1,36)],L,'haar');
end
for i=1:348
    signal3(:,i)=[Xrecon(i,:),Xrecon(i+1,:)];
end
