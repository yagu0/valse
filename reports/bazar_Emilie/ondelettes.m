%Dessin ondelettes
%script
figure(1)
courbe = X2(:,1)-mean(X2(:,1));
[C,L] = wavedec(courbe, 4,'haar');

subplot(6,1,1)
plot(1:7,courbe(1:7),'b','LineWidth',2)
axis([0 48 -1.6 1.6])
hold on
plot(8:48,courbe(8:48),'b','LineWidth',2)
%plot(waverec([C(1:12)',zeros(1,36)],L,'haar'),'b--','LineWidth',2)
%courbe2 = X2(:,1);
%[C2,L2] = wavedec(courbe2,4,'haar');
%plot(waverec([zeros(1,3),C2(4:12)',zeros(1,36)],L2,'haar'),'b-.','LineWidth',2)
hold off
ylabel('z','FontSize', 30)
set(gca, 'FontSize', 20)
subplot(6,1,2)
A4=zeros(1,48);
Coeff = zeros(5,48);
for r=1:16
    A4(r)=C(1)/power(2,5/2);
    A4(r+16)=C(2)/power(2,5/2);
    A4(r+32) = C(3)/power(2,5/2);
end
Coeff(5,8)=C(1);
Coeff(5,24) = C(2); 
Coeff(5,40) = C(3);
stairs(A4,'b','LineWidth',2)
axis([0 48 -0.5 0.5])
ylabel('A_4','FontSize', 30)
set(gca, 'FontSize', 20)
subplot(6,1,3)
D4=zeros(1,48);
for r=1:8
    D4(r) = C(4)/power(2,4/2);
    D4(r+8) = -C(4)/power(2,4/2);
    D4(r+16) = C(5)/power(2,4/2);
    D4(r+24) = -C(5)/power(2,4/2);
    D4(r+32) = C(6)/power(2,4/2);
    D4(r+40) = -C(6)/power(2,4/2);
end
Coeff(4,8) = C(4);
Coeff(4,24) = C(5);
Coeff(4,40) = C(6);
stairs(D4,'b','LineWidth',2)
axis([0 48 -0.9 0.9])
ylabel('D_4','FontSize', 30)
        set(gca, 'FontSize', 20)
subplot(6,1,4)
D3=zeros(1,48);
for k=1:12
    for r=1:4
        D3(r+4*(k-1)) = (-1)^(k+1) *C(7+floor((k-1)/2))/power(2,3/2);
    end
end
stairs(D3,'b','LineWidth',2)
ylabel('D_3','FontSize', 30)
axis([0 48 -0.5 0.5])
        set(gca, 'FontSize', 20)
subplot(6,1,5)
D2=zeros(1,48);
for k=1:24
    for r=1:2
        D2(r+2*(k-1)) = (-1)^(k+1) *C(13+floor((k-1)/2))/power(2,2/2);
    end
end
stairs(D2,'b','LineWidth',2)
ylabel('D_2','FontSize', 30)
axis([0 48 -0.8 0.8])
        set(gca, 'FontSize', 20)
subplot(6,1,6)
D1=zeros(1,48);
for k=1:48
    for r=1
        D1(r+1*(k-1)) = (-1)^(k+1) *C(25+floor((k-1)/2))/power(2,1/2);
    end
end
plot(D1,'b','LineWidth',2);
axis([0 48 -0.9 0.9])
ylabel('D_1','FontSize', 30)
set(gca, 'FontSize', 20)
