for i = 1:100
        courbeX(i,:) = data(((i-1)*96+1):((i-1)*96+48));
        courbeY(i,:) = data(((i-1)*96+49):((i-1)*96+96));
end
plot((courbeX)')
figure(2)
plot((courbeY)')

X2 = (courbeX)'