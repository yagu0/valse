%niveau de rouge
cmapRouge = ones(100,3);
for i=1:100
    cmapRouge(i,2) = i/100;
    cmapRouge(i,3) = i/100;
end
%niveau de vert
cmapVert = ones(100,3);
for i=1:100
    cmapVert(i,1) = i/100;
    cmapVert(i,3) = i/100;
end
%niveau de bleu
cmapBleu = ones(100,3);
for i=1:100
    cmapBleu(i,1) = i/100;
    cmapBleu(i,2) = i/100;
end