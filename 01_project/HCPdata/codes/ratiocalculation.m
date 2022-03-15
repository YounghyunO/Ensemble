gradient1 = zeros(size(EC,1),1);
for i = 1:size(EC,1)
    pos_of = sum(max(EC(:,i), 0));
    pos_if = sum(max(EC(i,:),0));
    gradient1(i,1) = pos_of/(pos_if+pos_of);
end

gradient2 = zeros(size(EC,1),1);
for i = 1:size(EC,1)
    neg_of = -sum(min(EC(:,i), 0));
    neg_if = -sum(min(EC(i,:),0));
    gradient2(i,1) = neg_of/(neg_if+neg_of);
end

gradient3 = zeros(size(EC,1),1);
for i = 1:size(EC,1)
    pos_of = sum(max(EC(:,i), 0));
    pos_if = sum(max(EC(i,:),0));
    gradient3(i,1) = pos_if/(pos_if+pos_of);
end

gradient4 = zeros(size(EC,1),1);
for i = 1:size(EC,1)
    neg_of = -sum(min(EC(:,i), 0));
    neg_if = -sum(min(EC(i,:),0));
    gradient4(i,1) = neg_if/(neg_if+neg_of);
end

posvec = max(EC(:),0);
pos1p = prctile(posvec,90);
for i = 1:length(posvec)
    if posvec(i,1) < pos1p
        posvec(i,1) = 0;
    end
end

negvec = min(EC(:),0);
neg1p = prctile(negvec,10);
for i = 1:length(negvec)
    if negvec(i,1) > neg1p
        negvec(i,1) = 0;
    end
end

pos1pm = reshape(posvec,100,100);
neg1pm = reshape(negvec,100,100);

EC1p = pos1pm+neg1pm;
figure; imagesc(EC1p);