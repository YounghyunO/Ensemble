%% script for significance test
nullECs = [];

for i = 1:1000
    rDCM = cell2mat(rDCMs(randi(202)));
    VAR = cell2mat(VARs(randi(202)));
    GC = cell2mat(GCs(randi(202)));
    FASK = cell2mat(FASKs(randi(202)));
    
    Data = [rDCM(:) VAR(:) GC(:) FASK(:)];
    
    temp = predict(EnsembleMdl,Data);
    EC = reshape(temp,100,100);
    nullECs = cat(3,nullECs,EC);
end

save('nullEC.mat','nullEcs')

figure; hist(squeeze(nullECs(1,1,:)),1000)

av_EC = mean(nullECs,3);
