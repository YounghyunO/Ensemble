function temp3 = subsample_wrapper(BOLD,i,oAlpha,alpha)
    FASKtemp = cell(100,1);
    parfor j = 1:100
        FASKtemp{j} = subsampleFASK(BOLD,j,i,oAlpha,alpha)
    end
    temp2 = [];
    for k = 1:100
        temp2 = cat(3,temp2,cell2mat(FASKtemp(k)));
    end
    temp3 = mean(temp2,3);
%     temp3(temp3<0.5) = 0;
    delete BOLD* FASK* FAS*
end