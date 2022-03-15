function temp = subsampleFASK(BOLD,j,i,oAlpha,alpha)
        % select subsamples
        nobs = size(BOLD,1);
        BOLDsignal = datasample(BOLD,nobs/2,1,'Replace',false);


        % write the BOLD signal out as a txt file
        writematrix(BOLDsignal,sprintf('BOLDsignal%d%d',j,i))

        % run FASK algorithm
        command = sprintf('java  -jar causal-cmd-1.3.0-jar-with-dependencies.jar --algorithm fask --test fisher-z-test --orientationAlpha %d --faskAdjacencyMethod 1 --alpha %d --data-type continuous --dataset BOLDsignal%d%d.txt --delimiter comma --no-header --prefix FASK%d%d', oAlpha,alpha,j,i,j,i);
        system(command)
        temp = Tetrad2Matrix(sprintf('FASK%d%d.txt',j,i),'directed')';
end