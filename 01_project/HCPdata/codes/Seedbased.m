%% seed-based analysis
%afferent to default
N = size(Ensembles,1);
vistodm = cell(N,1); smtodm = cell(N,1);
datodm = cell(N,1); vatodm = cell(N,1);
limtodm = cell(N,1); fptodm = cell(N,1);
dmtodm = cell(N,1);
%efferent from default
visfromdm = cell(N,1); smfromdm = cell(N,1);
dafromdm = cell(N,1); vafromdm = cell(N,1);
limfromdm = cell(N,1); fpfromdm = cell(N,1);
dmfromdm = cell(N,1);

for i = 1:N
    EC = cell2mat(Ensembles(i));
    EC(abs(EC)<0.03) = 0;
    %calculate afferent
    vistodm{i} = EC([38:50,90:100],[1:9,51:58]);
    smtodm{i} = EC([38:50,90:100],[10:15,59:66]);
    datodm{i} = EC([38:50,90:100],[16:23,67:73]);
    vatodm{i} = EC([38:50,90:100],[24:30,74:78]);
    limtodm{i} = EC([38:50,90:100],[31:33,79:80]);
    fptodm{i} = EC([38:50,90:100],[34:37,81:89]);
    dmtodm{i} = EC([38:50,90:100],[38:50,90:100]);
    
    %calculate efferent
    visfromdm{i} = EC([1:9,51:58],[38:50,90:100]);
    smfromdm{i} = EC([10:15,59:66],[38:50,90:100]);
    dafromdm{i} = EC([16:23,67:73],[38:50,90:100]);
    vafromdm{i} = EC([24:30,74:78],[38:50,90:100]);
    limfromdm{i} = EC([31:33,79:80],[38:50,90:100]);
    fpfromdm{i} = EC([34:37,81:89],[38:50,90:100]);
    dmfromdm{i} = EC([38:50,90:100],[38:50,90:100]);
end
    
%concat across subjets
vistodmall = []; smtodmall = [];
datodmall = []; vatodmall = [];
limtodmall = []; fptodmall = [];
dmtodmall = [];

visfromdmall = []; smfromdmall = [];
dafromdmall = []; vafromdmall = [];
limfromdmall = []; fpfromdmall = [];
dmfromdmall = [];

for i = 1:N
    vistodmall = cat(2,vistodmall,cell2mat(vistodm(i)));
    smtodmall = cat(2,smtodmall,cell2mat(smtodm(i)));
    datodmall = cat(2,datodmall,cell2mat(datodm(i)));
    vatodmall = cat(2,vatodmall,cell2mat(vatodm(i)));
    limtodmall = cat(2,limtodmall,cell2mat(limtodm(i)));
    fptodmall = cat(2,fptodmall,cell2mat(fptodm(i)));
    dmtodmall = cat(2,dmtodmall,cell2mat(dmtodm(i)));
    
    visfromdmall = cat(2,visfromdmall,cell2mat(visfromdm(i)));
    smfromdmall = cat(2,smfromdmall,cell2mat(smfromdm(i)));
    dafromdmall = cat(2,dafromdmall,cell2mat(dafromdm(i)));
    vafromdmall = cat(2,vafromdmall,cell2mat(vafromdm(i)));
    limfromdmall = cat(2,limfromdmall,cell2mat(limfromdm(i)));
    fpfromdmall = cat(2,fpfromdmall,cell2mat(fpfromdm(i)));
    dmfromdmall = cat(2,dmfromdmall,cell2mat(dmfromdm(i)));
end

vistodmpos = nonzeros(max(vistodmall,0));
vistodmneg = nonzeros(-min(vistodmall,0));

smtodmpos = nonzeros(max(smtodmall,0));
smtodmneg = nonzeros(-min(smtodmall,0));

datodmpos = nonzeros(max(datodmall,0));
datodmneg = nonzeros(-min(datodmall,0));

vatodmpos = nonzeros(max(vatodmall,0));
vatodmneg = nonzeros(-min(vatodmall,0));

limtodmpos = nonzeros(max(limtodmall,0));
limtodmneg = nonzeros(-min(limtodmall,0));

fptodmpos = nonzeros(max(fptodmall,0));
fptodmneg = nonzeros(-min(fptodmall,0));

dmtodmpos = nonzeros(max(dmtodmall,0));
dmtodmneg = nonzeros(-min(dmtodmall,0));


data = [vistodmpos;vistodmneg;smtodmpos;smtodmneg;datodmpos;datodmneg;...
    vatodmpos;vatodmneg;limtodmpos;limtodmneg;fptodmpos;fptodmneg;...
    dmtodmpos;dmtodmneg;];

visn = size([vistodmpos;vistodmneg],1);
smn = size([smtodmpos;smtodmneg],1);
dan = size([datodmpos;datodmneg],1);
van = size([vatodmpos;vatodmneg],1);
limn = size([limtodmpos;limtodmneg],1);
fpn = size([fptodmpos;fptodmneg],1);
dmn = size([dmtodmpos;dmtodmneg],1);

color1 = cell(visn,1);
origin1 = cell(visn,1);
for i = 1:visn
    origin1(i) = {'Visual'};
end
for i = 1:size(vistodmpos,1)
    color1{i} = {'Positive'};
end
for i = size(vistodmpos,1)+1:size(vistodmpos,1)+size(vistodmneg,1)
    color1{i} = {'negative'};
end


origin2 = cell(smn,1);
color2 = cell(smn,1);
for i = 1:smn
    origin2(i) = {'Somatomotor'};
end
for i = 1:size(smtodmpos,1)
    color2{i} = {'Positive'};
end
for i = size(smtodmpos,1)+1:size(smtodmpos,1)+size(smtodmneg,1)
    color2{i} = {'negative'};
end


origin3 = cell(dan,1);
color3 = cell(dan,1);
for i = 1:dan
    origin3(i) = {'Dorsal attention'};
end
for i = 1:size(datodmpos,1)
    color3{i} = {'Positive'};
end
for i = size(datodmpos,1)+1:size(datodmpos,1)+size(datodmneg,1)
    color3{i} = {'negative'};
end


origin4 = cell(van,1);
color4 = cell(van,1);
for i = 1:van
    origin4(i) = {'Ventral attention'};
end
for i = 1:size(vatodmpos,1)
    color4{i} = {'Positive'};
end
for i = size(vatodmpos,1)+1:size(vatodmpos,1)+size(vatodmneg,1)
    color4{i} = {'negative'};
end



origin5 = cell(limn,1);
color5 = cell(limn,1);
for i = 1:limn
    origin5(i) = {'Limbic'};
end
for i = 1:size(limtodmpos,1)
    color5{i} = {'Positive'};
end
for i = size(limtodmpos,1)+1:size(limtodmpos,1)+size(limtodmneg,1)
    color5{i} = {'negative'};
end


origin6 = cell(fpn,1);
color6 = cell(fpn,1);
for i = 1:fpn
    origin6(i) = {'Fronto Parietal'};
end
for i = 1:size(fptodmpos,1)
    color6{i} = {'Positive'};
end
for i = size(fptodmpos,1)+1:size(fptodmpos,1)+size(fptodmneg,1)
    color6{i} = {'negative'};
end


origin7 = cell(dmn,1);
color7 = cell(dmn,1);
for i = 1:dmn
    origin7(i) = {'Default'};
end
for i = 1:size(dmtodmpos,1)
    color7{i} = {'Positive'};
end
for i = size(dmtodmpos,1)+1:size(dmtodmpos,1)+size(dmtodmneg,1)
    color7{i} = {'negative'};
end


origins = [origin1;origin2;origin3;origin4;origin5;origin6;origin7];
colors = [color1;color2;color3;color4;color5;color6;color7];
result = [array2table(data),cell2table(colors),cell2table(origins)];


