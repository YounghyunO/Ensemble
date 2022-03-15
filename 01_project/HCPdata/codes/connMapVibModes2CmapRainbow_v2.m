%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Selen Atasoy
%
% maps the limits of a vibration mode to a color map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Selen Atasoy
%
% maps the limits of a vibration mode to a color map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cMap] = connMapVibModes2CmapRainbow_v2(cmin, cmax, resolution)


%%
if cmin<0
    if abs(cmin) < cmax && cmax>0

            % |--------|---------|--------------------|    
          % -cmax      cmin       0                  cmax         [cmin,cmax]


        nr_blues    = round(resolution* (abs(cmin) / (abs(cmin) + cmax))); 
        nr_reds     = round(resolution* (abs(cmax) / (abs(cmin) + cmax))); 



    elseif abs(cmin) >= cmax && cmax>0

         % |------------------|------|--------------|    
       %  cmin                0     cmax          -cmin         [cmin,cmax]
       %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))

        nr_blues    = round(resolution* (abs(cmin) / (abs(cmin) + cmax))); 
        nr_reds     = round(resolution* (abs(cmax) / (abs(cmin) + cmax))); 

    %%%% modified by Katharina Glomb, CHUV, March 2019 %%%%%
    elseif cmax<0
        nr_blues    = round(resolution);%* (abs(cmin) / (abs(cmin) + abs(cmax))));
        nr_reds = 0;

    end
elseif (cmin==0) && (cmax==0)
        
        nr_blues        = round(resolution*0.01);
        nr_reds         = round(resolution*0.01);

        
        
else
         % |------------------|------|--------------|    
       %                      0     cmin          -cmax         [cmin,cmax]
       %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))

        nr_blues    = round(resolution* (abs(cmin) / max(cmin, cmax))); 
        nr_reds     = round(resolution* (abs(cmax) / max(cmin, cmax))); 

end




if nr_blues == 0 
    cMap = black2yellow2(nr_reds);
elseif nr_reds == 0 
    cMap = cyan2black(nr_blues); 
else
%     if nr_blues>nr_reds
%         cMap = [cyan2black(nr_blues); black2yellow(nr_reds)];
%     else
%         cMap = [flipud(black2yellow(nr_reds)); flipud(cyan2black(nr_blues))];
%     end
     %%% modified by Katharina Glomb
    % old version: always uses the same range of colors
    cMap =[cyan2black(nr_blues); black2yellow(nr_reds)];
    % new version: adjusts the color range to the number of blues/reds
%     use_nr = max([nr_blues,nr_reds]);
%     cMap_blues = cyan2black(use_nr); 
%     cMap_reds = black2yellow(use_nr);
%     cMap = [cMap_blues((use_nr-nr_blues)+1:end,:);cMap_reds(1:nr_reds,:)];
end


function map = black2yellow2(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
   
% 'da'; 'da'; 'da'; %grey
% 'da'; 'da'; 'da'; %grey
% % 'da'; 'da'; 'da'; %grey
% % 'da'; 'da'; 'da'; %grey
'6a';'67';'67';
'6a';'67';'67';
'6a';'67';'67';
'6a';'67';'67';
'6a';'67';'67';
'6a';'67';'67';
% % 'ff'; '00'; '00';  % red
'ff'; '00'; '00';  % red
'ff'; '3a'; '00';  % red
'ff'; '5d'; '00';  % orange
'ff'; '70'; '00';  % orange
'ff'; '83'; '00';  % oran-yell
'ff'; '99'; '00';  % oran-yell
'ff';'ae';'00';
'ff'; 'b9'; '00';  % yellow
'ff'; 'c9'; '00';  % yellow
'ff'; 'ff'; '00'; 
'ff'; 'ff'; '00'; 
'ff'; 'ff'; '00'; 


% '2a';'2a';'86';
% '29';'39';'fe';
% '92';'97';'ff';
% '65';'2a';'a1';
% '44';'92';'94';
% '4c';'dc';'29';
% 'c1';'fa';'7d';
% '7f';'fa';'81';
% '7E';'F9';'06';
% 'F5';'1A';'00';
% '7F';'07';'00';
% 'F8';'A0';'C1';
% '42';'2E';'40';
% '1F';'1E';'1D';
% 'C2';'60';'22';
% '55';'FB';'C4';
% '81';'81';'81';
% '7E';'1C';'5D';
% '61';'41';'20';
% 'C2';'C2';'C3';
% '30';'20';'0C';
% 'E2';'C0';'C0';

]; 

% '00';'00';'00'; % black
% '66';'00';'00'; % dark red
% '66';'00';'00'; % dark red
% % 'ff';'00';'00'; % red
% 'ff';'00';'00'; % red
% 'ff'; '69'; '00'; % orange
% 'ff'; '69'; '00'; % orange
% 'ff'; '99'; '00'; % oran-yell
% 'ff'; '99'; '00'; % oran-yell
% % 'ed'; 'b1'; '20'; 
% % 'ed'; 'b1'; '20'; 
% 'ff'; 'ff'; '00';
% 'ff'; 'ff'; '00';

% % 'bb'; 'bb'; 'bb';
% '00'; '00'; '00';  % black
% % '66'; '00'; '33';  % purple2
% '66'; '00'; '00' ; % dark red
% % '66'; '00'; '00' ; % dark red
% % 'ff'; '38'; '8d';  % hotpink
% 'ff'; '00'; '00';  % red
% % 'ff'; '00'; '00';  % red
% 'ff'; '69'; '00';  % orange
% 'ff'; '69'; '00';  % orange
% 'ff'; '99'; '00';  % oran-yell
% 'ff'; '99'; '00';  % oran-yell
% 'ff'; 'ff'; '00';  % yellow
% 'ff'; 'ff'; '00'];  % yellow

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');


function map = cyan2black(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...

% '10'; 'b0'; '10';  % limegreen
% '00'; 'ff'; '00';  % green
'00'; 'ff'; 'ff';  % cyan
'00'; 'ff'; 'ff';  % cyan

% '00'; 'ff'; '00';  % green
% '00'; 'ff'; '00';  % green
% '10'; 'b0'; '10';  % limegreen
% '10'; 'b0'; '10';  % limegreen

'00'; 'cc'; 'ff'; % light blue
'00'; 'cc'; 'ff'; % light blue
'33'; '99'; 'ff'; % blue
'33'; '99'; 'ff'; % blue

% 'e2'; '51'; 'e2';  % violet
% %'ff'; '38'; '8d';  % hotpink
% '66'; '00'; '33';  % purple2
% '00'; '00'; '66';  % blue
% '00'; '00'; '66';  % blue
%'7f'; '7f'; 'cc';  % blue_videen7
%'4c'; '4c'; '7f';  % blue_videen9
%'33'; '33'; '4c';  % blue_videen11
%'00'; '00'; '00'];  % black
%'ff'; '00'; '00';  % red
% 'ff'; '69'; '00';  % orange
% 'ff'; '99'; '00';  % oran-yell
% 'ff'; 'ff'; '00';  % yellow
% 
% 
% 'ff'; 'ff'; 'ff';  % white
% 'dd'; 'dd'; 'dd';  % gry-dd
% 'bb'; 'bb'; 'bb';  % gry-bb
'6a';'67';'67';
'6a';'67';'67';
]; % black

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');

function map = black2yellow(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'6a';'67';'67';
'6a';'67';'67';  % black
%'66'; '00'; '33';  % purple2
'66'; '00'; '00' ; % dark red
'66'; '00'; '00' ; % dark red
%'ff'; '38'; '8d';  % hotpink
'ff'; '00'; '00';  % red
'ff'; '00'; '00';  % red
'ff'; '69'; '00';  % orange
'ff'; '69'; '00';  % orange
'ff'; '99'; '00';  % oran-yell
'ff'; '99'; '00';  % oran-yell
'ff'; 'ff'; '00';  % yellow
'ff'; 'ff'; '00'];  % yellow

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');




