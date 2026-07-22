close all
clear
clc

%% paths
dataPath = 'C:\Data\Demuse Data\Luigi\SMA03\W3_20230608\Decomposition Output\';
%% load data
[files_name,files_path] = uigetfile(strcat(dataPath,'\','*.mat'));
decOutPath = (strcat(string(files_path),string(files_name)));
load(decOutPath);
%% global variables
fsamp = 4000;
singleDifferential = true;
%% collect MUAPs
MUAPs = []; MUAPrms = []; MUAPp2p = []; MU_ID = [];
%% Prepare signal
%filter SIG low-pass at 20 Hz and High-pass at 500 Hz
for r=1:size(SIG,1)
    for c=1:size(SIG,2)
        try
            SIG{r,c} = lpbutter(SIG{r,c},500,fsamp);
            SIG{r,c} = hpbutter(SIG{r,c},20,fsamp);
        catch
        end
    end
end

%flip signal if needed

%flipSIG= flip(flip(SIG,2),1); %match with grid pic SMA03 bicep fem
flipSIG = SIG; %match with DEMUSE --> no flip

%perform single differential if needed
if singleDifferential
    newSig = cell(size(flipSIG,1)-1,size(flipSIG,2));
    %SD across row
    for r=1:size(flipSIG,1)-1
        for c=1:size(flipSIG,2)
            if ~isempty(flipSIG{r,c}) && ~isempty(flipSIG{r+1,c})
                newSig{r,c} = flipSIG{r,c} - flipSIG{r+1,c};
            end
        end
    end

    %SD across column
    %     for r=1:size(flipSIG,1)
    %         for c=1:size(flipSIG,2)-1
    %             if ~isempty(flipSIG{r,c+1}) && ~isempty(flipSIG{r,c})
    %                 newSig{r,c} = flipSIG{r,c} - flipSIG{r,c+1};
    %             end
    %         end
    %     end

else
    newSig = flipSIG;
end

%% Extract MUAPs
%remove empy cells
emptyCells = cellfun(@isempty,MUPulses);
MUPulses(emptyCells) = [];
MUIDs(emptyCells) = [];

for m = 1:length(MUPulses) %  for all the MUs
    [MUAPs{end+1}, MUAPrms{end+1}, MUAPp2p{end+1}, ~] = extractMUAPs_Luigi(newSig,MUPulses{m},fsamp);

    if length(MUIDs{m}) < 10
        MU_ID{end+1} = MUIDs{m};
    else
        MU_ID{end+1} = 'MUID error';
    end
end
%% plot MUAPs
plotAllMUAPs(MUAPs,MU_ID,fsamp,0,[1,1,1],[0,0,0]);
set(gcf, 'Renderer', 'painters');
%saveas(gcf,'PE_MUAPS.svg');
%% Cross-Correlation Analysis for a Specific Motor Unit
% Assuming you have a specific MU in mind and it has been identified
targetMU = 'MU5'; % Example MU ID
targetIndex = find(strcmp(MU_ID, targetMU));  % Find the index of target MU

% Check if the MU ID is found and process
if isempty(targetIndex)
    disp(['Motor Unit ' targetMU ' not found.']);
else
    % Assuming that 'newSig' holds the processed signal matrices for each MU
    % Extract the signal matrix for the target MU
    targetSignals = newSig{targetIndex};

    % Number of channels (columns in the signal matrix for the target MU)
    numChannels = size(targetSignals, 2);

    % Initialize matrix to store maximum correlation coefficients
    maxCorrCoeffs = zeros(numChannels, numChannels);

    % Calculate cross-correlation for each pair of channels
    for i = 1:numChannels
        for j = i+1:numChannels % Only compare each pair once
            [cc, lags] = xcorr(targetSignals(:, i), targetSignals(:, j), 'coeff');
            % Find the maximum correlation coefficient
            maxCorrCoeffs(i, j) = max(abs(cc));
            fprintf('Max correlation between channels %d and %d of %s is %f\n', i, j, targetMU, maxCorrCoeffs(i, j));
        end
    end
end

%% prepare subplot territory
%define rows and columns
switch length(MUAPs)
    case {1,2,3}
        R=1; C=length(MUAPs);
    case {4,5,6}
        R=2; C=ceil(length(MUAPs)/2);
    case {7,8,9}
        R=3; C=3;
    case {10,11,12}
        R=3; C=4;
    case {13,14,15}
        R=3; C=5;
    case {16,17,18}
        R=3; C=6;
    case {19,20}
        R=4; C=5;
    case {21,22,23,24}
        R=4; C=6;
    case {25}
        R=5; C=5;
    case {26,27,28,29,30}
        R=5; C=6;
    case {31,32,33,34,35}
        R=5; C=7;
    case {36,37,38,39,40}
        R=5; C=8;
    otherwise
        R=8; C=8;
end
%% convex hull measurements
thrs = 20/100; %20% source
[X,Y] = meshgrid(1:7,1:8);
grid_crds = [reshape(X,56,1) reshape(Y,56,1)];
fillColor = [0.96,0.76,0.56];

iMUAPp2p = MUAPp2p;
iMUAPs = MUAPs;
iNMUs = length(iMUAPs);
area = nan(1,length(iMUAPs));
areaCh = nan(1,length(iMUAPs));

h = figure;
for j = 1:iNMUs
    subplot(R,C,j),
    plot(grid_crds(:,2),grid_crds(:,1),'b*'),hold on, title(strcat('MU', " ",MU_ID{j}));
    % flip stats to match with MUAP figures
    jMUAP = flip(iMUAPs{j},1);
    jMUAPp2p = flip(iMUAPp2p{j},1);
    jMUAPrms = flip(iMUAPp2p{j},1);
    [sourceValue,sourcePos] = max(jMUAPp2p,[],'all','linear');
    %jnormMUAPp2p = jMUAPp2p./sourceValue;

    %identify MUAPs above threshold
    [sRows,sColumns] = find(jMUAPp2p >= thrs*sourceValue);
    sRows = 1+(size(jMUAPp2p,1)-sRows); %space plot same as matrix
    plot(sColumns,sRows,'r*');
    %set(gca,'YTick',1:size(jMUAPp2p,1),'YTickLabel',flip(1:size(jMUAPp2p,1)));
    set(gca,'XLim',[0,9],'YLim',[0,8]);

    %convex hull
    P = [sColumns,sRows];
    [k,av] = convhull(P);
    %[k,av] = boundary(P);
    plot(P(k,1),P(k,2));
    fill( P(k,1),P(k,2),fillColor,'FaceColor',fillColor,'EdgeColor','none','FaceAlpha',0.5);
    area(j) = av;
    areaCh(j) = length(sRows);
    text(6.5,7.5,strcat('Area = ',num2str(av)));
end
set(h,'Color',[1,1,1],'WindowState','maximized');


