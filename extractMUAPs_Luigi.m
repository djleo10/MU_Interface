function [MUAP, MUAPrms, MUAPp2p, MUAPduration] = extractMUAPs_Luigi(SIG,MUPulses,fsamp)
%new added luigi
cutWindowEdge = 170; %measure peak to peak up to 170
%figure,hold on,
for r=1:size(SIG,1)
    for c=1:size(SIG,2)
        if ~isempty(SIG{r,c})
            
            tmpY = SIG{r,c}; % SD (DD when SIG is SD)
            tmpMUAP=extractMUAPocurrences(tmpY,MUPulses,round(0.03*fsamp));
            MUAP{r,c} = mean(tmpMUAP{1});
            %old
%             MUAPrms(r,c) = sqrt(mean(MUAP{r,c}.^2));
%             MUAPp2p(r,c) = max(MUAP{r,c}) - min(MUAP{r,c});
            %new
            %plot(MUAP{r,c}),xline(cutWindowEdge,'r');
            MUAPrms(r,c) = sqrt(mean(MUAP{r,c}(1:cutWindowEdge).^2));
            MUAPp2p(r,c) = max(MUAP{r,c}(1:cutWindowEdge)) - min(MUAP{r,c}(1:cutWindowEdge));
            
            ind = find(abs(MUAP{r,c}) >= 0.8 * max(abs(MUAP{r,c})));
            MUAPduration(r,c) = max(ind) - min(ind);
        end
    end
end

function MUAP=extractMUAPocurrences(Y,MU_IND,WinSemiLength,winMode,TransSec)
% Rutine for spike triggered extraction of MUAPs from multichannel measurements
% MUAP=extractMUAPocurrences(Y,MU_IND,WinSemiLength,winMode,TransSec)
%
% INPUTS:
%       Y - row-array of all the measurements (each row corresponds to one measurement)
%       MU_IND - MU's discharge times (in samples)
%       WinSemiLength - semilength of window averaging function
%       winMode - window type (optional): if not specified, rectangular window is used
%       TransSec - optional parameter for the Tukey window
% OUTPUTS:
%       MUAP - cell structure containing the row-array of all MUAP ocurrences for the specific measurements
%              MUAP{i} for i-th measurement (i.e. i-th row in Y).

if iscell(Y)
    sigLength = 0;
    for r = 1:size(Y,1)
        for c = 1:size(Y,2)
            if ~isempty(Y{r,c})
                sigLength = max([sigLength,length(Y{r,c})]);
            end
        end
    end
    MU_IND=MU_IND(find( MU_IND>WinSemiLength+1 & MU_IND<sigLength-WinSemiLength-1));
else
    MU_IND=MU_IND(find( MU_IND>WinSemiLength+1 & MU_IND<length(Y)-WinSemiLength-1));
end
rY=size(Y,1);

if (nargin > 3)
    WinSemiLength=WinSemiLength+mod(WinSemiLength,2); % make WinSemiLength even
    if (strcmp('BLACKMAN',upper(winMode)))
        wind=blackman(WinSemiLength+1)'; % Use Blackman window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('BLACKMAN-HARRIS',upper(winMode)))
        wind=blackmanharris(WinSemiLength+1)'; % Use Blackman-Harris window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('NUTTALLWIN',upper(winMode)))
        wind=nuttallwin(WinSemiLength+1)'; % Use nutall window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('HANN',upper(winMode)))
        wind=hann(WinSemiLength+1)'; % Use Hanning window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('CHEBWIN',upper(winMode)))
        wind=chebwin(WinSemiLength+1,100)'; % Use Chebyshev window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('TUKEYWIN',upper(winMode)))
        wind=tukeywin(2*WinSemiLength+1,TransSec)'; % Use Tukey window
    else
        disp('ERROR: window not found');
        return;
    end
else
    wind=ones(1,2*WinSemiLength+1);
end

if isempty(MU_IND)
    for k=1:rY
        MUAP{k}=zeros(1,2*WinSemiLength+1);
    end
else
    
    if ~iscell(Y)
        for k=1:rY
            tmp = cutMUAP(MU_IND,WinSemiLength,Y(k,:));
            for k2 = 1:size(tmp,1)
                MUAP{k}(k2,:) = wind.*tmp(k2,:);
            end
        end
    else
        for r=1:size(Y,1)
            for c=1:size(Y,2)
                if ~isempty(Y{r,c})
                    tmp = cutMUAP(MU_IND,WinSemiLength,Y{r,c});
                    for k2 = 1:size(tmp,1)
                        MUAP{r,c}(k2,:) = wind.*tmp(k2,:);
                    end
                end
            end
        end
    end
    
end

function MUAPs=cutMUAP(MUPulses,len,Y)
% extracts cosecutive MUAPs out of signal Y and stores them row-wise in the output variable MUAPs
% INPUTS:
%   - MUPulses      triggering positions (in samples) for rectangular window used in extraction of MUAPs (MU firring patterns);
%   - len           radius of rectangular window (window length = 2*len+1)
%   - Y             single signal channel (raw vector containing a single channel of a recorded signals)
%
% OUTPUTS:
%   - MUAPs         row-wise matrix of extracted MUAPs (aligned signal intervals of length 2*len+1)

while ~isempty(MUPulses) && MUPulses(end)+2*len>length(Y)
    %disp(['Error: cutMUAP, Index exceeds the length of the array: ' num2str(MUPulses(end)) '+' num2str(len) ' > ' num2str(length(Y)) ]);
    MUPulses = MUPulses(1:end-1);
end
c=length(MUPulses);

edgeLen = round(len/5);
tmp = gausswin(2*edgeLen)';
win = [tmp(1:edgeLen) ones(1,2*len-2*edgeLen+1) tmp(edgeLen+1:end) ];
MUAPs=zeros(c,1+2*len);
for k=1:c
    MUAPs(k,1:1+2*len)= win .*[ zeros(1,max(MUPulses(k)-len,1)-(MUPulses(k)-len)) Y(max(MUPulses(k)-len,1):min(MUPulses(k)+len,length(Y))) zeros(1,MUPulses(k)+len-min(MUPulses(k)+len,length(Y)))];
end

