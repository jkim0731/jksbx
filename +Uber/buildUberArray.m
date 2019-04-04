function u = buildUberArray(mouse, session)

mouseName = sprintf('JK%03d', mouse);
sessionName = sprintf('S%02d',session);
caSessionName = sprintf('%03d',session);

savefn = sprintf('Uber%s%s.mat',mouseName, sessionName); 

bDirBase = 'Y:\Whiskernas\JK\SoloData\';
% bDirBase = 'C:\Data\SoloData\';
bDir = [bDirBase, mouseName];

wDirBase = 'Y:\Whiskernas\JK\whisker\tracked\';
% wDirBase = 'C:\Data\WhiskerVideo\';
wDir = [wDirBase, mouseName, sessionName];

caDirBase = 'Y:\Whiskernas\JK\suite2p\';
% caDirBase = 'C:\Data\suite2p\';
caDir = sprintf('%s%03d\\',caDirBase, mouse);

cd(caDir)
% if exist(savefn, 'file')
%     load(savefn)
% else
%     if exist('b','var') && iscell(b)
%         if isprop(b{1}, 'mouseName') && strcmp(b{1}.mouseName, mouseName)
%             if exist('bSession','var')
%                 if isprop(bSession, 'sessionName')
%                     if strcmp(bSession.sessionName, sessionName)
%                         disp('running with existing bSession')
%                     else
%                         disp('finding right bSession')
%                         bSession = b{cellfun(@(x) strcmp(x.sessionName,sessionName), b)};
%                     end
%                 end
%             else
%                 disp('assigning bSession')
%                 bSession = b{cellfun(@(x) strcmp(x.sessionName,sessionName), b)};
%             end
%         else
%             disp('building new behavior array')
%             cd(bDir)
%             load(['behavior_',mouseName]) % loading b{}
%             bSession = b{cellfun(@(x) strcmp(x.sessionName,sessionName), b)};
%         end
%     else
%         disp('building new behavior array')
        cd(bDir)
        load(['behavior_',mouseName]) % loading b{}
        bSession = b{cellfun(@(x) strcmp(x.sessionName,sessionName), b)};
%     end


%     if exist('wfArray', 'var') && strcmp(wfArray.mouseName, mouseName) && strcmp(wfArray.sessionName, sessionName)
%         disp('using existing w3Array')
%     else
%         disp('building new whisker trial lite array')
        cd(wDir)
        wfArray = Whisker.WhiskerFinal_2padArray(wDir);
%     end
    
    
%     if exist('caArray', 'var') && strcmp(caArray.mouseName, mouseName(3:end)) && strcmp(caArray.sessionName, caSessionName)
%         disp('using existing caArray')
%     else
%         disp('building new calcium array')
        cd(caDir)
        caArray = Calcium.CalciumDataArray(mouse,session);
%     end
    
    
    if strcmp(bSession.mouseName, wfArray.mouseName) && strcmp(bSession.sessionName, wfArray.sessionName) && strcmp(bSession.sessionName, ['S', caArray.sessionName(2:3)])
        u = Uber.Uber_2padArray(bSession, wfArray, caArray);
    else
        disp('mouseName or sessionName mismatch')
    end
    
%     if strcmp(bSession.mouseName, w3Array.mouseName) && strcmp(bSession.mouseName, w3Array.mouseName) ...            
%             && strcmp(bSession.sessionName, w3Array.sessionName) && strcmp(bSession.sessionName, ['S', caArray.sessionName(2:3)])
%         u = Uber.Uber_2padArray(bSession, w3Array, caArray);
%     else
%         disp('mouseName or sessionName mismatch')
%     end

    save(savefn, 'u')
% end