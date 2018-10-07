function u = buildUberArray(mouse, session)

mouseName = sprintf('JK%03d', mouse);
sessionName = sprintf('S%02d',session);
caSessionName = sprintf('%03d',session);

savefn = sprintf('Uber%s%s.mat',mouseName, sessionName); 

bDirBase = 'C:\Data\SoloData\';
bDir = [bDirBase, mouseName];

wDirBase = 'C:\Data\WhiskerVideo\';
wDir = [wDirBase, mouseName, sessionName];

caDirBase = 'C:\Data\suite2p\';
caDir = sprintf('%s%03d\\',caDirBase, mouse);

cd(caDir)
% if exist(savefn, 'file')
%     load(savefn)
% else
    if exist('b','var') && iscell(b)
        if isprop(b{1}, 'mouseName') && strcmp(b{1}.mouseName, mouseName)
            if exist('bSession','var')
                if isprop(bSession, 'sessionName')
                    if strcmp(bSession.sessionName, sessionName)
                        disp('running with existing bSession')
                    else
                        disp('finding right bSession')
                        bSession = b{cellfun(@(x) strcmp(x.sessionName,sessionName), b)};
                    end
                end
            else
                disp('assigning bSession')
                bSession = b{cellfun(@(x) strcmp(x.sessionName,sessionName), b)};
            end
        else
            disp('building new behavior array')
            cd(bDir)
            load(['behavior_',mouseName]) % loading b{}
            bSession = b{cellfun(@(x) strcmp(x.sessionName,sessionName), b)};
        end
    else
        disp('building new behavior array')
        cd(bDir)
        load(['behavior_',mouseName]) % loading b{}
        bSession = b{cellfun(@(x) strcmp(x.sessionName,sessionName), b)};
    end

    if exist('wlArray', 'var') && strcmp(wlArray.mouseName, mouseName) && strcmp(wlArray.sessionName, sessionName)
        disp('using existing wlArray')
    else
        disp('building new whisker trial lite array')
        cd(wDir)
        wlArray = Whisker.WhiskerTrialLite_2padArray(mouseName, sessionName);
    end

    if exist('caArray', 'var') && strcmp(caArray.mouseName, mouseName(3:end)) && strcmp(caArray.sessionName, caSessionName)
        disp('using existing caArray')
    else
        disp('building new calcium array')
        cd(caDir)
        caArray = Calcium.CalciumDataArray(mouse,session);
    end

    if strcmp(bSession.mouseName, wlArray.mouseName) && strcmp(bSession.mouseName, wlArray.mouseName) ...
            && strcmp(bSession.sessionName, wlArray.sessionName) && strcmp(bSession.sessionName, ['S', caArray.sessionName(2:3)])
        u = Uber.Uber_2padArray(bSession, wlArray, caArray);
    else
        disp('mouseName or sessionName mismatch')
    end

    save(savefn, 'u')
% end