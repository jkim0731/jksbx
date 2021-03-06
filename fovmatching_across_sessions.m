% Matching ROIs between different sessions
% Into a reference FOV (usually, naive fine angle session)
% Using uber array. there are averaged FOV images (mimg{x}) and cell map
% (cellmap{x}) in different planes.
% save roi match to a separate files.
% struct: 
% roiMatch.mouseName
% roiMatch.refSession
% roiMatch.sessionName(i)
% roiMatch.toRef{i}
% roiMatch.loss{i}
% roiMatch.gain{i}

% updates:
% 2019/10/10. (1) match histograms between images. (2) Constraints on the amount of rotation, scaling, and
% translation, (3) apply affine -> rigid 

clear
close all
baseDir = 'Y:\Whiskernas\JK\suite2p\';
mice = [25,27,30,36,39,52];
sessions = {[4,19,22],[3,10,17],[3,21,22],[1,17,18],[1,23,24],[3,21,26]}; % each second one is the reference session (expert session)

% mice = 27;
% sessions = {[3,8,14]}; % each second one is the reference session (expert session)
refSessionInd = 2;
[optimizer, metric] = imregconfig('monomodal');
overlapThreshold = 0.1;
% for mi = 1 : length(mice)
for mi = 2
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseDir, mouse))
    
    numSessions = length(sessions{mi}); % including ref session
    if numSessions > 1
        savefn = sprintf('cellIDmatchTest_JK%03d',mouse);
        us = struct;    
        for si = 1 : numSessions
            session = sessions{mi}(si);
            ufn = sprintf('UberJK%03dS%02d', mouse, session);
            load(ufn, 'u')
            us.sessions(si).mimg = u.mimg;
            us.sessions(si).cellmap = u.cellmap;
            us.sessions(si).sessionName = u.sessionName;
            us.sessions(si).cellID = u.cellNums;
            if mouse == 27
                for pi = 1 : 4
                    us.sessions(si).mimg{pi} = [];
                    us.sessions(si).cellmap{pi} = [];
                end
                us.sessions(si).cellID(find(us.sessions(si).cellID < 5000)) = [];
            end
            if si == refSessionInd % reference session
                us.refSession = u.sessionName;
            end
        end
        us.mouseName = u.mouseName;
        

        for si = 1 : 2 : numSessions
            session = sessions{mi}(si);
            % assume all # of planes are the same
            numPlanes = length(us.sessions(refSessionInd).mimg);
            
            tform = cell(numPlanes,1);
            for pi = 1 : numPlanes
                fprintf('Processing JK%03d S%02d plane #%d\n', mouse, session, pi)
                if ~isempty(us.sessions(refSessionInd).mimg{pi})
                    ref = mat2gray(us.sessions(refSessionInd).mimg{pi});
                    moving = mat2gray(us.sessions(si).mimg{pi});
                    ref2 = adapthisteq(ref);
                    moving2 = imhistmatch(moving, ref2);
                    tform{pi} = imregtform(moving2, ref2, 'affine', optimizer, metric);
                    T = tform{pi}.T;
                    if abs(1-T(1,1)) > 0.1 || abs(1-T(2,2)) > 0.1 || ...
                            abs(T(1,2)) > 0.1 || abs(T(2,1)) > 0.1 || ...
                            abs(T(3,1)) > 20 || abs(T(3,2)) > 20

                        moving2 = adapthisteq(moving);
                        ref2 = imhistmatch(ref2, moving);
                        tform{pi} = imregtform(moving2, ref2, 'affine', optimizer, metric);
                    end
                    T = tform{pi}.T;
                    if abs(1-T(1,1)) > 0.1 || abs(1-T(2,2)) > 0.1 || ...
                            abs(T(1,2)) > 0.1 || abs(T(2,1)) > 0.1 || ...
                            abs(T(3,1)) > 20 || abs(T(3,2)) > 20

                        ref2 = adapthisteq(ref);
                        moving2 = imhistmatch(moving, ref2);
                        tform{pi} = imregtform(moving2, ref2, 'rigid', optimizer, metric);
                    end
                    T = tform{pi}.T;
                    if abs(1-T(1,1)) > 0.1 || abs(1-T(2,2)) > 0.1 || ...
                            abs(T(1,2)) > 0.1 || abs(T(2,1)) > 0.1 || ...
                            abs(T(3,1)) > 20 || abs(T(3,2)) > 20

                        moving2 = adapthisteq(moving);
                        ref2 = imhistmatch(ref, moving2);
                        tform{pi} = imregtform(moving2, ref2, 'rigid', optimizer, metric);
                    end
                    T = tform{pi}.T;
                    if abs(1-T(1,1)) > 0.1 || abs(1-T(2,2)) > 0.1 || ...
                            abs(T(1,2)) > 0.1 || abs(T(2,1)) > 0.1 || ...
                            abs(T(3,1)) > 20 || abs(T(3,2)) > 20

                        moving2 = adapthisteq(moving);
                        ref2 = adapthisteq(ref);
                        tform{pi} = imregtform(moving2, ref2, 'rigid', optimizer, metric);
                    end 
                    T = tform{pi}.T;
                    if abs(1-T(1,1)) > 0.1 || abs(1-T(2,2)) > 0.1 || ...
                            abs(T(1,2)) > 0.1 || abs(T(2,1)) > 0.1 || ...
                            abs(T(3,1)) > 20 || abs(T(3,2)) > 20

                        moving2 = moving;
                        ref2 = ref;
                        tform{pi} = imregtform(moving2, ref2, 'rigid', optimizer, metric);
                    end
    %                 if abs(1-T(1,1)) > 0.1 || abs(1-T(2,2)) > 0.1 || ...
    %                         abs(1-T(1,2)) > 0.1 || abs(1-T(2,1)) > 0.1
    %                     ref2 = ref;
    %                     moving2 = moving;
    %                     tform{pi} = imregtform(moving2, ref2, 'rigid', optimizer, metric);
    %                 end
    %                 T = tform{pi}.T;
    %                 if abs(1-T(1,1)) > 0.1 || abs(1-T(2,2)) > 0.1 || ...
    %                         abs(1-T(1,2)) > 0.1 || abs(1-T(2,1)) > 0.1
    %                     ref2 = adapthisteq(ref);
    %                     moving2 = adapthisteq(moving);
    %                     tform{pi} = imregtform(moving2, ref2, 'translation', optimizer, metric);
    %                 end
    %                 T = tform{pi}.T;
    %                 if abs(1-T(1,1)) > 0.1 || abs(1-T(2,2)) > 0.1 || ...
    %                         abs(1-T(1,2)) > 0.1 || abs(1-T(2,1)) > 0.1
    %                     ref2 = ref;
    %                     moving2 = moving;
    %                     tform{pi} = imregtform(moving2, ref2, 'translation', optimizer, metric);
    %                 end
                end
            end
            us.sessions(si).tform = tform;
%             figure
            us.sessions(si).matchedRefCellID = zeros(length(us.sessions(si).cellID),1);
            refCellmap = cell(numPlanes,1);
            movedCellmap = cell(numPlanes,1);
            for pi = 1 : numPlanes
                if ~isempty(us.sessions(refSessionInd).cellmap{pi})
                
                    refCellmap{pi} = us.sessions(refSessionInd).cellmap{pi};
                    movingCellmap = us.sessions(si).cellmap{pi};
                    movedCellmap{pi} = imwarp(movingCellmap, tform{pi}, 'OutputView', imref2d(size(ref)));
                    subplot(3,3,pi)
                    imshowpair(refCellmap{pi}, movedCellmap{pi})

                    movCurrInd = find(floor(us.sessions(si).cellID/1000) == pi);
                    movCellID = us.sessions(si).cellID(movCurrInd);
                    matchingID = zeros(length(movCellID),1);

                    refCurrInd = find(floor(us.sessions(refSessionInd).cellID/1000) == pi);
                    refCellID = us.sessions(refSessionInd).cellID(refCurrInd);

                    overlapMat = zeros(length(movCurrInd), length(refCurrInd));
                    for movi = 1 : length(movCurrInd)
                        % pixel indices for each cell ID's of moved cell map
                        movPixi = find(movedCellmap{pi} == movCellID(movi));
                        for refi = 1 : length(refCurrInd)
                            % pixel indices for each cell ID's of reference cell map
                            refPixi = find(refCellmap{pi} == refCellID(refi));
                            % calculate the overlap index, and allocate to the matrix
                            tempOverlap = length(intersect(movPixi, refPixi)) / length(union(movPixi, refPixi));
                            if tempOverlap > overlapThreshold
                                overlapMat(movi, refi) = tempOverlap;
                            end
                        end
                    end

                    [~, overlapInd] = sort(overlapMat(:), 'descend');                
                    overlapInd = overlapInd(1:length(find(overlapMat(:))));

                    [overlapi, overlapj] = ind2sub(size(overlapMat), overlapInd);
                    for jj = 1 : length(overlapj)
                        if overlapj(jj)
                            matchingID(overlapi(jj)) = refCellID(overlapj(jj));
                            tobeDeletedInd = find(overlapj == overlapj(jj));
                            overlapj(tobeDeletedInd) = 0;
                        end
                    end

                    % checking multiple allocations
                    mIDlist = setdiff(unique(matchingID), 0);
                    for midi = 1 : length(mIDlist)
                        if length(find(matchingID == mIDlist(midi))) > 1
                            error('More than one cell allocated')
                        end
                    end
                    us.sessions(si).matchedRefCellID(movCurrInd) = matchingID;
                end
            end
            
            figure
            for pi = 1 : numPlanes
                if ~isempty(refCellmap{pi})
                    currMovInd = find(floor(us.sessions(si).cellID/1000) == pi);
                    tempRefCellID = us.sessions(si).matchedRefCellID(currMovInd);
                    tempMovedCellID = us.sessions(si).cellID(currMovInd);
                    matchedInd = find(tempRefCellID);

                    selectedMoved = zeros(size(movedCellmap{pi}));
                    selectedRef = zeros(size(refCellmap{pi}));

                    for i = 1 : length(matchedInd)
                        selectedMoved(movedCellmap{pi}==tempMovedCellID(matchedInd(i))) = 1;
                        selectedRef(refCellmap{pi}==tempRefCellID(matchedInd(i))) = 1;
                    end

                    subplot(3,3,pi)
                    imshowpair(selectedRef, selectedMoved)
                end
            end
            drawnow
        end
        save(savefn, 'us')
    end
    
end

