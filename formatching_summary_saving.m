baseDir = 'Y:\Whiskernas\JK\suite2p\';
mice = [25,27,30,36,39,52];
match = cell(length(mice),3); % 1: naive cID, 2: matched ref cID, 3: ref cID
for mi = 1 : length(mice)
    load(sprintf('%s%03d\\cellIDmatch_JK%03d',baseDir,mice(mi),mice(mi)), 'us')
    match{mi,1} = us.sessions(1).cellID;
    match{mi,2} = us.sessions(1).matchedRefCellID;
    match{mi,3} = us.sessions(2).cellID;
end
%%
saveFn = 'cellMatching_beforeNafter';
save(sprintf('%s%s', baseDir, saveFn), 'match')
