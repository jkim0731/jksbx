% # of trials, # of touches in each angle in 2 different volumes
% Considering just unimodal (or monotonic) tuning for now 2018/06/27 JK

% load u and, if it exists, ANOVA tuning file
% clear
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19,22],[3,16,17],[3,21,22],[1,17,18],[7],[2],[1,23,24],[3],[3,21,26],[3],[3],[3]};
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19,22],[3,10,17],[3,21,22],[1,17,18],[7],[2],[1,23,24],[3],[3,21,26],[3],[3],[3]};
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

for mi = 1 : length(mice)
    for si = 1 : length(sessions{mi})
        u = Uber.buildUberArray(mice(mi), sessions{mi}(si));
    end
end
