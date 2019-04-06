% # of trials, # of touches in each angle in 2 different volumes
% Considering just unimodal (or monotonic) tuning for now 2018/06/27 JK

% load u and, if it exists, ANOVA tuning file
clear
mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};  
        
for mi = 6 : length(mice)
    for si = 1 : length(sessions{mi})
        u = Uber.buildUberArray(mice(mi), sessions{mi}(si));
    end
end
