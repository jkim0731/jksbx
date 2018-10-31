% # of trials, # of touches in each angle in 2 different volumes
% Considering just unimodal (or monotonic) tuning for now 2018/06/27 JK

% load u and, if it exists, ANOVA tuning file
clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
        
for mi = 1 : length(mice)
    for si = 1 : length(sessions{mi})
        u = Uber.buildUberArray(mice(mi), sessions{mi}(si));
    end
end
