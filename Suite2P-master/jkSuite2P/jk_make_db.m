% UPDATE Christmas 2016: number of clusters determined automatically, but
% do specify the "diameter" of an average cell for best results. You can do this with either
% db(iexp).diameter, or ops0.diameter. 
mice = {'025'};
sessions = 3;
i = 0;


siglist = [0.2,0.5,1];
nsvdlist = [100,300,1000];
nframeslist = [1000,5000,10000];
diameterlist = [10, 20, 40];

for mi = 1 : length(mice)
    for si = 1 : length(sessions)  
        for siglisti = 1 : length(siglist)
            for nsvdlisti = 1 : length(nsvdlist)
                for nframeslisti = 1 : length(nframeslist)
                    for dlisti = 1 : length(diameterlist)
        i = i + 1;
        db(i).mouse_name    = mice{mi};
        db(i).session       = sessions(si);      
        db(i).nonrigid      = 0;
        db(i).doDeconvolution = 0;
        db(i).ResultsSavePath = sprintf('D:\TPM\JK\suite2ptest\celltest\%02d\',i);
        db(i).signalExtraction = 'surround';
        
        
        
        db(i).sig                    = siglist(siglisti);  % spatial smoothing length in pixels; encourages localized clusters
        db(i).nSVDforROI             = nsvdlist(nsvdlisti); % how many SVD components for cell clustering
        db(i).NavgFramesSVD          = nframeslist(nframeslisti); % how many (binned) timepoints to do the SVD based on
        db(i).diameter               = diameterlist(dlisti);
                    end
                end
            end
        end
    end
end

% mice = {'027'};
% sessions = 0:25;
% for mi = 1 : length(mice)
%     for si = 1 : length(sessions)
%         i = i + 1;
%         db(i).mouse_name    = mice{mi};
%         db(i).session       = sessions(si);        
%     end
% end

% mice = {'039'};
% sessions = [1:28,901];
% for mi = 1 : length(mice)
%     for si = 1 : length(sessions)        
%         i = i + 1;
%         db(i).mouse_name    = mice{mi};
%         db(i).session       = sessions(si);        
%     end
% end




%% examples
% db(i).mouse_name    = '039';
% db(i).session       = [3,4,5,6,9,10,11,18,19];        
% db(i).targetSession = 4; %targetSession determines which session to be reference for registration
% 

% %
% i = i+1;
% db(i).mouse_name    = 'M150329_MP009';
% db(i).date          = '2015-04-10';
% db(i).expts         = [5 6 7 8 9 10 11];
% db(i).diameter      = 12;
% 
% i = i+1;
% db(i).mouse_name    = 'M150824_MP019';
% db(i).date          = '2015-12-19';
% db(i).expts         = [4];
% db(i).diameter      = 6;
% 
% % example of datasets, which consist of several sessions - use cell arrays
% % will be treated as subsets of experiment with the same FOV, with
% % different names/dates (for one reason or another), analyzed together
% i = i+1;
% db(i).mouse_name    = {'MK020', 'M150416_MK020'};
% db(i).date          = {'2015-07-30', '2015-07-30'};
% db(i).expts         = {[2010 2107], [1 2 3]};
% db(i).diameter      = 12;
% 
% % example for datasets without folder structure
% db(i).mouse_name    = 'notImportant';
% db(i).date          = '2016';
% db(i).expts         = []; % leave empty, or specify subolders as numbers
% db(i).diameter      = 12;
% db(i).RootDir       = 'F:\DATA\neurofinder\neurofinder.01.00\images'; % specify full path to tiffs here
% 
% % example extra entries
% % db(i).AlignToRedChannel= 1;
% % db(i).BiDiPhase        = 0; % adjust the relative phase of consecutive lines
% % db(i).nSVD             = 1000; % will overwrite the default, only for this dataset
% % db(i).comments      = 'this was an adaptation experiment';
% % db(i).expred        = [4]; % say one block which had a red channel 
% % db(i).nchannels_red = 2; % how many channels did the red block have in total (assumes red is last)
