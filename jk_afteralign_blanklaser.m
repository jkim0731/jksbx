mouse = {'648', '650''653'};
dir_base = 'd:\2p\';
for i = 1 : length(mouse)
    cd([dir_base,mouse{i},'\'])
    flist = dir('*.sbx');
    for j = 1 : length(flist)
        tic;
        fn = flist(j).name(1:end-4);
        disp(['processing ', dir_base, mouse{i},' ' fn, '.sbx'])
        
        sbxaligndir({fn});
        maxidx = jkget_maxidx(fn);
        flevel = zeros(1,maxidx);
        for k = 1 : maxidx
            flevel(k) = mean(mean(squeeze(sbxread(fn,(k-1),1))));
        end

        blankidx = find(flevel < min(flevel) + std(flevel));

        load([fn,'.align'],'-mat')
        T = T - repmat(T(1,:),size(T,1),1);
        for k = 1 : length(blankidx)
            T(blankidx(k),:) = [0 0];
        end
        tempim = squeeze(sbxread(fn,0,1));
        m = zeros(size(tempim));
        for k = 1 : maxidx
            m = m + double(circshift(squeeze(sbxread(fn,(k-1),1)),T(k,:)))/maxidx;
        end        
        timepassed = toc;
        disp(['done in ',num2str(timepassed),' seconds'])
        save([fn,'.align'],'T','m','-mat')
    end
end

% mouse = {'648', '650' ,'653'};
% dir_base = '\\whisker-nas\Storage\JK-temp\2p\';
% for i = 1 : length(mouse)
%     cd([dir_base,mouse{i},'\'])
%     flist = dir('*.sbx');
%     for j = 1 : length(flist)
%         tic;
%         fn = flist(j).name(1:end-4);
%         disp(['processing ', dir_base, mouse{i},' ' fn, '.sbx'])
%         
%         if ~exist([fn,'.align'],'file')
%             sbxaligndir({fn});
%             maxidx = jkget_maxidx(fn);
%             flevel = zeros(1,maxidx);
%             for k = 1 : maxidx
%                 flevel(k) = mean(mean(squeeze(sbxread(fn,(k-1),1))));
%             end
%             blankidx = find(flevel < min(flevel) + std(flevel));
% 
%             load([fn,'.align'],'-mat')
%             T = T - repmat(T(1,:),size(T,1),1);
%             for k = 1 : length(blankidx)
%                 T(blankidx(k),:) = [0 0];
%             end
%             tempim = squeeze(sbxread(fn,0,1));
%             m = zeros(size(tempim));
%             for k = 1 : maxidx
%                 m = m + double(circshift(squeeze(sbxread(fn,(k-1),1)),T(k,:)))/maxidx;                
%             end   
%             save([fn,'.align'],'T','m','-mat')
%             timepassed = toc;
%             disp(['done in ',num2str(ceil(timepassed/60)),' mins'])
%         else
%             disp(['skipping ', fn])
%         end
%     end
% end
% 
