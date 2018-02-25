% loops over data and computes registration offsets
% inputs:
%%% data = [ly, lx, nframes]
%%% j is which tiff
%%% iplane0 is first position of each plane in tiff
%%% ops and ops1 are options
function [dsall, ops1] = jkRigidOffsets(data, ops, ops1)

% split into subsets (for high scanning resolution recordings)
[xFOVs, yFOVs] = get_xyFOVs(ops);


dat = data(yFOVs(:,l),xFOVs(:,l),:);
if ~isempty(ops.smooth_time_space)
    dat = smooth_movie(dat, ops);
end
% align all loaded frames to target image of current plane
% (get registration offsets)
if ops.kriging
    [dsall, Corr]  = regoffKriging(dat, ops1, 0);
else
    [dsall, Corr]  = regoffLinear(dat, ops1, 0);
end

%ds          = RemoveBadShifts(ds);

% collect ds
%         if j==1
%             ds(1,:,:) = 0;
%         end
ops1.DS          = cat(1, ops1.DS, dsall);
ops1.CorrFrame   = cat(1, ops1.CorrFrame, Corr);

    
    
    % check if there was a sharp drop in fluorescence
%     lbright = sq(mean(data(:,:,indframes),2));
%     mlbright = mean(ops1{i}.mimg, 2);
%     
%     lbright = bsxfun(@rdivide, lbright, mlbright);
%     badi = max(abs(lbright(1:end-4,:) - lbright(5:end,:)), [], 1) > .5;
%     badi = find(badi);
%     
%     ops1{i}.badframes(sum(ops1{i}.Nframes) + badi) = true;    
end