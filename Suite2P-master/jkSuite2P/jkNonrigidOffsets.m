% computes registration offsets for data split into blocks
% loops over blocks and returns offsets dsall
function [dsall,ops1] = jkNonrigidOffsets(data, ops, ops1)

% Don't bother aligning across planes for now JK 2018/02/23
% alignAcrossPlanes  = getOr(ops, {'alignAcrossPlanes'}, false);
% planesToInterpolate = getOr(ops, {'planesToInterpolate'}, 1:nplanes);

nblocks = ops.numBlocks(1)*ops.numBlocks(2);

dsall = zeros(size(data,3), 2, nblocks,'double');
Corr = zeros(size(data,3), nblocks,'double');
for ib = 1:nblocks
    % collect ds
    ops1.mimg = ops1.mimgB{ib};
    if ops.kriging
      [dsall(:,:,ib), Corr(:,ib)]  = ...
          regoffKriging(data(ops1.yBL{ib},ops1.xBL{ib},:),ops1, 0);
    else
      [dsall(:,:,ib), Corr(:,ib)]  = ...
          regoffLinear(data(ops1.yBL{ib},ops1.xBL{ib},:),ops1,0);
    end
end
% Don't know what this is for, and it is allocated to dsall for
% nonrigid but not for rigid... what is going on? This is for the very
% first frame of the whole frames combined, so not going to affect that
% much. Ignore for now.
%     if j==1
%         ds(1,:,:) = 0;
%     end

ops1.DS                 = cat(1, ops1.DS, dsall);
ops1.CorrFrame          = cat(1, ops1.CorrFrame, Corr);

end
    