function [m,disp] = sbxalign_nonrigid(fname,idx)

switch length(idx)
    case 1      % align one frame 
        A = gpuArray(squeeze(sbxread(fname,idx(1),1)));   % just one frame... easy!
        m = A;
        disp = {zeros([size(A) 2],'gpuArray')};
    case 2      % only two frames
        A = gpuArray(squeeze(sbxread(fname,idx(1),1)));   % read the frames
        B = gpuArray(squeeze(sbxread(fname,idx(2),1)));
        [D,Ar] = imregdemons(A,B,[32 16 8 4],'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4,'DisplayWaitBar',false);
        m = (Ar/2+B/2);
        disp = {D zeros([size(A) 2],'gpuArray')};     
    otherwise   % recursion
        I{1} = idx(1:floor(end/2));             % split dataset in two
        I{2} = idx(floor(end/2)+1 : end);       % recursive alignment
        D = cell(1,2);
        A = cell(1,2);
        parfor m=1:2
            [A{m},D{m}] = sbxalign_nonrigid(fname,I{m});
        end
        [Dnew,Ar] = imregdemons(A{1},A{2},[32 16 8 4],'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4,'DisplayWaitBar',false);
        m = (Ar/2+A{2}/2);
        D{1} = cellfun(@(x) (x+Dnew),D{1},'UniformOutput',false);  % concatenate distortions
        disp = [D{1} D{2}];
end
