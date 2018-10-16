%% Load u and singleCell_F_ANOVA first
%% Load ops1 too. ('regops_xxx' file)

% using JK025S19 plane 1
% JK025S04 plane 5 

% zoom 2x. 0.7 um per pix


v = VideoWriter('test');
v.FrameRate = u.frameRate; 
open(v)
colors = {[0.2 0.2 0.7], [0.2 0.7 0.2], [0.7 0.2 0.2]};
colorsPane = {[0.2 0.2 0.7], [0.2 0.7 0.2], [0.7 0.2 0.2]};
trialchoice = [7 7 6]; % [1 1 2];[2 4 3];
plane = 5; % or 5
paneWidth = 0;
headerHeight = 70;
fn = '027_016_000';
figure,
for grandI = 1 : 3
%% sort out indices
    angle = 45:45:135;
    depths =  269 : 165/3 : 269+165; 
    % 110 : 165/3 : 110+165; 244 : 165/3 : 244+165; % for JK025
    % 110 : 165/3 : 110+165,  269 : 165/3 : 269+165] for JK027
    umppix = 0.7;
    
    touchInd = find(cellfun(@(x) length(x.touchChunks), u.trials));
    tempDrinkInd = find(cellfun(@(x) x.touchChunks{1}(1) < 2, u.trials(touchInd)));
    tempPlaneInd = find(cellfun(@(x) ismember(plane,x.planes), u.trials));
    inds = cell(length(angle),1);
    for ai = 1 : length(angle)
        tempAngleInd = find(cellfun(@(x) x.angle == angle(ai), u.trials));
        inds{ai} = intersect(touchInd(tempDrinkInd), intersect(tempAngleInd, tempPlaneInd));
    end

    %% play one by one
    a = angle(grandI);
    i = trialchoice(grandI);

    
    % load([fn, '.align'], '-mat')

    tempinds = inds{angle==a};

    ut = u.trials{tempinds(i)};
    tnum = ut.trialNum;
    frames = ops1{1}.trials(find([ops1{1}.trials.trialnum] == tnum)).frames(1):ops1{1}.trials(find([ops1{1}.trials.trialnum] == tnum)).frames(2);

    if mod(frames(1),4)
        frames = frames(mod(frames(1),4):end);
    else
        frames = frames(4:end);
    end

    %%

    z = jksbxreadframes(fn, frames);
    z = z(:,101:end,1:end-mod(size(z,3),4));
    zz = zeros(headerHeight + 2*size(z,1)+paneWidth*4, 2*size(z,2)+paneWidth*4, 3, size(z,3)/4);
    for i = 1 : size(z,3)/4
        zz(1+headerHeight+paneWidth : headerHeight+size(z,1)+paneWidth, 1+paneWidth : paneWidth+size(z,2), :, i) = repmat((z(:,:,1+(i-1)*4)), [1 1 3]);
        zz(1+headerHeight+paneWidth : headerHeight+size(z,1)+paneWidth, 1+paneWidth+size(z,2)+paneWidth*2 : end-paneWidth, :, i) = repmat((z(:,:,2+(i-1)*4)), [1 1 3]);
        zz(1+headerHeight+paneWidth*3+size(z,1) : end-paneWidth, 1+paneWidth : paneWidth+size(z,2), :, i) = repmat((z(:,:,3+(i-1)*4)), [1 1 3]);
        zz(1+headerHeight+paneWidth*3+size(z,1) : end-paneWidth, 1+paneWidth+size(z,2)+paneWidth*2 : end-paneWidth, :, i) = repmat((z(:,:,4+(i-1)*4)), [1 1 3]);
    end
    zz(1:headerHeight,:,:,:) = deal(max(max(max(max(zz)))));
    %%
    dTextPositions = {[headerHeight+paneWidth+30, paneWidth+600], [headerHeight+paneWidth+30, paneWidth*3+size(z,2)+600], ...
        [headerHeight+paneWidth*3+30+size(z,1), paneWidth+600], [headerHeight+paneWidth*3+30+size(z,1), paneWidth*3+size(z,2)+600]};
    rulerMargin = 20;
    rulerLength = 100; % in um
    times = cell2mat(ut.touchChunks);
    linewidth = 5;
    tpmFrames = unique(round(times*u.frameRate));
%     if length(tpmFrames) < 3
%         tpmFrames(2) = tpmFrames(1)+1;
%         tpmFrames(3) = tpmFrames(1)+2;
%     end
    %%
    
    for i = 1 : size(zz,4)

        imshow(mat2gray(zz(:,:,:,i))), hold on    
        for j = 1 : 4
            text(dTextPositions{j}(2), dTextPositions{j}(1), [num2str(round(depths(j))),'\mum'], 'color', 'white', 'fontweight', 'bold', 'fontsize', 15)
        end
        plot([size(zz,2) - rulerMargin - rulerLength/umppix, size(zz,2) - rulerMargin], [size(zz,1)-rulerMargin, size(zz,1)-rulerMargin], 'linewidth', 3, 'color', [1 1 1])
        text(10, 20, [num2str(a), '\circ'], 'fontweight', 'bold', 'fontsize', 30, 'color', colors{grandI})
        if i == 1
            text(size(zz,2)-120, 20, '0 s', 'fontweight', 'bold', 'fontsize', 20, 'color', colors{grandI})
        else
            text(size(zz,2)-120, 20, [num2str(round(ut.tpmTime{1}(i-1),2)), ' s'], 'fontweight', 'bold', 'fontsize', 20, 'color', colors{grandI})
        end
        
        plot([200, size(zz,2)-200], [headerHeight - 10, headerHeight - 10], 'lineWidth', linewidth, 'color', colorsPane{grandI})
        
        for j = 1 : length(tpmFrames)
            temppoint = 200 + (size(zz,2) - 400) / size(zz,4) * tpmFrames(j);
            plot([temppoint, temppoint], [headerHeight - 10, headerHeight - 70], 'lineWidth', linewidth, 'color', colors{grandI})
        end
        timepoint = 200 + (size(zz,2) - 400)/size(zz,4) * i;
        plot([timepoint, timepoint], [headerHeight - 5, headerHeight - 60], 'lineWidth', linewidth, 'color', [0 0 0])
%         if ismember(i, tpmFrames)
%             plot([1,size(zz,2)], [1+headerHeight, 1+headerHeight], 'lineWidth', linewidth, 'color', colorsPane{grandI})
%             plot([1,size(zz,2)], [size(zz,1), size(zz,1)], 'lineWidth', linewidth, 'color', colorsPane{grandI})
%             plot([1,size(zz,2)], [1+headerHeight+size(z,1)+paneWidth, 1+headerHeight+size(z,1)+paneWidth], 'lineWidth', linewidth, 'color', colorsPane{grandI})
%             plot([1,1], [1+headerHeight, size(zz,1)], 'lineWidth', linewidth, 'color', colorsPane{grandI})
%             plot([1+size(z,2)+paneWidth*2,1+size(z,2)+paneWidth*2], [1+headerHeight, size(zz,1)], 'lineWidth', linewidth, 'color', colorsPane{grandI})
%             plot([size(zz,2),size(zz,2)], [1+headerHeight, size(zz,1)], 'lineWidth', linewidth, 'color', colorsPane{grandI})
%         end
        hold off
        drawnow;
        tempMov = getframe(gcf);
        writeVideo(v,tempMov)
    end    
end
close(v)


