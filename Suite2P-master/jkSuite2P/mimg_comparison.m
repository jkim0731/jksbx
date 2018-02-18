base_dir = 'D:\TPM\JK\suite2p_results_tif\';
mouse = '036';
plane = 8;
sessions = 1:21;
mimgs = cell(1,length(sessions));
for i = 1 : length(sessions)
    session_dir = sprintf('%s%s\\%s_%03d_000\\%d',base_dir,mouse,mouse,sessions(i),plane);
    cd(session_dir)
    load_fn = dir('F*.mat');
    load(load_fn(1).name)
    mimgs{i} = ops.mimg;
end

%%
figure('units','normalized','outerposition',[0 0 0.7 0.7]), subplot(221), ...
    imagesc(mimgs{4}, [min(min(mimgs{4}(:,:))),max(max(mimgs{4}(:,:)))*0.7]), axis image, axis off, title('S01')
%     imagesc(mimgs{4}), axis image, axis off, title('S04')
j = 5;
for k = 1 : 3
    if j <= length(mimgs)
        subplot(2,2,k+1), imagesc(mimgs{j}, [min(min(mimgs{j}(:,:))),max(max(mimgs{j}(:,:)))*0.7]),...
%         subplot(2,2,k+1), imagesc(mimgs{j}),...            
            axis image, axis off, title(sprintf('S%02d',j))
    end
    j = j + 1;
    if j == 1
        j = j + 1;
    end
end

%%

base_dir = 'D:\TPM\JK\suite2p_results_tif\';
mouse = '027';
planes = 4:-1:1;
sessions = [3,16];
mimg = cell(2,4);
for i = 1 : 2    
    for j = 1 : 4
        session_dir = sprintf('%s%s\\%s_%03d_000\\%d',base_dir,mouse,mouse,sessions(i),planes(j));
        cd(session_dir)
        load_fn = dir('F*.mat');
        load(load_fn(1).name)
        mimgs{i,j} = ops.mimg;
    end
end
%%
for j = 1 : 4
    figure('units','normalized','outerposition',[0 0 0.3 0.7]),
    subplot(2,1,1), imagesc(mimgs{1,j}), axis image, axis off, title(sprintf('S%02d plane %d', sessions(1),4-planes(j)+1))
    subplot(2,1,2), imagesc(mimgs{2,j}), axis image, axis off, title(sprintf('S%02d plane %d', sessions(2),4-planes(j)+1))
end