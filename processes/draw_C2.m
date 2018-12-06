function draw_C2(fn)
mouse = fn(1:3);
if str2double(mouse) < 50
    diffim = show_C2_old(fn);
else
    diffim = show_C2(fn);
end
%
figure, imagesc(diffim{4}), axis image, axis off, hold on
% %%
% cd([sbxdir, mouse])
% totalDiffim = cell(8,1);
% for i = 1 : 8
%     fn = sprintf('036_9998_10%d',i);
%     totalDiffim{i} = show_C2_old(fn);
% end
% 
% %%
% temp = zeros(size(totalDiffim{1}{1}));
% for i = 1 : 8
%     for j = 3
%         temp = temp + totalDiffim{i}{j}/8;
%     end
% end
% figure, imagesc(temp, [0 0.4]), axis image, axis off, hold on


%%

xpoints = [];
ypoints = [];
while true
    [x, y] = ginput(1);
    if isempty(x) % enter
        break
    end
    xpoints = [xpoints, x];
    ypoints = [ypoints, y];
    plot(xpoints, ypoints, 'r-', 'linewidth', 5)
end
plot([xpoints, xpoints(1)], [ypoints, ypoints(1)], 'r-', 'linewidth', 5)

% cd([suite2pdir,mouse])
savefn = sprintf('JK%sC2.mat',mouse);
save(savefn, 'xpoints', 'ypoints')

