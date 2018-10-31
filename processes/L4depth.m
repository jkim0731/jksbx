zstep = mean(diff([knobbyInfo.z]))/cosd(35); % in um
zoom = 1;
pixpum = 1/(1.4/zoom);
%%
z = zeros(size(zstack),'like',zstack);
for i = 1 : size(zstack,3)
    z(:,:,i) = adapthisteq(zstack(:,:,i));
end
%%
% a = squeeze(mean(z(180:190,100:end,:),1));
% a = flip(a');
% binsize = round(pixpum * zstep);
% b = zeros([size(a,1), ceil(size(a,2)/binsize)]);
% for i = 1 : size(b,2)-1
%     b(:,i) = mean(a(:,(i-1)*binsize+1:i*binsize),2);
% end
% b(:,end) = mean(a(:,i*binsize+1:end),2);
%     
% figure, imshow(mat2gray(b)), axis on
%%
implay(z)

%%
implay(zstack)
%%
% for mi = 45:48
%     mouse = sprintf('%03d',mi);
%     cd(['D:\TPM\JK\',mouse])
%     tempflist = dir([mouse, '_*.sbx']);
%     sessionlist = zeros(length(tempflist),1);
%     for j = 1 : length(tempflist)
%         a = split(tempflist(j).name, '_');
%         sessionlist(j) = str2num(a{2});
%     end
%     sessionlist = unique(sessionlist);
%     for j = 1 : length(sessionlist)
%         make_zstack(sprintf('%s_%03d',mouse,sessionlist(j)));
%     end
% end
%%
septumFrames = 46:76;
numLines = 3;
pia = 258;
tempim = mean(zstack(:,:,septumFrames),3);
points = zeros(numLines,4); % x1 y1 x2 y2
figure, imshow(mat2gray(tempim)), hold on
lineProfiles = cell(3,1);
H = fspecial('gaussian');

for i = 1 : numLines
    for j = 1 : 2
        [points(i,(j-1)*2 + 1), points(i,(j-1)*2 + 2)] = ginput(1);
    end
    plot(points(i,[1,3]), points(i,[2,4]), 'r-')
    lineLength = round(sqrt((points(i,3) - points(i,1))^2 + (points(i,4) - points(i,2))^2));
    lineProfiles{i} = zeros(pia, lineLength);
    xbin = (points(i,3)-points(i,1))/lineLength;
    ybin = (points(i,4)-points(i,2))/lineLength;
    for j = 1 : pia
        currPlane = imfilter(zstack(:,:,j), H);
        for k = 1 : lineLength            
%             lineProfiles{i}(j,k) = zstack(round(points(i,2)+(k-1)*ybin), round(points(i,1)+(k-1)*xbin), j);
            lineProfiles{i}(j,k) = currPlane(round(points(i,2)+(k-1)*ybin), round(points(i,1)+(k-1)*xbin));
        end
    end
end

%%
stdVals = zeros(numLines,pia);
for i = 1:3
    line = zeros(size(lineProfiles{i}));
    for j = 1 : pia
        line(j,:) = smooth(lineProfiles{i}(j,:));
        line(j,:) = (line(j,:) - min(line(j,:)))/ (max(line(j,:)) - min(line(j,:)));
    end
    stdVals(i,:) = std(line,0,2);
end

signal = smooth(flip(mean(stdVals)));
figure, plot(1:size(stdVals,2), signal)
xlabel('Planes from the pia'), ylabel('Standard deviation')

