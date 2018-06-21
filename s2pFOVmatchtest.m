mice = [25,27,30,36,37,38,39,41];
sessions = {[1:25],[0:25],[1:13,15,17:25],[901,1:21],[901,902,1:24],[901,1:31],[901,1:28],[901,1:30]};
FOVcorr = cell(length(mice),1);
for mi = 3 : length(mice)
    mouse = sprintf('%03d',mice(mi));
    cd(['D:\TPM\JK\',mouse])
    % initialization
    FOVcorr{mi} = cell(8,1);
    mimgs = cell(8,1);
    for pi = 1 : 8
        FOVcorr{mi}{pi} = zeros(length(sessions{mi}));
        mimgs{pi} = cell(length(sessions{mi}));
    end
    
    % collecting all mimg into mimgs
    for si = 1 : length(sessions{mi})
        session = sprintf('%03d',sessions{mi}(si));
        load(['regops_',mouse,'_',session]) % loading ops1, 8x1 cell
        for pi = 1 : 8
            mimgs{pi}{si} = ops1{pi}.mimg1;
        end
    end
    
    % calculate correlation between all the sessions of the same plane
    for pi = 1 : 8
        disp(['Plane ', num2str(pi), ' of mouse ', mouse, ' in process.'])
        for si = 1 : length(sessions{mi})
            refMimg = mimgs{pi}{si};
            for i = si : length(sessions{mi})
                targetMimg = mimgs{pi}{i};
                [u, v] = jkfftalign(targetMimg, refMimg);
                alignedMimg = circshift(targetMimg, [u, v]);
                if u >= 0
                    if v >= 0
                        ref = refMimg(u+1:end,v+1:end);
                        tar = alignedMimg(u+1:end,v+1:end);
                    else
                        ref = refMimg(u+1:end,1:end+v);
                        tar = alignedMimg(u+1:end,1:end+v);
                    end
                else
                    if v >= 0
                        ref = refMimg(1:end+u,v+1:end);
                        tar = alignedMimg(1:end+u,v+1:end);
                    else
                        ref = refMimg(1:end+u,1:end+v);
                        tar = alignedMimg(1:end+u,1:end+v);
                    end
                end
                if size(ref) == size(tar)
                    rho = corrcoef(ref(:),tar(:));                    
                    FOVcorr{mi}{pi}(si,i) = rho(2);
                end
            end
        end
    end
end
                        
                        
                        
                        
                        