function jksbxwbdcb(src,callbackdata)

global bgimg mask
global cell_list cellpoly status segmenttool_h data mode_h
global plane_num im_plane info patch_h cm
global zimg patch_zh nm

action_type = get(gcbf, 'SelectionType');
m = findobj(segmenttool_h,'tag','method');
lb = findobj(segmenttool_h,'tag','listbox1');
lb.Max = 2;

if(mode_h.Value ==1) % when in segment mode
    if strcmp(action_type, 'normal') % left click        
        if info.volscan         
            p = gca;
            z = round(p.CurrentPoint);
            z = z(1,1:2);

            if(z(1)>0 && z(2)>0 && z(1)< size(bgimg.CData,2) && z(2)<size(bgimg.CData,1))
                switch m.Value
                    case 1
                        bw = bgimg.CData(:,:,3);                                 
                        if ~isempty(intersect(sub2ind(size(bw), find(bw ~= 0)),sub2ind(size(mask{im_plane}),find(mask{im_plane}~=0))))
                            status.String = sprintf('ROI overlapping: cannot draw');
                            drawnow;
                        else
                            B = bwboundaries(bw);
                            xy = B{1};
                            hold(bgimg.Parent,'on');
                            if isempty(cell_list)
                                cell_list = 1;
                                ncell = 1;
                            else
                                for ii = 1 : max(cell_list)+1
                                    if isempty(find(cell_list == ii,1))
                                        ncell = ii;
                                        if ii == 1
                                            cell_list = [1, cell_list];
                                        elseif ii == (max(cell_list) + 1)
                                            cell_list = [cell_list, ii];
                                        else
                                            cell_list = [cell_list(1:ii-1), ii, cell_list(ii:end)];
                                        end
                                        break;
                                    end
                                end
                            end
                            lb.String = cell_list;
                            patch_h{ncell} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.7,'edgecolor',[1 1 1],'parent',bgimg.Parent,'FaceLighting','none');
%                             cellpoly{ncell} = patch_h{ncell};
                            [xx,yy] = meshgrid(1:size(mask{im_plane},2),1:size(mask{im_plane},1));
                            bw=(inpolygon(xx,yy,xy(:,2),xy(:,1)));
                            mask{im_plane}(bw) = ncell;
                            hold(bgimg.Parent,'off');

                            % rescale...
                            r = single(bgimg.CData(:,:,1));
                            r(mask{im_plane}>0) = 0;
                            r = (r-min(r(:)))/(max(r(:))-min(r(:)));
                            bgimg.CData(:,:,1) = uint8(r*255);
                            status.String = sprintf('Segmented %d cells',length(cell_list));
                            
                            %do the same for zimg
                            hold(zimg.Parent,'on');
                            zimg.CData(:,:,1) = uint8(255*nm);
                            patch_zh{ncell} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.2,'edgecolor',[1 1 1],'parent',zimg.Parent,'FaceLighting','none');
                            hold(zimg.Parent,'off');
                            drawnow;
                        end
                    case 2
                        global nhood_h
                        nh = str2double(nhood_h.String);
                        img = data(z(2)-nh:z(2)+nh,z(1)-nh:z(1)+nh,:);
                        N = size(img,1);
                        [xx,yy] = meshgrid(1:N,1:N);
                        img = reshape(img,[N*N size(img,3)]);
                        lambda = 150;
                        img(:,end+1) = lambda*xx(:);
                        img(:,end+1) = lambda*yy(:);
                        [uu,~,~] = svd(img,0);
                        idx = kmeans(uu(:,1:8),4);
                        idx = reshape(idx,[N N]);

                        idx = (idx==idx(nh+1,nh+1));
                        idx = imclearborder(gather(idx));
                        idx = bwdistgeodesic(idx,nh+1,nh+1);
                        idx = imopen(isfinite(idx),strel('disk',1));

                        if(sum(idx(:))>40)

                            [ii,jj] = find(idx);
                            ii = ii+z(2)-(nh+1); jj = jj+z(1)-(nh+1);
                            z = zeros(size(bgimg.CData,1),size(bgimg.CData,2),'uint8');
                            idx = sub2ind(size(z),ii,jj);
                            z(idx) = 1;
                            if ~isempty(intersect(sub2ind(size(bw), find(bw ~= 0)),sub2ind(size(mask{im_plane}),find(mask{im_plane}~=0))))
                                status.String = sprintf('ROI overlapping: cannot draw');
                                drawnow;
                            else
                                bgimg.CData(:,:,3) = z;

                                B = bwboundaries(z);
                                xy = B{1};
                                hold(bgimg.Parent,'on');                                
                                patch_h{ncell} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.7,'edgecolor',[1 1 1],'parent',bgimg.Parent,'FaceLighting','none');
                                if isempty(cell_list)
                                    cell_list = 1;
                                    ncell = 1;
                                else
                                    for ii = 1 : max(cell_list)+1
                                        if isempty(find(cell_list == ii,1))
                                            ncell = ii;
                                            if ii == 1
                                                cell_list = [1, cell_list];
                                            elseif ii == (max(cell_list) + 1)
                                                cell_list = [cell_list, ii];
                                            else
                                                cell_list = [cell_list(1:ii-1), ii, cell_list(ii:end)];
                                            end
                                            break;
                                        end
                                    end                            
                                end
                                lb.String = cell_list;
%                                 cellpoly{ncell} = patch_h{ncell};                                
                                mask{im_plane}(z==1) = ncell;

                                % rescale...
                                r = single(bgimg.CData(:,:,1));
                                r(mask{im_plane}>0) = 0;
                                r = (r-min(r(:)))/(max(r(:))-min(r(:)));
                                bgimg.CData(:,:,1) = uint8(r*255);
                                status.String = sprintf('Segmented %d cells',length(cell_list));
                                hold(bgimg.Parent,'off');
                                
                                %do the same for zimg
                                hold(zimg.Parent,'on');
                                zimg.CData(:,:,1) = uint8(255*nm);
                                patch_zh{ncell} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.2,'edgecolor',[1 1 1],'parent',zimg.Parent,'FaceLighting','none','userdata',ncell);
                                hold(zimg.Parent,'off');
                                drawnow;
                            end
                        end
                end
            end        
        else        
            p = gca;
            z = round(p.CurrentPoint);
            z = z(1,1:2);
            m = findobj(segmenttool_h,'tag','method');
            if(z(1)>0 && z(2)>0 && z(1)< size(bgimg.CData,2) && z(2)<size(bgimg.CData,1))
                switch m.Value
                    case 1
                        bw = bgimg.CData(:,:,3);
                        B = bwboundaries(bw);
                        if ~isempty(intersect(sub2ind(size(bw), find(bw ~= 0)),sub2ind(size(mask),find(mask~=0))))
                            status.String = sprintf('ROI overlapping: cannot draw');
                            drawnow;
                        else
                            xy = B{1};
                            hold(bgimg.Parent,'on');                            
                            if isempty(cell_list)
                                cell_list = 1;
                                ncell = 1;
                            else
                                for ii = 1 : max(cell_list)+1
                                    if isempty(find(cell_list == ii,1))
                                        ncell = ii;
                                        if ii == 1
                                            cell_list = [1, cell_list];
                                        elseif ii == (max(cell_list) + 1)
                                            cell_list = [cell_list, ii];
                                        else
                                            cell_list = [cell_list(1:ii-1), ii, cell_list(ii:end)];
                                        end
                                        break;
                                    end
                                end
                            end
                            lb.String = cell_list;
                            patch_h{ncell} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.7,'edgecolor',[1 1 1],'parent',bgimg.Parent,'FaceLighting','none');
%                             cellpoly{ncell} = patch_h{ncell};
                            [xx,yy] = meshgrid(1:size(mask,2),1:size(mask,1));
                            bw=(inpolygon(xx,yy,xy(:,2),xy(:,1)));
                            mask(bw) = ncell;

                            % rescale...
                            r = single(bgimg.CData(:,:,1));
                            r(mask>0) = 0;
                            r = (r-min(r(:)))/(max(r(:))-min(r(:)));
                            bgimg.CData(:,:,1) = uint8(r*255);
                            status.String = sprintf('Segmented %d cells',length(cell_list));
                            hold(bgimg.Parent,'off');
                            
                            %do the same for zimg
                            hold(zimg.Parent,'on');
                            zimg.CData(:,:,1) = uint8(255*nm);
                            patch_zh{ncell} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.2,'edgecolor',[1 1 1],'parent',zimg.Parent,'FaceLighting','none');
                            hold(zimg.Parent,'off');
                            drawnow;
                        end
                    case 2
                        global nhood_h
                        nh = str2double(nhood_h.String);
                        img = data(z(2)-nh:z(2)+nh,z(1)-nh:z(1)+nh,:);
                        N = size(img,1);
                        [xx,yy] = meshgrid(1:N,1:N);
                        img = reshape(img,[N*N size(img,3)]);
                        lambda = 150;
                        img(:,end+1) = lambda*xx(:);
                        img(:,end+1) = lambda*yy(:);
                        [uu,~,~] = svd(img,0);
                        idx = kmeans(uu(:,1:8),4);
                        idx = reshape(idx,[N N]);

                        idx = (idx==idx(nh+1,nh+1));
                        idx = imclearborder(gather(idx));
                        idx = bwdistgeodesic(idx,nh+1,nh+1);
                        idx = imopen(isfinite(idx),strel('disk',1));

                        if(sum(idx(:))>40)

                            [ii,jj] = find(idx);
                            ii = ii+z(2)-(nh+1); jj = jj+z(1)-(nh+1);
                            z = zeros(size(bgimg.CData,1),size(bgimg.CData,2),'uint8');
                            idx = sub2ind(size(z),ii,jj);
                            z(idx) = 1;
                            if ~isempty(intersect(sub2ind(size(bw), find(bw ~= 0)),sub2ind(size(mask),find(mask~=0))))
                                status.String = sprintf('ROI overlapping: cannot draw');
                                drawnow;
                            else
                                bgimg.CData(:,:,3) = z;

                                B = bwboundaries(z);
                                xy = B{1};
                                hold(bgimg.Parent,'on');                                
                                if isempty(cell_list)
                                    cell_list = 1;
                                    ncell = 1;
                                else
                                    for ii = 1 : max(cell_list)+1
                                        if ~isempty(find(cell_list == ii,1))
                                            ncell = ii;
                                            if ii == 1
                                                cell_list = [1, cell_list];
                                            elseif ii == (max(cell_list) + 1)
                                                cell_list = [cell_list, ii];
                                            else
                                                cell_list = [cell_list(1:ii-1), ii, cell_list(ii:end)];
                                            end
                                            break;
                                        end
                                    end
                                end
                                patch_h{ncell} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.7,'edgecolor',[1 1 1],'parent',bgimg.Parent,'FaceLighting','none');
%                                 cellpoly{ncell} = patch_h{ncell};
                                hold(bgimg.Parent,'off');
                                mask(z==1) = ncell;
                                lb.String = cell_list;

                                % rescale...
                                r = single(bgimg.CData(:,:,1));
                                r(mask>0) = 0;
                                r = (r-min(r(:)))/(max(r(:))-min(r(:)));
                                bgimg.CData(:,:,1) = uint8(r*255);
                                status.String = sprintf('Segmented %d cells',length(cell_list));
                                
                                %do the same for zimg
                                hold(zimg.Parent,'on');
                                zimg.CData(:,:,1) = uint8(255*nm);
                                patch_zh{ncell} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.2,'edgecolor',[1 1 1],'parent',zimg.Parent,'FaceLighting','none');
                                hold(zimg.Parent,'off');
                                drawnow;
                            end
                        end
                end
            end
        end    
    elseif strcmp(action_type, 'alt') % right click
        if info.volscan            
            p = gca;
            z = round(p.CurrentPoint);
            z = z(1,1:2);
            if(z(1)>0 && z(2)>0 && z(1)< size(bgimg.CData,2) && z(2)<size(bgimg.CData,1))
                ncell = mask{im_plane}(z(2),z(1));
                if ncell == 0
                    status.String = sprintf('cannot select roi');
                else
                    lb.Value = find(cell_list == ncell);
                    answer = questdlg(sprintf('Delete this ROI? (%d)',ncell));
                    if strcmp(answer, 'Yes')
                        cell_list_idx = find(cell_list == ncell);
                        cell_list = cell_list(setdiff([1:length(cell_list)],cell_list_idx));                    
%                         lb.String = '';
                        lb.String = cell_list;
                        mask{im_plane}(sub2ind(size(mask{im_plane}),find(mask{im_plane} == ncell))) = 0;
                        delete(patch_h{ncell})
                        delete(patch_zh{ncell})
                        bgimg.CData(:,:,1) = uint8(255*nm); 
                    end
                    status.String = sprintf('Segmented %d cells',length(cell_list));
                    lb.Value = [];
                    drawnow;
                end
            end
       
        else
            p = gca;
            z = round(p.CurrentPoint);
            z = z(1,1:2);
            if(z(1)>0 && z(2)>0 && z(1)< size(bgimg.CData,2) && z(2)<size(bgimg.CData,1))
                ncell = mask(z(2),z(1));
                if ncell == 0
                    status.String = sprintf('cannot select roi');
                else
                    lb.Value = find(cell_list == ncell);
                    answer = questdlg(sprintf('Delete this ROI? (%d)',ncell));
                    if strcmp(answer, 'Yes')
                        cell_list_idx = find(cell_list == ncell);
                        cell_list = cell_list(setdiff([1:length(cell_list)],cell_list_idx)); 
    %                         lb.String = '';
                        lb.String = cell_list;
                        mask(sub2ind(size(mask),find(mask == ncell))) = 0;      
                        delete(patch_h{ncell})
                        delete(patch_zh{ncell})
                        bgimg.CData(:,:,1) = uint8(255*nm);
                    end
                    status.String = sprintf('Segmented %d cells',length(cell_list));
                    lb.Value = [];
                    drawnow;
                end
            end
        end    
    end
end
