function varargout = jksbxsegmenttool(varargin)
% jksbxsegmenttool MATLAB code for jksbxsegmenttool.fig
%      jksbxsegmenttool, by itself, creates a new jksbxsegmenttool or raises the existing
%      singleton*.
%
%      H = jksbxsegmenttool returns the handle to a new jksbxsegmenttool or the handle to
%      the existing singleton*.
%
%      jksbxsegmenttool('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in jksbxsegmenttool.M with the given input arguments.
%
%      jksbxsegmenttool('Property','Value',...) creates a new jksbxsegmenttool or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before jksbxsegmenttool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to jksbxsegmenttool_OpeningFcn via guidevarargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help jksbxsegmenttool

% Last Modified by GUIDE v2.5 19-Jun-2017 18:25:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @jksbxsegmenttool_OpeningFcn, ...
                   'gui_OutputFcn',  @jksbxsegmenttool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before jksbxsegmenttool is made visible.
function jksbxsegmenttool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to jksbxsegmenttool (see VARARGIN)

% Choose default command line output for jksbxsegmenttool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes jksbxsegmenttool wait for user response (see UIRESUME)
% uiwait(handles.figure1);
axes(handles.axz);
global zimg hline vline
zimg = imagesc(zeros(512,796,3,'uint8'));
hold on
hline = plot([ 0 0],[0 0],'m--');
vline = plot([0 0 ],[0 0],'m--');
colormap gray
axis off
hold off

axes(handles.ax);
global bgimg
bgimg = imagesc(zeros(512,796,3,'uint8'));
colormap gray
axis off

end

% --- Outputs from this function are returned to the command line.
function varargout = jksbxsegmenttool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global bgimg data segmenttool_h nframes cell_list cellpoly mask rfn pathname
global plane_num im_plane info patch_h cm
global zimg hline vline patch_zh nm
global manual_h
handles.status.String = 'Resetting/Clearing GPU';
drawnow;
gpuDevice(1);

global sig;
delete(handles.axes3.Children);
sig = [];

[fn,pathname] = uigetfile({'*.sbx; *.align'});
if isequal(fn,0)
    return
end
if ~strcmpi(pwd,pathname(1:end-1))
    cd(pathname)
end

rfn = strtok(fn,'.');
% idx = max(strfind(rfn,'_'));
% rfnx = rfn(1 : (idx-1));
rfnx = rfn;

try
    load('-mat',[rfnx '.align']);     
catch
    return
end
axis off

% global info
z = sbxread(rfnx,0,1); 
if info.volscan
    plane_num = length(info.otwave_um);
    handles.popupmenu3.String = 1:plane_num;
else
    plane_num = 1;
    handles.popupmenu3.String = 1;
end
drawnow;

handles.status.String = 'Loading alignment data';

if iscell(m)
    im_plane = get(handles.popupmenu3,'Value');     
    m = m{im_plane};
    T = T{im_plane};
else
    im_plane = 1;
end

if(exist('mnr','var'))
    m = gather(mnr);
else
    m = double(m); % 2016/11/01 JK
end

nm = (m-min(m(:)))/(max(m(:))-min(m(:))); % normalized m
nm = adapthisteq(nm);
nm = single(nm);
% x = single(x);
nm = (nm-min(nm(:)))/(max(nm(:))-min(nm(:)));

bgimg = image(zeros(size(nm,1),size(nm,2),3,'uint8'));

bgimg.CData(:,:,1) = uint8(255*nm);
bgimg.CData(:,:,2) = bgimg.CData(:,:,1);
bgimg.CData(:,:,3) = bgimg.CData(:,:,1);

a = zimg.Parent;
idx = [];
for k = 1:length(a.Children)
    if(isa(a.Children(k),'matlab.graphics.primitive.Patch'))
        idx = [idx k];
    end
end
delete(a.Children(idx));

zimg.CData(:,:,1) = uint8(255*nm);
zimg.CData(:,:,2) = uint8(255*nm);
zimg.CData(:,:,3) = uint8(255*nm);
set(zimg.Parent,'xlim',[size(m,2)/2-32 size(m,2)/2+32],'ylim',[size(m,1)/2 - 32 size(m,1)/2+32]);
set(hline,'xdata',[size(m,2)/2-32 size(m,2)/2+32],'ydata',size(m,1)/2 *[ 1 1]);
set(vline,'ydata',[size(m,1)/2-32 size(m,1)/2+32],'xdata',size(m,2)/2 *[ 1 1]);

% if(~isempty(cellpoly))
%     cellfun(@delete,cellpoly);
% end

drawnow;

handles.status.String = 'Loading spatio-temporal data';

nframes = str2double(handles.frames.String);
if exist([rfnx,'.trials'],'file')
    try
        load('-mat',[rfnx '.trials']);     
    catch
        return
    end
    frames = [];
    for i = 1 : length(trials)
        frames = [frames,trials(i).frames(1)+1 : trials(i).frames(2)-1];
    end
    data = single(gpuArray(jksbxreadrand(rfn,nframes,frames, im_plane)));
else
    skip = floor(size(T,1)/nframes);
    data = single(gpuArray(jksbxreadskip(rfn,nframes,skip, im_plane)));
end
data = zscore(data,[],3);

% compute and display correlation map...

% handles.status.String = 'Computing correlation map';
% drawnow;

% corrmap = zeros([size(data,1), size(data,2)],'single','gpuArray');
%     
% for(m=-1:1)
%     for(n=-1:1)
%         if(m~=0 || n~=0)
%             corrmap = corrmap+squeeze(sum(data.*circshift(data,[m n 0]),3));
%         end
%     end
% end
% corrmap = corrmap/8/size(data,3);
% 
% cm = gather(corrmap);
% cm = (cm-min(cm(:)))/(max(cm(:))-min(cm(:)));

% bgimg.CData(:,:,1) = uint8(255*cm);
% bgimg.CData(:,:,2) = uint8(0);
% bgimg.CData(:,:,3) = uint8(0);

% drawnow;

cell_list = [];
% cellpoly = {};
if exist([rfnx,'.segment'],'file')        
    load([rfnx,'.segment'],'-mat')
    if info.volscan
        temp_mask = mask{im_plane};
        ncell = max(temp_mask(:));        
    else
        temp_mask = mask;
        ncell = max(mask(:));
    end
    
    if ncell > 0
        hold(bgimg.Parent,'on');
        for i = 1 : ncell
            temp_ind = ind2sub(size(temp_mask), find(temp_mask == i));            
            if ~isempty(temp_ind)
                cell_list = [cell_list, i];
                temp_mat = zeros(size(temp_mask));
                temp_mat(temp_ind) = 1;
                bw = bwlabel(temp_mat);
                B = bwboundaries(bw);
                xy = B{1};
                patch_h{i} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.7,'edgecolor',[1 1 1],'parent',bgimg.Parent,'FaceLighting','none');
%                 cellpoly{i} = patch_h{i};
                drawnow;
            end
        end
    end
    hold(bgimg.Parent,'off'); 
    
    if ncell > 0
        hold(zimg.Parent,'on');
        for i = 1 : ncell
            temp_ind = ind2sub(size(temp_mask), find(temp_mask == i));            
            if ~isempty(temp_ind)
                cell_list = [cell_list, i];
                temp_mat = zeros(size(temp_mask));
                temp_mat(temp_ind) = 1;
                bw = bwlabel(temp_mat);
                B = bwboundaries(bw);
                xy = B{1};
                patch_zh{i} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.7,'edgecolor',[1 1 1],'parent',zimg.Parent,'FaceLighting','none');                
                drawnow;
            end
        end
    end    
    hold(zimg.Parent,'off');
else
    if info.volscan
        mask = cell(1,plane_num);
        for i = 1 : plane_num
            mask{i} = zeros([size(data,1), size(data,2)]);
        end
    else
        mask = zeros([size(data,1), size(data,2)]);
    end
end
handles.listbox1.String = cell_list;
handles.status.String = 'Showing correlation map. Start segmenting';
m = findobj(segmenttool_h,'tag','method');
if m.Value < 3
    set(segmenttool_h,'WindowButtonMotionFcn',@jksbxwbmcb)
    set(segmenttool_h,'WindowScrollWheelFcn',@jksbxwswcb)
    set(segmenttool_h,'WindowButtonDownFcn',@jksbxwbdcb)
else
    if ~isempty(manual_h)
        delete(manual_h);        
    end
    manual_h = imfreehand;
end
end   

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rfn mask cellpoly pathname cell_list

save([rfn '.segment'],'mask');
handles.status.String = sprintf('Saved %d cells in %s.segment',length(cell_list),rfn);
end

function frames_Callback(hObject, eventdata, handles)
% hObject    handle to frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frames as text
%        str2double(get(hObject,'String')) returns contents of frames as a double

global nframes
nframes = str2double(hObject.String);
end

% --- Executes during object creation, after setting all properties.
function frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function ax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax

global bgimg

bgimg = imagesc(zeros(512,805,3,'uint8'));
colormap gray
axis off
end

function nhood_Callback(hObject, eventdata, handles)
% hObject    handle to nhood (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nhood as text
%        str2double(get(hObject,'String')) returns contents of nhood as a double

global nhood
nhood = str2double(hObject.String);
end

% --- Executes during object creation, after setting all properties.
function nhood_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nhood (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% 'NHOOD'
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global segmenttool_h frames th_corr zs ps
segmenttool_h = hObject;
zs = 0;
ps = 0;

frames = 300;
th_corr = 0.2;
end


% --- Executes during object creation, after setting all properties.
function bgimg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax
end

% --- Executes during object creation, after setting all properties.
function status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global status

status = hObject;
end

% --- Executes on selection change in method.
function method_Callback(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method

global method segmenttool_h
global manual_h

method = hObject;

m = findobj(segmenttool_h,'tag','method');
if m.Value < 3
    handles.addbutton.Visible = 'off';
    set(segmenttool_h,'WindowButtonMotionFcn',@jksbxwbmcb)
    set(segmenttool_h,'WindowScrollWheelFcn',@jksbxwswcb)
    set(segmenttool_h,'WindowButtonDownFcn',@jksbxwbdcb)
else
    handles.addbutton.Visible = 'on';
    manual_h = imfreehand;
end

end

% --- Executes during object creation, after setting all properties.
function method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function nhsize_Callback(hObject, eventdata, handles)
% hObject    handle to nhsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nhsize as text
%        str2double(get(hObject,'String')) returns contents of nhsize as a double

global nhood 
nhood = str2double(hObject.String);
end

% --- Executes during object creation, after setting all properties.
function nhsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nhsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global  nhood_h
nhood_h = hObject;
end

% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


global segmenttool_h

switch get(hObject,'Value')
    case 1
        pan(segmenttool_h,'off');
        zoom(segmenttool_h,'off');
    case 2
        pan(segmenttool_h,'off');
        zoom(segmenttool_h,'on');
    case 3
        pan(segmenttool_h,'on');
        zoom(segmenttool_h,'off');
end
end




% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global mode_h

mode_h = hObject;
end


% --- Executes on button press in extract.
function extract_Callback(hObject, eventdata, handles)
% hObject    handle to extract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rfn sig info im_plane
if info.volscan
    sig = jksbxpullsignals(rfn, im_plane);
else
    sig = jksbxpullsignals(rfn);
end
handles.status.String = sprintf('Signals extracted and saved');
if info.volscan
    plot(handles.axes3,zscore(sig{im_plane}));
else
    plot(handles.axes3,zscore(sig));
end
handles.axes3.Visible = 'off';
handles.axes3.YLim = [-2 10];
end

% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3

axis off;
end

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
end

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

global bgimg data segmenttool_h nframes cell_list cellpoly mask rfn pathname
global plane_num im_plane info patch_h cm

previous_plane = im_plane;
im_plane = get(handles.popupmenu3,'Value');
if previous_plane == im_plane
    return
else
    handles.status.String = 'Resetting/Clearing GPU';
    drawnow;
    gpuDevice(1);

    global sig;
    delete(handles.axes3.Children);
    sig = [];
    rfnx = rfn;

    load('-mat',[rfnx '.align']); 
    axis off

    % global info
    handles.status.String = 'Loading alignment data';

    if(exist('mnr','var'))
        m = gather(mnr);
    else
        m = double(m); % 2016/11/01 JK
    end

    m = (m-min(m(:)))/(max(m(:))-min(m(:)));
    x = adapthisteq(m);
    x = single(x);
    x = (x-min(x(:)))/(max(x(:))-min(x(:)));

    bgimg = image(zeros(size(x,1),size(x,2),3,'uint8'));

    bgimg.CData(:,:,1) = uint8(255*x);
    bgimg.CData(:,:,2) = bgimg.CData(:,:,1);
    bgimg.CData(:,:,3) = bgimg.CData(:,:,1);

    % if(~isempty(cellpoly))
    %     cellfun(@delete,cellpoly);
    % end

    drawnow;

    handles.status.String = 'Loading spatio-temporal data';

    nframes = str2double(handles.frames.String);
    skip = floor(size(T,1)/nframes);
    data = single(gpuArray(jksbxreadskip(rfn,nframes,skip, im_plane))); % jksbxreadskip is modified to have imaging plane as the last input argument
    data = zscore(data,[],3);

    % compute and display correlation map...

    handles.status.String = 'Computing correlation map';
    drawnow;

    corrmap = zeros([size(data,1), size(data,2)],'single','gpuArray');

    for(m=-1:1)
        for(n=-1:1)
            if(m~=0 || n~=0)
                corrmap = corrmap+squeeze(sum(data.*circshift(data,[m n 0]),3));
            end
        end
    end
    corrmap = corrmap/8/size(data,3);

    cm = gather(corrmap);
    cm = (cm-min(cm(:)))/(max(cm(:))-min(cm(:)));

    bgimg.CData(:,:,1) = uint8(255*cm);
    bgimg.CData(:,:,2) = uint8(0);
    bgimg.CData(:,:,3) = uint8(0);

    drawnow;

    set(segmenttool_h,'WindowButtonMotionFcn',@jksbxwbmcb)
    set(segmenttool_h,'WindowScrollWheelFcn',@jksbxwswcb)
    set(segmenttool_h,'WindowButtonDownFcn',@jksbxwbdcb)

    cell_list = [];
    % cellpoly = {};
    if exist([rfnx,'.segment'],'file')        
        load([rfnx,'.segment'],'-mat')
        if info.volscan
            temp_mask = mask{im_plane};
            ncell = max(temp_mask(:));        
        else
            temp_mask = mask;
            ncell = max(mask(:));
        end

        if ncell > 0
            hold(bgimg.Parent,'on');
            for i = 1 : ncell
                temp_ind = ind2sub(size(temp_mask), find(temp_mask == i));            
                if ~isempty(temp_ind)
                    cell_list = [cell_list, i];
                    temp_mat = zeros(size(temp_mask));
                    temp_mat(temp_ind) = 1;
                    bw = bwlabel(temp_mat);
                    B = bwboundaries(bw);
                    xy = B{1};
                    patch_h{i} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.7,'edgecolor',[1 1 1],'parent',bgimg.Parent,'FaceLighting','none');
    %                 cellpoly{i} = patch_h{i};
                    drawnow;
                end
            end
        end
        hold(bgimg.Parent,'off');        
    else
        if info.volscan
            mask = cell(1,plane_num);
            for i = 1 : plane_num
                mask{i} = zeros([size(data,1), size(data,2)]);
            end
        else
            mask = zeros([size(data,1), size(data,2)]);
        end
    end
    handles.listbox1.String = cell_list;
    handles.status.String = 'Showing correlation map. Start segmenting';
end
end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
end

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in addbutton.
function addbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global info cell_list manual_h mask im_plane bgimg zimg patch_h patch_zh nm segmenttool_h
        
    lb = findobj(segmenttool_h,'tag','listbox1');
    lb.Max = 2;
    
    B = bwboundaries(createMask(manual_h,bgimg));
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
    hold(bgimg.Parent,'off');

    %do the same for zimg
    hold(zimg.Parent,'on');
    zimg.CData(:,:,1) = uint8(255*nm);
    patch_zh{ncell} = patch(xy(:,2),xy(:,1),'white','facecolor',[1 .7 .7],'facealpha',0.2,'edgecolor',[1 1 1],'parent',zimg.Parent,'FaceLighting','none');
    hold(zimg.Parent,'off');
    drawnow;

    if info.volscan
        hold(bgimg.Parent,'on');
        [xx,yy] = meshgrid(1:size(mask{im_plane},2),1:size(mask{im_plane},1));
        bw=(inpolygon(xx,yy,xy(:,2),xy(:,1)));
        mask{im_plane}(bw) = ncell;    

        % rescale...
        r = single(bgimg.CData(:,:,1));
        r(mask{im_plane}>0) = 0;
        r = (r-min(r(:)))/(max(r(:))-min(r(:)));
        bgimg.CData(:,:,1) = uint8(r*255);
        status.String = sprintf('Segmented %d cells',length(cell_list));
        hold(bgimg.Parent,'off');
    else
        hold(bgimg.Parent,'on');
        [xx,yy] = meshgrid(1:size(mask,2),1:size(mask,1));
        bw=(inpolygon(xx,yy,xy(:,2),xy(:,1)));
        mask(bw) = ncell;

        r = single(bgimg.CData(:,:,1));
        r(mask>0) = 0;
        r = (r-min(r(:)))/(max(r(:))-min(r(:)));
        bgimg.CData(:,:,1) = uint8(r*255);
        status.String = sprintf('Segmented %d cells',length(cell_list));
        hold(bgimg.Parent,'off');
    end
    
    delete(manual_h)
    manual_h = imfreehand;
    
end
