function z = jksbxreadframes(fname, frames, channels)

% read desired frames (in a vector form) from sbx file

% z = sbxread(fname,0,1);
% z = zeros([size(z,2) size(z,3) length(frames)]);
% 
% for i = 1:length(frames)
%     temp = sbxread(fn,frames(i),1);
%     z(:,:,i) = squeeze(temp(1,:,:));
% end


% modified sbxread to select multiple non-consecutive frames 2020/12/03
% img = sbxread(fname,k,N,varargin)
%
% Reads from frame k to k+N-1 in file fname
% 
% fname - the file name (e.g., 'xx0_000_001')
% k     - the index of the first frame to be read.  The first index is 0.
% N     - the number of consecutive frames to read starting with k.
%
% If N>1 it returns a 4D array of size = [#pmt rows cols N] 
% If N=1 it returns a 3D array of size = [#pmt rows cols]
%
% #pmts is the number of pmt channels being sampled (1 or 2)
% rows is the number of lines in the image
% cols is the number of pixels in each line
%
% The function also creates a global 'info' variable with additional
% informationi about the file

global info_loaded info

% check if already loaded...

if(isempty(info_loaded) || ~strcmp(fname,info_loaded)) % because of
% frequent error in info variable 2017/07/14 JK
    
%     if(~isempty(info_loaded))   % try closing previous...
%         try
%             fclose(info.fid);
%         catch
%             error('could not close info.fid')
%         end
%     end

    load(fname);
    
    if(exist([fname ,'.align'])) % aligned?
        info.aligned = load([fname ,'.align'],'-mat');
    else
        info.aligned = [];
    end   
    
    info_loaded = fname;
    
    if(~isfield(info,'sz'))
        info.sz = [512 796];    % it was only sz = .... 
    end
    
    if(~isfield(info,'scanmode'))
        info.scanmode = 1;      % unidirectional
    end
    
    if(info.scanmode==0)
        info.recordsPerBuffer = info.recordsPerBuffer*2;
    end
    
    switch info.channels
        case 1
            info.nchan = 2;      % both PMT0 & 1
            factor = 1;
        case 2
            info.nchan = 1;      % PMT 0
            factor = 2;
        case 3
            info.nchan = 1;      % PMT 1
            factor = 2;
    end
    
    info.fid = fopen([fname '.sbx']);
    d = dir([fname '.sbx']);
    info.nsamples = (info.sz(2) * info.recordsPerBuffer * 2 * info.nchan);   % bytes per record 
    %Edit Patrick: to maintain compatibility with new version
    
    if isfield(info,'scanbox_version') && info.scanbox_version >= 2
        info.max_idx =  d.bytes/info.recordsPerBuffer/info.sz(2)*factor/4 - 1;
        info.nsamples = (info.sz(2) * info.recordsPerBuffer * 2 * info.nchan);   % bytes per record 
    else
        info.max_idx =  d.bytes/info.bytesPerBuffer*factor - 1;
    end
end

if(isfield(info,'fid') && info.fid ~= -1)
    if channels == 1 || channels == 2 || channels == [1,2]
    %     nsamples = info.postTriggerSamples * info.recordsPerBuffer;

        N = 1;
%         try
%             % initialize
%             z = zeros([length(channels), info.sz(1), info.sz(2), length(frames)],'uint16');
%             % loop through frames
%             for fi = 1 : length(frames)
%                 fseek(info.fid,frames(fi)*info.nsamples,'bof');
%                 x = fread(info.fid,info.nsamples/2 * N,'uint16=>uint16');
%                 x = reshape(x,[info.nchan info.sz(2) info.recordsPerBuffer  N]);
%                 x = intmax('uint16')-permute(x,[1 3 2 4]);
%                 x = x(channels, :,:,:);
%                 z(:,:,:,fi) = x;
%             end
%         catch
            try         
%                 info.fid = fopen([fname '.sbx']); % for some reason can't be specified, fseek does not work, and it could be recovered by opening
                % initialize
                z = zeros([length(channels), info.sz(1), info.sz(2), length(frames)],'uint16');
                % loop through frames
                for fi = 1 : length(frames)
                    info.fid = fopen([fname '.sbx']);
                    fseek(info.fid,frames(fi)*info.nsamples,'bof');
                    x = fread(info.fid,info.nsamples/2 * N,'uint16=>uint16');
                    x = reshape(x,[info.nchan info.sz(2) info.recordsPerBuffer  N]);
                    x = intmax('uint16')-permute(x,[1 3 2 4]);
                    x = x(channels, :,:,:);
                    z(:,:,:,fi) = x;
                    fclose(info.fid); % test remedy 2017/07/14 JK
                end
            catch
                error('Cannot read frame.  Index range likely outside of bounds.');
            end
%         end

    %     x = intmax('uint16')-permute(x,[1 3 2 4]);

        
    else
        error('Input ''channels'' should be either [1], [2], or [1,2]')
    end
    
else
    z = [];
end