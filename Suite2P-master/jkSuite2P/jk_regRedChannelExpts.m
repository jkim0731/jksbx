function mimgR = jk_regRedChannelExpts(ops)


ops.rchannel = getOr(ops, 'rchannel', 2);
% numPlanes = length(ops.planesToProcess);

%% build file list with red channel
flag = 0;
if (isfield(ops, 'expred') && ~isempty(ops.expred))
    if sum(ismember(ops.session, ops.expred)) > 0
        disp('red/green channel experiments found');
    else
        flag = 1;
    end
else
        flag = 1;
end
if flag
    error('no red/green channel experiment supplied')
end

% green file list
fs        = ops.fsroot;

%% was GREEN channel registered?
DS = cell(ops.nplanes, 1);
try
    root = ops.ResultsSavePath;
    fname = sprintf('regops_%s_%03d.mat', ops.mouse_name, ops.session);
    load(fullfile(root, fname))
    for j = 1:ops.nplanes
        DS{j} = ops1{j}.DS;
        ops.BiDiPhase = ops1{j}.BiDiPhase;
    end
catch       
    error('registration not run for green channel');
end


%% loop over experiments -- register ones with red channel (expred)
ired = 0;
ntf0 = 0;
for k = 1:length(fs)
    % is this an experiment with a red channel?
    if ismember(ops.session(k), ops.expred)
        ired = ired + 1;
        numPlanes = ops.nplanes;
        iplane0 = 1:1:ops.nplanes;
        clear DSexp;
        for iplane = 1:ops.nplanes
            csumNframes = [0 cumsum(ops1{iplane}.Nframes)];
            DSexp{iplane} = DS{iplane}(csumNframes(k)+1:csumNframes(k+1),:);
        end
        
        ix0 = zeros(ops.nplanes,1);
        nbytes = 0;
        for j = 1:length(fs{k})
            if abs(nbytes - fs{k}(j).bytes)>1e3
                nbytes = fs{k}(j).bytes;
                nFr = nFramesTiff(fs{k}(j).name);
            end
            iplane0 = mod(iplane0-1, ops.nplanes) + 1;
    
            ichanset = [ops.rchannel; nFr; ops.nchannels_red];
            % only load frames of RED channel
            data = jkloadFramesBuff(fs{k}(j).name, ops., ops.rchannel);
            
            [Ly, Lx, ~] = size(data);
            if abs(ops.BiDiPhase) > 0
                data = ShiftBiDi(ops.BiDiPhase, data, Ly, Lx);
            end
            
            if ired == 1
                mimgR = zeros(Ly, Lx, ops.nplanes);
            end
            
            for iplane = 1:length(ops.planesToProcess)
                ifr0 = iplane0(ops.planesToProcess(iplane));
                indframes = ifr0:ops.nplanes:size(data,3);
                dataR0    = data(:,:,indframes);
                nt        = size(dataR0,3);
                pframes   = ops1{iplane}.Nframes;
                
                ds        = DSexp{iplane}(ix0(iplane)+[1:nt], :);
                
                dreg      = register_movie(dataR0, ops, ds);
                ix0(iplane) = ix0(iplane) + nt;
                
                mimgR(:,:,iplane) = mimgR(:,:,iplane) + mean(dreg, 3);            
            end
            
            iplane0 = iplane0 - nFr/ops.nchannels;
             
            ntf0 = ntf0 + 1;
%             disp(ntf0);
        end
    end
end

mimgR = mimgR/ntf0;