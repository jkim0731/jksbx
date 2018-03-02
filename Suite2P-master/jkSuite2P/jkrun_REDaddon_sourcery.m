function jkrun_REDaddon_sourcery(db, ops0)
ops = jk_build_ops3(db, ops0);
% red channel addon to already processed data
mimgG = [];

if sum(ismember(ops.session, ops.expred))>0
    mimgR = jk_regRedChannelExpts(ops);
else
    if ops.nchannels_red==2
        % return green mean image from short red/green recording
        [mimgR, mimgG] = jk_regRedGreenChannel(ops);
    else
        mimgR = jk_regRedChannelOnly(ops);
    end
end

save(fullfile(ops.ResultsSavePath,'redchannel.mat'), 'mimgR','mimgG');

jk_add_red_channel_sourcery(ops);
