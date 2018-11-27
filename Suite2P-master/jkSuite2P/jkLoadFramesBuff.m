function data = jkLoadFramesBuff(ops, plane, channel)

curr_dir = pwd;
cd(ops.RootDir)

data = jksbxreadframes_multifile(ops.sbxfnlist, ops.frame_to_use{plane}, ops.max_idx);
data = squeeze(data(channel,:,:,:));

cd(curr_dir)