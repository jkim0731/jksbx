function ops = jk_build_ops3(db, ops)


ops.nplanes = getOr(ops, 'nplanes', 1);
ops.nchannels = getOr(ops, 'nchannels', 1);
ops.readTiffHeader = getOr(ops,'readTiffHeader',1);


% ops = db;
if ~iscell(db.fn)
    % this is the usual case where we have a simple single session recording
    ops = addfields(ops, db);
    
    if ~isfield(ops, 'RootDir')
        ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name);
    end
    
else
    % here we might have multiple sessions, which we want to be analyzed
    % together (exactly the same FOV)
    nSessions = length(db.fn);
    % a backwards compatible version of db
    dbCompat = db;
    dbCompat.fn = db.fn{1};
    ops = addfields(ops, dbCompat);
    ops.db_orig = db;
    
    % this line to be backward compatible (just in case)
    ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name);
end

if ~(isfield(ops, 'planesToProcess') && ~isempty(ops.planesToProcess))
    ops.planesToProcess = 1:ops.nplanes;
% else
%     % planesToProcess is not working right now
%     ops.planesToProcess = 1:ops.nplanes;
end

curr_dir = pwd;
cd(ops.RootDir);
global info
sbxread(db.fn,0,1);
ops.info = info;
clear info
fclose all;
cd(curr_dir)