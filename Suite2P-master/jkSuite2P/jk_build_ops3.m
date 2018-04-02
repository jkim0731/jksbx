function ops = jk_build_ops3(db, ops)


ops.nplanes = getOr(ops, 'nplanes', 1);
ops.nchannels = getOr(ops, 'nchannels', 1);
if ~isfield(ops, 'RootDir')
    ops.RootDir = fullfile(db.RootStorage, db.mouse_name);
end
ops = addfields(ops, db);
    
% if length(ops.session) == 1
    % this is the usual case where we have a simple single session recording
    ops.fsroot{1} = dir(fullfile(ops.RootDir, sprintf('%s_%03d_*.sbx',ops.mouse_name,ops.session))); % leave this as a cell to comply with other codes
        
% else
%     % here we might have multiple sessions, which we want to be analyzed
%     % together (exactly the same FOV)
%     nSessions = length(ops.session);
% %     % a backwards compatible version of db
% %     dbCompat = db;
% %     dbCompat.fn = db.fn{1};
% %     ops = addfields(ops, dbCompat);
% %     ops.db_orig = db;
%     ops.fsroot = cell(0);
%     for iSession = 1:nSessions        
%         ops.fsroot{end+1} = dir(fullfile(ops.RootDir, sprintf('%s_%03d_*.sbx',ops.mouse_name,ops.session(iSession))));
%     end 
%     % this line to be backward compatible (just in case)
% %     ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name);
% end

if ~(isfield(ops, 'planesToProcess') && ~isempty(ops.planesToProcess))
    ops.planesToProcess = 1:ops.nplanes;
end

ops.ResultsSavePath = sprintf('%s//%s//', ops.ResultsSavePath, ops.mouse_name);