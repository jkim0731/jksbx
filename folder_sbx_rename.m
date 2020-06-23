% renaming multiple files (.sbx and derivatives)
% using FileRename

flist = dir('052_7002_*.*');
newMid = '9999';
lastIncrement = 200;
for fi = 1 : length(flist)
    fn = flist(fi).name;
    fnstr = strsplit(fn,'_');
    [fnstr{3},fext] = strtok(fnstr{3},'.');
    
    newfnstr = cell(1,3);
    newfnstr{1} = fnstr{1};
    newfnstr{2} = newMid;
    if str2double(fnstr{3}) < 1000
        newfnstr{3} = num2str(str2double(fnstr{3})+lastIncrement);
    else
        newfnstr{3} = num2str(str2double(fnstr{3})+lastIncrement*10);
    end
    newfn = [newfnstr{1}, '_', newfnstr{2}, '_', newfnstr{3}, fext];

    FileRename(fn, newfn)
end