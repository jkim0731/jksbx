function jkdir_func(func_str)
% make every kind of function working for all files in a dir. skip methods should have been declared in each function. 

d = dir('*.sbx');
for(i=1:length(d))
%     try
        fn = strtok(d(i).name,'.');
        eval([func_str,'(''',fn,''')'])
%     catch
%         sprintf('Could not perform %s on %s',func_str, fn)
%     end
end
end
            