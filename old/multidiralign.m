list = [436,437,438,450,452,467,468];

for i = 1 : length(list)
    cd(num2str(list(i)));
    sbxaligndir
    cd('..')
end
