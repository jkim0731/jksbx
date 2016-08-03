N = length(info.event_id);
temp = 0;
prevframe = 0;
prevline = 0;

dframe = [];
dline = [];
for i = 1 : N
    switch temp
        case 0
            if info.event_id(i) == 1
                temp = 1;
                prevframe = info.frame(i);
                prevline = info.line(i);
            end
        case 1
            if info.event_id(i) == 1
                dframe(end+1) = info.frame(i) - prevframe;
                dline(end+1) = info.line(i) - prevline;
                prevframe = info.frame(i);
                prevline = info.line(i);
            elseif info.event_id(i)== 2
                temp = 0;
            end
    end
end
        
for i = 1 : length(dframe)
    if dframe(i) == 1
        dline(i) = dline(i) + 780;
    end
end

dline_test = round(dline/114); % 114 for 30 Hz. maybe 57 for 15 Hz. 