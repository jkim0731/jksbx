%%
%% Manual info correction
%%
fn = '039_9999_208';
load([fn, '.mat'], 'info')
tempMessages = cell(10,1);
tempMessages(2:end) = info.messages(1:9);
tempMessages{1} = '152';
info.messages = tempMessages;
%%
save([fn,'.mat'], 'info')