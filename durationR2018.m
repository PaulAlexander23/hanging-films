function secondsOut = durationR2018(stringIn)
    secondsOut = duration(cell2mat(cellfun(@str2num,split(char(stringIn),':'),'UniformOutput',false)'));
end
