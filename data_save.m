%DATA_SAVE

filename = replace(sprintf('data-theta-%f-Re-%f-We-%f-C-%f-xL-%f-yL-%f-T-%f',[params(2:end),xL,tL]),'.','_');
save(filename,'x','t','y','params');