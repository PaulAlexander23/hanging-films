%DATA_SAVE

filename = replace(sprintf('data-theta-%g-Re-%g-C-%g-xL-%g-yL-%g-T-%g-interface-%s-xN-%g-yN-%g-AbsTol-%g', ...
    theta,Re,C,x_length,y_length,t_final,func2str(interface),xN(1),xN(2),AbsTol),'.','_');
save(filename,'y','params','t','x','timeTaken');