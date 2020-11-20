function mainExtractFinalStatesFromFolder(folder)

    d = dir(folder + "/*/*.mat");

    for n = 1:length(d)
        filename = strcat(d(n).folder,'/',d(n).name);

        data = extractFinalState(filename);

        save(strcat(d(n).folder,'/',"data-out.mat"),"data")
    end

end
