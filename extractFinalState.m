function data = extractFinalState(filename)

    data = load(filename);

    data.solution.y = data.solution.y(:,:,end);
end
