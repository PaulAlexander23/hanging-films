
files = dir('*.m');
testFiles = extractTestFiles(files);
testResults = [];

for n = 1:length(testFiles)
    testResults = [testResults, runtests(testFiles(n).name)];
end

if any([testResults.Failed])
    display(table(testResults));
else
    display(testResults);
end

function testFiles = extractTestFiles(files)
    testFiles = [];
    for n = 1:length(files)
        if strncmpi(files(n).name, "test", 4) && files(n).name ~= "test.m"
            testFiles = [testFiles, files(n)];
        end
    end
end
