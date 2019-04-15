


estimate_file = "data-theta-2_74889-Re-1-We-0-C-0_01-xL-32-yL-32-T-300-xN-256-yN-256-AbsTol-1e-06.mat";

load(estimate_file);

te = t;
xe = x;
ye = y;
[XE,YE,TE] = meshgrid(xe{1},xe{2},te);

files = dir("*.mat");
number_of_files = length(files);

res = zeros(number_of_files,1);
error = zeros(number_of_files,1);

for file_number = 1:number_of_files
    load(files(file_number).name);
    
    res(file_number) = length(x{1});
    
    [X,Y,T] = meshgrid(x{1},x{2},t);
    
    temp = interp3(XE,YE,TE,ye,X,Y,T);
    
    figure
    plot(t,squeeze(max(abs(y-temp),[],[1,2])))
    
    error(file_number) = max(abs(y(:,:,end)-temp),[],'all');
end

[res,I] = sort(res);
plot(res,error(I));