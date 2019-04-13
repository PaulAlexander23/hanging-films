

files = dir("*.mat");
number_of_files = length(files);

Re = zeros(number_of_files,1);
initial_min_h = zeros(number_of_files,1);
final_min_h = zeros(number_of_files,1);

initial_max_h = zeros(number_of_files,1);
final_max_h = zeros(number_of_files,1);

for file_number = 1:number_of_files
    load(files(file_number).name);
    
    Re(file_number) = params(3);
    initial_min_h(file_number) = min(y(:,:,1),[],'all');
    final_min_h(file_number) = min(y(:,:,end),[],'all');
    
    initial_max_h(file_number) = max(y(:,:,1),[],'all');
    final_max_h(file_number) = max(y(:,:,end),[],'all');
end

plot([0,1.5],[1,1]);
hold on
plot(Re,initial_min_h);
plot(Re,final_min_h);
plot(Re,initial_max_h);
plot(Re,final_max_h);