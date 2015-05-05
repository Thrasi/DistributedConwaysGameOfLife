%data = load('testvis');
file = fopen('testvis', 'r');
header = textscan(file, '%d %d %d %d\n', 1);
it = header{1};
w = header{2};
h = header{3};
d = header{4};
data = textscan(file, '%s', iter, 'Delimiter', '\n');
data = data{1};

for j = 1:it
    row = str2num(data{j});
    l = length(row);
    p = l/3;
    m = reshape(row, [3, p]);
    
    clf;
    for i = 1:p
        voxel(m(:,i), [1,1,1], 'r', 0.7);
    end
    axis([0 w 0 h 0 d]);
    pause(0.3);
end