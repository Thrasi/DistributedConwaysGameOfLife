[R, it, w, h, d] = load('metainfodata.txt');

data = cell(R, it);
for i = 1:R
    file = fopen([num2str(i-1) 'data.txt'], 'r');
    f_data = textscan(file, '%s', iter, 'Delimiter', '\n');
    f_data = f_data{1};
    for j = 1:it
        data{i,j} = str2num(f_data{j});
    end
end

%{
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
%}