metadata = num2cell(load('metainfodata.txt'));
[R, it, w, h, d] = deal(metadata{:});

data = cell(R, it);
for i = 1:R
    file = fopen([num2str(i-1) 'data.txt'], 'r');
    f_data = textscan(file, '%s', it, 'Delimiter', '\n');
    f_data = f_data{1};
    for j = 1:it
        row = str2num(f_data{j});
        l = length(row);
        p = l/3;
        data{i,j} = reshape(row, [3, p]);
    end
    fclose(file);
end

for j = 1:it
    clf;
    for k = 1:R
        m = data{k,j};
        for i = 1:size(m, 2)
            voxel(m(:,i), [1,1,1], 'r', 0.7);
        end
    end
    axis([0 w 0 h 0 d]);
    pause(0.3);
end
