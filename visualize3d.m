metadata = num2cell(load('data/metainfodata.txt'));
[R, it, w, h, d] = deal(metadata{:});

data = cell(R, it);
for i = 1:R
    file = fopen(['data/' num2str(i-1) 'data.txt'], 'r');
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
    title(sprintf('Iteration %d',j));
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    for k = 1:R
        m = data{k,j};
        for i = 1:size(m, 2)
            voxel(m(:,i), [1,1,1], 'r', 0.5);
        end
    end
    grid on
    axis([0 w 0 h 0 d]);
    pause(0.5);
end
