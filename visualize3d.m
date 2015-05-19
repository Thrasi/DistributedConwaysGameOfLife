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
kit = 18:23
for j = 1:it
    clf;
    %subplot(3,2,j);
    title(sprintf('Iteration %d',j));
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    for k = 1:R
        m = data{k,j};
        for i = 1:size(m, 2)
            voxel(m(:,i), [1,1,1], 'r', 0.7);
        end
    end
    voxel([0 0 0], [5,5,5], 'g', 0.1);
    voxel([0 0 5], [5,5,5], 'b', 0.1);
    voxel([0 5 0], [5,5,5], 'y', 0.1);
    voxel([0 5 5], [5,5,5], 'c', 0.1);
    voxel([5 0 0], [5,5,5], 'y', 0.1);
    voxel([5 0 5], [5,5,5], 'c', 0.1);
    voxel([5 5 0], [5,5,5], 'g', 0.1);
    voxel([5 5 5], [5,5,5], 'b', 0.1);
    grid on
    axis([0 w 0 h 0 d]);
    pause(0.2);
end
