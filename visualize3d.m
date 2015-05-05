data = load('testvis');
w = 5;
h = 5;
d = 5;
it = 30;

for j = 1:it
    offset = (j-1) * h * d;
    data3d = zeros([w, h, d]);
    for i = 1:d
        ids = offset + (i-1)*h+1 : offset + i*h;
        data3d(:,:,i) = data(ids,:);
    end

    ids = find(data3d);
    N = size(ids, 1);
    xs = zeros(N, 1);
    ys = zeros(N, 1);
    zs = zeros(N, 1);
    clf;
    for i = 1:N
        id = ids(i);
        zs(i) = floor(id / (w*h));
        id = mod(id, w*h);
        ys(i) = floor(id / w);
        xs(i) = mod(id, w);
        voxel([xs(i), ys(i), zs(i)], [1,1,1], 'r', 0.7);
    end

    %scatter3(xs, ys, zs, 225*ones([N, 1]), 's', 'filled');
    axis([0 w 0 h 0 d]);
    %axis equal
    pause(0.3);
    %figure;
end