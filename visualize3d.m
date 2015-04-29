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
    for i = 1:N
        id = ids(i);
        zs(i) = floor(id / (w*h));
        id = mod(id, w*h);
        ys(i) = floor(id / w);
        xs(i) = mod(id, w);
    end

    scatter3(xs, ys, zs, 225*ones([N, 1]), 's', 'filled');
    axis([0 w 0 h 0 d]);
    %axis image
    pause(0.3);
    %figure;
end