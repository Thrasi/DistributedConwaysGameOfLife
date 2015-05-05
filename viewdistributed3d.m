% DONT USE THIS

load 0data.txt

a = reshape(X0data(1,:),[5 5 5]);
ind = find(a);
[x,y,z] = ind2sub(size(a),ind);
l=[x,y,z];
figure();
hold on;
for i=1:4
    voxel(l(i,:),[1 1 1],'r',0.7);
end
axis equal;
axis([0 10 0 10 0 10]);