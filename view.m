load -ascii fullData.txt
iter=25;
N = 64;
data = cell(1,iter);
a = [];
P = size(fullData,1);
sqrtP = sqrt(P);
n = sqrt(N);
NP = N/P;
pp = sqrtP;

for i=1:iter
    data{i} = zeros(n);
    for p=0:P-1
        tmp = reshape(fullData(p+1, NP*(i-1)+1:NP*i),[sqrt(NP),sqrt(NP)])';
        col = rem(p,pp)*sqrt(NP)+1;
        row = fix(p/pp)*sqrt(NP)+1;
        data{i}(row:row+sqrt(NP)-1,col:col+sqrt(NP)-1) = tmp;
    end
    subplot(5,5,i);
    imagesc(data{i});
    colormap('gray');
    axis equal;
    xlim([0.5,8.5]);
    ylim([0.5,8.5])
    
end