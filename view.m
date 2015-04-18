load -ascii fullData.txt
<<<<<<< HEAD
iter=400;
N = 4096;
=======
iter=25;
N = 64;
>>>>>>> 367de5a92837f93d2eba5fa8a3eb98117f023247
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
<<<<<<< HEAD
    %subplot(5,5,i);
    %imagesc(data{i});
    %colormap('gray');
    %axis equal;
    %xlim([0.5,8.5]);
    %ylim([0.5,8.5])
    
end
nframe=iter;
%x=rand(100,nframe);
%mov(1:nframe)= struct('cdata',[],'colormap',[]);
%set(gca,'nextplot','replacechildren')
for k=1:nframe
    imagesc(data{k});
    pause(0.01);
%    mov(k)=getframe(gcf);
end
%movie2avi(mov, '1moviename.avi', 'compression', 'None');
=======
    subplot(5,5,i);
    imagesc(data{i});
    colormap('gray');
    axis equal;
    xlim([0.5,8.5]);
    ylim([0.5,8.5])
    
end
>>>>>>> 367de5a92837f93d2eba5fa8a3eb98117f023247
