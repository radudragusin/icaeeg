% animMCMC4fmri

% Loads dataset
load data2
% Loads data from MCMC4fmri
load resultMCMC_m4_363_300burn100

d = 3891;
N = size(res.hZ,1)/prs.m;
Z = reshape( median( res.hZ, 2 ), prs.m, N );
Wr = reshape( median( res.hW, 2 ), d, prs.m);

% standardizing
Energy = zeros(1,prs.m);
for j=1:prs.m,
    amplZ = std(Z(j,:));                % standard deviation
    signZ = sign(mean(Z(j,:)));         % scale with the sign of Z
    Z(j,:) = signZ*Z(j,:)/amplZ;
    Wr(:,j) = signZ*Wr(:,j)*amplZ;
    Energy(j) = sum(Wr(:,j).*Wr(:,j));  % energy of W
end

% sorting (arrange components by energy)
[Energy,idx] = sort(Energy,'descend');
SortZ = zeros(size(Z));
SortW = zeros(size(Wr));
for j = 1:prs.m
    SortZ(j,:) = Z(idx(j),:);
    SortW(:,j) = Wr(:,idx(j));
end
Zr = SortZ;
Wr = SortW;

% Time series
fs = 1;                         % Sampling frequency
t = 0 : 1/fs : 363/fs-1/fs;     % Time series

% Reference function
pa=zeros(N,1);
for j=1:3,
    pa((30+(j-1)*121):(60+(j-1)*121))=1;
end

scrsz = get(0,'ScreenSize');
fig = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);

for i=1:prs.m
    sp(i) = subplot( 6, 6, [6*(i+1)+3 6*(i+1)+4 6*(i+1)+5 6*(i+1)+6]);
    plot(t, Zr(i,:), 'b' );     % Time series
    hold on;
    plot(t, pa,'k');
    ylabel(['Comp. ' num2str(i)]);
    axis tight;
    set(gca,'xtick',[],'yaxislocation','right');
    hold off;
end
set(gca,'xtickMode', 'auto');

a1 = get(sp(1),'Position');
a2 = get(sp(2),'Position');
a3 = get(sp(3),'Position');
a4 = get(sp(4),'Position');
annot1 = annotation('line',[0.6 0.6],[0.1 0.8],'Color','red');
annot2 = annotation('line',[0.6 0.6],[0.1 0.8],'Color','red');
annot3 = annotation('line',[0.6 0.6],[0.1 0.8],'Color','red');
annot4 = annotation('line',[0.6 0.6],[0.1 0.8],'Color','red');

for t = 1:size(X1,2)
    imgslice = X1(:,t);
    m = min(imgslice)*ones(82*68,1);
    m(logical(mask1)) = imgslice;
    subplot(6,6,[13 14 19 20 25 26 31 32]);
    imagesc(reshape(m,[82 68]));
    title(strcat('Frame: ',int2str(t),'/',num2str(size(X1,2)),'   Time: ',int2str(floor(t/3)),'/121s'),'FontSize',16);
    axis off;
    
    m1=min(Wr(:,1)*min(Zr(1,:)))*ones(82*68,1);
    m1(logical(mask1))=Wr(:,1)*Zr(1,t);
    m11=reshape(m1,[82 68]);
    subplot(6,6,[3 9]);
    imagesc(m11);
    title('Component 1');
    axis off;
    
    m1=min(Wr(:,2)*min(Zr(2,:)))*ones(82*68,1);
    m1(logical(mask1))=Wr(:,2)*Zr(2,t);
    m2=reshape(m1,[82 68]);
    subplot(6,6,[4 10]);
    imagesc(m2);
    title('Component 2');
    axis off;
   
    m1=min(Wr(:,3)*min(Zr(3,:)))*ones(82*68,1);
    m1(logical(mask1))=Wr(:,3)*Zr(3,t);
    m3=reshape(m1,[82 68]);
    subplot(6,6,[5 11]);
    imagesc(m3);
    title('Component 3');
    axis off;
    
    m1=min(Wr(:,4)*min(Zr(4,:)))*ones(82*68,1);
    m1(logical(mask1))=Wr(:,4)*Zr(4,t);
    m4=reshape(m1,[82 68]);
    subplot(6,6,[6 12]);
    imagesc(m4);
    title('Component 4');
    axis off;
    
    set(annot1, 'Position', [a1(1)+t/363*a1(3), a1(2), 0, a1(4)]);
    set(annot2, 'Position', [a2(1)+t/363*a2(3), a2(2), 0, a2(4)]);
    set(annot3, 'Position', [a3(1)+t/363*a3(3), a3(2), 0, a3(4)]);
    set(annot4, 'Position', [a4(1)+t/363*a4(3), a4(2), 0, a4(4)]);
    
    M(t) = getframe(gcf);
end