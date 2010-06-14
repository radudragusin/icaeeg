function [Zr,Wr] = plotMCMC4eeg(result,prs)
% Run MCMC4eeg first 
% or
% load MCMC4eeg_result

d = 61;
N = size(result.hZ,1)/prs.m;
Z = reshape( median( result.hZ, 2 ), prs.m, N );
Wr = reshape( median( result.hW, 2 ), d, prs.m);

% Standardizing
Energy = zeros(1,prs.m);
for j=1:prs.m,
    amplZ = std(Z(j,:));                % standard deviation
    signZ = sign(mean(Z(j,:)));         % scale with the sign of Z
    Z(j,:) = signZ*Z(j,:)/amplZ;
    Wr(:,j) = signZ*Wr(:,j)*amplZ;
    Energy(j) = sum(Wr(:,j).*Wr(:,j));  % energy of W
end

% Sorting (arrange components by energy)
[Energy,idx] = sort(Energy,'descend');
SortZ = zeros(size(Z));
SortW = zeros(size(Wr));
for j = 1:prs.m
    SortZ(j,:) = Z(idx(j),:);
    SortW(:,j) = Wr(:,idx(j));
end
Zr = SortZ; % components activation (time course)
Wr = SortW; % components scalp map
WZ = Wr*Zr; % Project the components back

%% Plot of the waves per component

% figure;
% for i=1:prs.m
%     subplot(prs.m,1,i);
%     plot(1:256, WZ(i,:), 'b' );
%     set(gca,'xticklabel',[]);
%     axis tight;
% end

%% Topographic maps

% Prepare for plots
ti=-1.1:0.01:1.1;
[qx,qy] = meshgrid(ti,ti);
load electrodes % C = 61 x 5 cell
Xs = cell2mat(C(:,2)); Ys = cell2mat(C(:,3));

% Plot topographic map per component
figure;
for i=1:prs.m
    subplot(8,8,floor(i/8)*8+mod(i,8));
    Zs = Wr(:,i);
    F = TriScatteredInterp(Xs,Ys,Zs);
    qz = F(qx,qy);
    imagesc(qz);
    title(strcat('Component',int2str(i)));
    axis off;
end

% Create movie of topographic maps in time, per component
% figure
% for t=1:256
%     clf
%     for i=1:1 %prs.m
%         %subplot(1,4,i);
%         % Compute the projection of the ith component onto the 
%         % original data channels:
%         Zs = Wr(:,i)*Zr(i,t); 
%         F = TriScatteredInterp(Xs,Ys,Zs);
%         qz = F(qx,qy);
%         imagesc(qz);
%         title(strcat('Component',int2str(i)));
%         axis off;
%     end
%     M(t) = getframe;
% end

end