%function [ ] = plotEEGdata( )
%
% Compatible with the EEG data from
% http://archive.ics.uci.edu/ml/datasets/EEG+Database
%
% Reads the text file containing the EEG data for one second, plots the
% waves for the 64 electrods, and makes a movie with the values of the 
% electrodes for each frame of that one second.
%
% Uses coords.mat

% Read the file which returns: the trial number, sensor position (channel),
% sample number (0-255), and sensor value (in micro volts).
[~,C,~,V] = textread('co2c0000337.rd.000','%d %s %d %f','headerlines',5,'commentstyle','shell');

% Reshape the data for further usage
V = reshape(V,256,64);
V = V';
[~, I] = unique(C,'first'); 
C = C(sort(I));
clear I;

% % Plot the waves for each electrod
% figure;
% subplot(size(V,1),1,1);
% plot(1:size(V,2),V(1,:));
% ylabel(C(1));
% set(gca,'xaxislocation','top');
% axis tight;
% for i=2:size(V,1),
%     subplot(size(V,1),1,i);
%     plot(1:size(V,2),V(i,:));
%     ylabel(C(i));
%     set(gca,'xticklabel',[]);
%     axis tight;
% end
% set(gca,'xtickMode', 'auto');

% Load the file containing the XYZ location of each of the 61 channels and
% the name of those channels
load coords

Xs = Coordinates(:,1);
Ys = Coordinates(:,2);
Zs = zeros(length(Channels),1);
ti=-1.1:0.01:1.1;
scrsz = get(0,'ScreenSize');

channels_text = cell(length(Channels),1);
indxs = zeros(length(Channels),1);
for i=1:length(Channels),
    for j=1:length(C),
        if all(strcmp(C{j}, Channels{i})) 
            indxs(i) = j; 
            break
        end
    end
    channels_text{i} = char(Channels(i));
end

figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
[qx,qy] = meshgrid(ti,ti);

values_text = cell(length(Channels),1);
for s=1:size(V,2),
    for i=1:length(Channels),
        value = V(indxs(i),s);
        Zs(i) = value;
        values_text{i} = mat2str(value);
    end
    F = TriScatteredInterp(Xs,Ys,Zs);
    qz = F(qx,qy);
    clf;
    imagesc(qz,[min(min(V)) max(max(V))]);
    hold on;
    scatter((Xs+1.1)*100,(Ys+1.1)*100,20,'o','filled','white');
    text((Xs+1.1)*100,(Ys+1.1)*100,channels_text);
    text((Xs+1.1)*100,(Ys+1.1)*100+4,values_text);
    hold off;
    M(s) = getframe;
end
