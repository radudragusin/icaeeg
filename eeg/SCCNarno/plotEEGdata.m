function [] = plotEEGdata(EEG)
%
% Compatible with the EEG data from
% http://sccn.ucsd.edu/~arno/fam2data/publicly_available_EEG_data.html
%
% Use savefromCNTdata first on the .CNT and .DAT files from the dataset. 
% This will create a .MAT file that can be used by this function.
%
% Creates a movie of the data contained in the 31 channels of the EEG
% structure given as argument.
%
% >>> Important! The colorbar changes at each frame.

% Load the file containing channel locations and the mask of the brain
load coords31;
cx = chloc.X; cy=chloc.Y;
chx = [cx ; mask.X];
chy = [cy ; mask.Y];

% Prepare the image 
x = 1:392;
y = 1:397;
[xx,yy] = meshgrid(x,y);
clear x y;
map=colormap(jet(256));
map(1,:)=[0 0 0];

% Create the movie frame by frame
for s=1:size(EEG.data,2),
    mag = [double(EEG.data(:,s)) ; zeros(1089,1)];

    % Set image values equal to input values
    I = zeros(397,392);
    for i=1:length(chloc.X),
        I(chloc.X(i),chloc.Y(i)) = mag(i);
    end

    z = griddata(chy,chx,mag,xx,yy,'cubic');
    clear mag;

    % Adjust colorbar to include black background
    ma = max(z(:));
    mi = min(z(:));
    step = (ma-mi)/(256-2);
    black = mi + (-1)*step;
    z(isnan(z)) = black;

    % Set electrode sites equal to black
    for i=1:length(cx),
        z(cx(i),cy(i))=black;
        z(cx(i)+1,cy(i))=black;
        z(cx(i)-1,cy(i))=black;
        z(cx(i)-1,cy(i)-1)=black;
        z(cx(i),cy(i)-1)=black;
        z(cx(i)+1,cy(i)-1)=black;
        z(cx(i)-1,cy(i)+1)=black;
        z(cx(i),cy(i)+1)=black;
        z(cx(i)+1,cy(i)+1)=black;
    end

    % Plot the brain 
    clf;
    imagesc(z);
    colormap(map);
    colorbar;
    M(s) = getframe;
end
clear all;
end
