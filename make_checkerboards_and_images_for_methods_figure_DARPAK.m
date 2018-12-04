
% Make checkerboard stimul for figure 1 of Bennett et al. paper
% R Maloney, Nov 2018

%%
sz1 = 560; %the original size in your stimulus code, 48*48 pixels, for example
sz2 = 1024; %the resolution you want the image to be, 1024*1024 pixels, for example

%rescaleMe = sz2/sz1; %now work out the rato between those 2 resolutions

NumCheckSquares = 8; % The total number of squares in a row of the checkerboard (both colours)

% number of squares in row of the checkerboard.
NumPixCheckSq = floor(sz1/NumCheckSquares);   % make it a proportion of the screen width.
DispCheckerboard = double(checkerboard(NumPixCheckSq, ...
    NumCheckSquares/2) > 0.5);                                % make the checkerboard of just 1s and 0s (white/black)
CheckerboardSizePix = size(DispCheckerboard);            % Set aside x,y size of checkerboard in pixels

%%
% We want to place the checkerboard against its grey background. 
% INPUT: B: Bigger matrix % b: small matrix that needs to put inside bigger matrix, B %OUTPUT: R: Resultant matrix % Example: % B=zeros(10,10); b=ones(5,5); % R=insertMatrix(B,b);
B = zeros(sz2);
[P,Q]=size(B); % The bigger matrix

fx=floor(P/2)-floor(size(DispCheckerboard,1)/2); 
fy=floor(Q/2)-floor(size(DispCheckerboard,2)/2);
R=B;
for p=1:size(DispCheckerboard,1)
    for q=1:size(DispCheckerboard,2)
        R(fx+p,fy+q)=DispCheckerboard(p,q);
    end
end

R1 = R;
% Make the 'full contrast' image
R1(~R) = 0.5; %make zeros=0.5. these stay fixed across images
R1(logical(R)) = 0.625; %25% diff against background
R2 = R1; 
R2(logical(R)) = 0.75; % 50% diff against background
R3 = R1; 
R3(logical(R)) = 0.85; % 85% diff against background

% Convert to RGB image
R1 = repmat(R1, 1, 1, 3);
R2 = repmat(R2, 1, 1, 3); 
R3 = repmat(R3, 1, 1, 3); 

figure
imshow(R1)
figure
imshow(R2)
figure
imshow(R3)


%%
% At this point we want to make another matrix with only the top and bottom rows 
% of the checkerboard showing, for the display of written feedback.

[P,Q,~]=size(R1); % The bigger matrix

% The smaller matrix: we want it to be 4 rows tall, all rows long
C = ones(NumPixCheckSq*4, size(DispCheckerboard,1)).*0.5; % all at background contrast

fx=floor(P/2)-floor(size(C,1)/2); 
fy=floor(Q/2)-floor(size(C,2)/2);
R=squeeze(R1(:,:,1));
for p=1:size(C,1)
    for q=1:size(C,2)
        R(fx+p,fy+q)=C(p,q);
    end
end

R4 = repmat(R,1,1,3);
figure
imshow(R4)

% Finally, make the blank image for start of trial:
Bl = repmat(ones(sz2),1,1,3).*0.5;
figure
imshow(Bl)

%% Now make reward mapping functions.
% Do monetary function:
r = 0; %target contrast
t = -15:0.5:15; % chosen contrasts

%t2 = ceil((0:30)/2);
%t2 = [fliplr(-t2(2:end)), t2];

f = 25 - 5*(abs(r-t))/3; % reward function

% Put in extra range (where no feedback given):
t2 = [-50:0.5:-15.5, t, 15.5:0.5:50];
e = length(t2)-length(f);
f2 = [zeros(1,e/2), f, zeros(1,e/2)];

figure
bar(t2,f2,1, 'EdgeColor', 'r', 'FaceColor', 'r', 'LineWidth', 2.5)
hold on
bar(t,f-0.1,1, 'EdgeColor', 'w', 'FaceColor', 'w', 'LineStyle','none')
% Cover up the red line at the bottom:
newt = -14.75:0.25:14.75;
plot(newt, zeros(1,length(newt)), 'w-', 'LineWidth', 3)
ylim([0 27])
xlim([-45 45])
% Do final formatting:
xlabel('Contrast difference from target (%)')
ylabel('Reward (cents)')
set(gca, 'TickDir','out', ...
    'box', 'off', ...
    'FontSize', 14)

%%
% Now make the same function, but for the instructive condition
figure
bar(t2,f2,1, 'EdgeColor', 'k', 'FaceColor', 'k', 'LineWidth', 2.5)
hold on
bar(t,f-0.1,1, 'EdgeColor', 'w', 'FaceColor', 'w', 'LineStyle','none')
% Cover up the red line at the bottom:
newt = -14.75:0.25:14.75;
plot(newt, zeros(1,length(newt)), 'w-', 'LineWidth', 3)
ylim([0 27])
xlim([-45 45])
% Do final formatting:
y_tick_lbl = fliplr(num2cell(linspace(0,15,6)));
y_tick_lbl{1} = '"Too far"';
xlabel('Contrast difference from target (%)')
ylabel('Instruction (%)')
set(gca, 'TickDir','out', ...
    'box', 'off', ...
    'YTickLabel', y_tick_lbl, ...
    'FontSize', 14)




