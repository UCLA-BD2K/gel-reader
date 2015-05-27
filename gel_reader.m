clear all;
close all;
gelImage = real(imread('raw_gel.jpg'));
bw = rgb2gray(gelImage);

%these are the known molecular weights of the ladder proteins (in kDa)
%the ladder is the leftmost column
knownLadderValues = [250 150 100 80 60 50 40 30 25 20 15 10];

[imageY, imageX] = size(bw);

%correct for background color
background = max(max(bw))
bw = background - bw;

imshow(bw); 
hold all;

%we plot the maximum value for each pixel-column (note that the pixel 0,0
%is in the top left corner of an image, so graphs correspond inversely to
%images. also note that this method of distinguishing between columns may
%not work if columns are significantly skewed
maximums = double(max(bw));
plot(maximums);

%find the significant peaks (i.e., the midpoints between the columns
midpoints = peakfinder(-1*maximums);

for i = 1:size(midpoints, 2)
    plot( midpoints(i), maximums(midpoints(i)), 'rx');
end

numberOfColumns = size(midpoints, 2) - 1

%let's say that the data of interest is the middle 1/3 between two
%midpoints (this is an arbitrary decision of mine)
bounds = zeros(2, numberOfColumns)
for i = 1:numberOfColumns
   onethird = (midpoints(i+1) - midpoints(i))/3;
   bounds(:,i) = [midpoints(i) + onethird, midpoints(i) + 2*onethird];
end

for i = 1:numberOfColumns
   plot([bounds(1,i) bounds(1,i)], [1 imageY]);
   plot([bounds(2,i) bounds(2,i)], [1 imageY]);
end

%the integrals are all the rows summed for each subsection
integrals = zeros(numberOfColumns, imageY);
for i = 1:numberOfColumns
    start = floor(bounds(1,i));
    stop = ceil(bounds(2,i));
    strip = bw(:,start:stop);
    %figure, imshow(strip);
    integrals(i,:) = sum(strip, 2)./(stop-start);
end

%the first column integrated is the ladder, so lets find the peaks. each
%one corresponds to a known molecular weight
ladderPoints = peakfinder(integrals(1,:));

%lets lower the selectivity until we have exactly one peak for every band
%in the ladder
selectivity = (max(integrals(1,:))-min(integrals(1,:)))/4;

while(size(ladderPoints,2) < size(knownLadderValues,2))
    selectivity = selectivity*0.9;
    ladderPoints = peakfinder(integrals(1,:), selectivity);
    size(ladderPoints);
end

%plot the peaks and first column integral
figure;
hold on;
plot(integrals(1,:));
for i = 1:size(ladderPoints, 2)
    plot( ladderPoints(i), integrals(1,ladderPoints(i)), 'rx');
end

figure;
hold on;
for i = 1:size(ladderPoints, 2)
    plot( ladderPoints(i), knownLadderValues(i), 'rx');
end

%transpose
ladderPoints = ladderPoints.';
knownLadderValues = knownLadderValues.';

%fit to 1/x
modelFunction =  @(c,t)(c(1).*(t.^c(2)));
coeffGuesses = [1000 -1];
[coeffEstimates,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(ladderPoints, knownLadderValues, modelFunction, coeffGuesses);

bestFitFunction = @(x)coeffEstimates(1)*x.^coeffEstimates(2);

x = ladderPoints(1):imageY;
y = bestFitFunction(x);
plot(x,y);

molecularWeights = bestFitFunction(1:imageY);

testSlice = integrals(2,:) - min(integrals(2,:));
%pick top X peaks, find the global max, move left and right until you find
%a minimum, fit the points on each side of the peak to a straight line,
%subtract from the full slice, do another iteration

compounds = peakfinder(testSlice);

figure;
hold on;
plot(testSlice);
for i = 1:size(compounds, 2)
    plot( compounds(i), testSlice(compounds(i)), 'rx');
end

% peaksToExtract = 10;
% for i = 1:10
%     [highestPeakY, highestPeakX] = max(testSlice);
%     lowestPointToRight = highestPeakX;
%     for ii = highestPeakX:imageY
%         if(testSlice(ii) > testSlice(lowestPointToRight))
%             break;
%         end
%         lowestPointToRight = testSlice(ii);
%     end
%     lowestPointToLeft = highestPeakX;
%     for ii = highestPeakX:-1:0
%         if(testSlice(ii) > testSlice(lowestPointToLeft))
%             break;
%         end
%         lowestPointToLeft = testSlice(ii);
%     end
%     
% end
%peakfit(testSlice,floor(imageY/2),imageY,50,21);