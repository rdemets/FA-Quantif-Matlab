function [Im,BW]=normalize(I, thres, average_level,show)
%%%thres is the threshold as a fraction of the background level: it should be
%adjusted around 1 (0.8 to 1.2 for example)
%%% average_level is the average pixel intensity in the cell compared to
% the maximum intensity for the class of images considered (for ex:1/3)
%%% show is a boolean to display the black and white mask or not
M=contrast(I,0.1,0);
cl=class(M);
%create histogram
nbins=50;
Imax=double(intmax(cl));
n=hist(double(M(:)),(Imax/nbins/2:Imax/nbins:Imax*(1-1/nbins/2)));
% smooth histogram
span = 4;
w = ones(span,1)/span; 
smoothed_n = conv(n,w,'same');
%figure, plot(smoothed_n)
[V,imax]=max(smoothed_n);
background=Imax/nbins/2+Imax/nbins*imax;
Mcorr=M-background;
%disp(background)
Mcorr(Mcorr<0)=0;
%determine threshold
BW1=Mcorr>thres*background;
%figure, imshow(BW1);
%keep only the largest feature in the binary image
BW=imclearborder(BW1);
CC=bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
for i=1:CC.NumObjects
    if i~=idx
BW(CC.PixelIdxList{i}) = 0;
    end
end
S=regionprops(BW,'BoundingBox','FilledImage');
start1=ceil(S.BoundingBox(2));
start2=ceil(S.BoundingBox(1));
width1=S.BoundingBox(4);
width2=S.BoundingBox(3);
BWfilled=false(size(Mcorr));
BWfilled(start1:start1+width1-1,start2:start2+width2-1)=S.FilledImage;
cellwin=Mcorr(BWfilled);
%  figure(3), subplot(1,2,1), imshow(Mcorr)
%  subplot(1,2,2), imshow(BWfilled)
if show
figure(15), imshow(BWfilled)
end
% hold on
% plot([S.BoundingBox(1) S.BoundingBox(1) S.BoundingBox(1)+S.BoundingBox(3) S.BoundingBox(1)+S.BoundingBox(3) S.BoundingBox(1)],[S.BoundingBox(2) S.BoundingBox(2)+S.BoundingBox(4) S.BoundingBox(2)+S.BoundingBox(4) S.BoundingBox(2) S.BoundingBox(2)],'g')
% hold off
Inorm=sum(cellwin)/(numel(cellwin)*Imax*average_level);
if strcmp(cl,'uint16')
    Im=uint16(double(Mcorr)/Inorm);
elseif strcmp(cl,'uint8')
    Im=uint8(double(Mcorr)/Inorm);
else
    Im=double(Mcorr)/Inorm;
end