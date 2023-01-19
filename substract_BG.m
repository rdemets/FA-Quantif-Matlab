function [Im]=substract_BG(I,show,ifcontrast)
%%%thres is the threshold as a fraction of the background level: it should be
%adjusted around 1 (0.8 to 1.2 for example)
%%% average_level is the average pixel intensity in the cell compared to
% the maximum intensity for the class of images considered (for ex:1/3)
%%% show is a boolean to display the black and white mask or not
if ifcontrast
    M=contrast(I,0.1,0);
    %disp('contrast');
else
    M=I;
    %disp('no contrast');
end
%cl=class(M);
%create histogram
nbins=50;
Imax=double(max(max(M)));
n=hist(double(M(:)),(Imax/nbins/2:Imax/nbins:Imax*(1-1/nbins/2)));
% smooth histogram
span = 4;
w = ones(span,1)/span; 
smoothed_n = conv(n,w,'same');
%figure, plot(smoothed_n)
[V,imax]=max(smoothed_n);
background=Imax/nbins/2+Imax/nbins*imax;
%disp(['Intensite soustraite :',num2str(background)]);
Im=M-background;

end