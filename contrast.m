function [Iout,Imax]=contrast(Iin,Fsat,Fzero)
%Iin is the input image
%Fsat in the percentage of saturated pixels
%Fzero is the percentage of pixels at zero intensity
cl=class(Iin);
Nbpix=numel(Iin);
Iin=double(intmax(cl))/double(max(Iin(:))-min(Iin(:)))*double(Iin-min(Iin(:)));
Imax=intmax(cl);
Imin=0;
if strcmp(cl,'uint16')
   X=Iin>Imax;
   while sum(X(:))<0.01*Fsat*Nbpix
       Imax=Imax-100;
       X=Iin>Imax;
   end
   Y=Iin<Imin;
   while sum(Y(:))<0.01*Fzero*Nbpix
       Imin=Imin+100;
       Y=Iin<Imin;
   end
Iout=uint16(double(intmax(cl))/double(Imax-Imin)*double(Iin-Imin));
elseif strcmp(cl,'uint8')
   X=Iin>Imax;
   while sum(X(:))<0.01*Fsat*Nbpix
       Imax=Imax-5;
       X=Iin>Imax;
   end
   Y=Iin<Imin;
   while sum(Y(:))<0.01*Fzero*Nbpix
       Imin=Imin+5;
       Y=Iin<Imin;
   end
Iout=uint8(double(intmax(cl))/double(Imax-Imin)*double(Iin-Imin));
end   
