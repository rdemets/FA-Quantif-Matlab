%%%%% Cette version propose une analyse des plaques d'adhesions
%%%%% par seuillage automatique,
%%%%% Organisation DATA : 1Cellule/Dossier, Fluo.tif et Mask.tif
%%%%% - d?termination de la distance ? la membrane
%%%%% - analyse des plaques (taille, forme, ...)
%%%%% Enregistrement des donn?es dans fichier dat
%%%%% Copie de graphs
%%%%% Threshold manuel, normalisation de l'intensite

%% Path
close all;
clc;
clear all;

% if ~exist('path','var')
%     path0='/Users/nus/Desktop/Papier in progress/Collaboration Lyon';
% end
%% Ouvre une image

% dirname = uigetdir(path0,'Choose directory ');
% files = dir(dirname);
% ind=cellfun(@find,{files(:).isdir},'UniformOutput',false);
% files=files(~cellfun(@isempty,ind));
% ind2=cellfun(@(s) isstrprop(s,'punct'),{files(:).name},'UniformOutput',false);
% files=files(cellfun(@sum,ind2)==0);
%file1=uipickfiles('FilterSpec','/Users/nus/Desktop/Papier in progress/Collaboration Lyon');

files=uipickfiles('FilterSpec','/Users/nus/Desktop/To analyse/Christelle');

Num_cell=length(files);
disp([num2str(Num_cell), ' cells found in directory.'])
premier_passage=0;
Tous=0;
path=files{1};
%%
% for i=1:Num_cell
%     a=files{i}.name;
%     
% end
%%
Conversion=input('Set Scale px-?m : 1pixel = ?\n');

%Pour la stat globale
Bilan=[];
Bilan_filtre=[];
Choix_fluo=0;


%% Boucle principale analyse en batch
for ii=2:Num_cell
    path0=files{ii};
    cd(path0)
    premier_passage=premier_passage+1;
    disp(['Cell n?',num2str(ii)]);
    tic
    NameFluo=fullfile(files(ii),'fluo.tif');
    Cellule=imread(NameFluo{1});
    if isa(Cellule,'uint8')
        Cellule=double(Cellule)/255;
    else
        Cellule=double(Cellule)/65535;
    end
    NameMask=fullfile(files(ii),'Mask.tif');
	BW=imread(NameMask{1});
	BW=logical(BW);
    %%
    
    
    % Binarisation de l'image selon une valeur de seuil
%     if Choix_fluo==0
%         Choix_fluo=input('Antibody labeling? Yes or No\n','s');
%     end
%     
%     
%     if strcmpi(Choix_fluo,'No')
%         avlevel=1/5;
%         [cell masque]=normalize(Cellule,0.6,avlevel,false);
%     else
%         
%         [cell]=substract_BG(Cellule,false,true);
%         [Cellule]=substract_BG(Cellule,false,false);
%         %3eme argument indique si on utilise la fonction 'contrast'
%         %constrast uniquement pour l'affichage, pas pour l'analyse
%         
%     end
    
    %%  Threshold
    
    if premier_passage==1
        threshold=thresh_tool(Cellule);        
    end
    Cell_bin=Cellule>threshold;
    %% Cherche le centroid et les extremas
    
    %trouve le centre g?om?trique
    s = regionprops(BW, 'centroid');
    centroids = cat(1, s.Centroid);
    
    %trouve les extremas
    s1 = regionprops(BW, 'extrema');
    extremas = cat(1, s1.Extrema);
    
    
    %% Dessine un disque centre au centroid et de rayon Rmax
    
    
    XCentre = centroids(1,1);
    YCentre = centroids(1,2);
    Centre=[XCentre YCentre];
    Nx = size(Cellule,2);
    Ny = size(Cellule,1);
    [x, y] = meshgrid(1:Nx, 1:Ny);
    R=zeros(size(extremas,1),1);
    % cherche le rayon Rmax du disque
    k=1;
    for j=1:size(extremas,1)
        R(k)= sqrt((extremas(j,1)-XCentre)^2+(extremas(j,2)-YCentre)^2);
        k=k+1;
    end
    Rmax = max (R);
    
    
    
    % dessine un disque de rayon Rmax pixels, centre en (XCentre,YCentre)
    cercleO = hypot(x-XCentre, y-YCentre)<Rmax;
    
    %% R?duction de la taille de l'image analys?e
    
    %Reduction de l'image centree sur le centroide pour limiter le calcul
    
    new_Cell_temp=imcrop(Cellule, [round(XCentre-Rmax) round(YCentre-Rmax) round(2*Rmax) round(2*Rmax)]);
    new_BW=imcrop(BW, [round(XCentre-Rmax) round(YCentre-Rmax) round(2*Rmax) round(2*Rmax)]);
    Cell_binaire=imcrop(Cell_bin, [round(XCentre-Rmax) round(YCentre-Rmax) round(2*Rmax) round(2*Rmax)]);
    
    %Elimination du signal provenant de l'exterieur de la membrane
    
    
    cl=class(new_Cell_temp);
    
    if strcmp(cl,'uint16')
        new_Cell=uint16(double(new_Cell_temp).*new_BW);
    elseif strcmp(cl,'uint8')
        new_Cell=uint8(double(new_Cell_temp).*new_BW);
    else
        new_Cell=double(double(new_Cell_temp).*new_BW);
    end
    
    Cell_binaire=double(Cell_binaire).*new_BW;
    
    
    
    %Filtre passe-haut => Detection de la membrane
    MembraneV = edge(new_BW,'sobel','horizontal');
    MembraneH = edge(new_BW,'sobel','vertical');
    Membrane=MembraneV+MembraneH;
    
    % Pour trouver le nombre de pixel composant la membrane
    Compteur=0;
    for i=1:size(Membrane,1)
        for j=1:size(Membrane,2)
            if Membrane(i,j)==1
                Compteur=Compteur+1;
            end
        end
    end
    
    %Recuperation des coordonnees de chaque pixel 'membrane'
    Coor_Membrane=zeros(Compteur,2);
    x=1;
    for i=1:size(Membrane,1)
        for j=1:size(Membrane,2)
            if Membrane(i,j)==1
                Coor_Membrane(x,1)=i;
                Coor_Membrane(x,2)=j;
                x=x+1;
            end
        end
    end
    
    
    % Detecte les 'objets connect?s' ? partir d'une image binaire, et calcule le
    % centroide de chaque objet
    CC_temp=bwconncomp(Cell_binaire);
    Taille_temp=regionprops(CC_temp,new_Cell,'MajorAxisLength');
    Taille_temp=cat(1,Taille_temp.MajorAxisLength);
    j=1;
    CC=CC_temp;
    
    %Suppression des plaques de 1pixel
    for i=1:size(Taille_temp,1)
        if Taille_temp(i)>min(min(Taille_temp))
            CC.PixelIdxList(1,j)=CC_temp.PixelIdxList(1,i);
            j=j+1;
        end
    end
    j=j-1;
    CC.PixelIdxList=CC.PixelIdxList(1,1:j);
    CC.NumObjects=j;
    
    Area=regionprops(CC,new_Cell,'Area');
    Area=cat(1,Area.Area);
    
    Shape=regionprops(CC,new_Cell,'Eccentricity'); % Shape factor
    Shape=cat(1,Shape.Eccentricity);
    
    Longueur=regionprops(CC,new_Cell,'MajorAxisLength');
    Longueur=cat(1,Longueur.MajorAxisLength);
    
    Intensity=regionprops(CC,new_Cell,'MeanIntensity');
    Intensity=cat(1,Intensity.MeanIntensity);
    
    %Matrice contenant la position de chaque plaque
    Plaques=regionprops(CC,'Centroid');
    Position_temp=cat(1,Plaques.Centroid);
    
    % Inversion des coordonn?es X,Y
    Position_Plaques=zeros(size(Position_temp,1),2);
    
    % Si des plaques sont d?tect?s
    if ~isempty(Position_temp)
        Position_Plaques(:,1)=Position_temp(:,2);
        Position_Plaques(:,2)=Position_temp(:,1);
        %% Calcule des distances
        
        %Catch Warning message
        warning('OFF', 'MATLAB:NonScalarInput');
        
        
        Centre = regionprops(new_BW, 'centroid');
        Centre = cat(1, Centre.Centroid);
        
        DistanceR_Centre=zeros(size(Coor_Membrane,1),1);
        % Calcul de Rmin pour la normalisation :
        for j=1:size(Coor_Membrane,1)
            DistanceR_Centre(j)=sqrt((Coor_Membrane(j,1)-Centre(2))^2+(Coor_Membrane(j,2)-Centre(1))^2);
        end
        %Distance minimal entre le centre et la membrane la plus proche
        Rmin= min(DistanceR_Centre);
        
        
        
        resultat=zeros(size(Position_Plaques,1),1);
        Resultat_temp=zeros(size(Coor_Membrane,1),1);
        for i=1:size(Position_Plaques)
            % Calcule de la distance entre la plaque et tout les pixels
            % 'membrane'
            for j=1:size(Coor_Membrane,1)
                Resultat_temp(j)=sqrt((Position_Plaques(i,1)-Coor_Membrane(j,1))^2+(Position_Plaques(i,2)-Coor_Membrane(j,2))^2);
            end
            %On sauvegarde la distance minimale
            resultat(i)=min(Resultat_temp);
        end
        
        % Prise en compte du facteur de conversion pour le calcul Aire et
        % Longueur
        for i=1:size(resultat)
            resultat(i)=resultat(i)*Conversion;
            temp=sqrt(Area(i));
            temp=temp*Conversion;
            Area(i)=temp^2;
            Longueur(i)=Longueur(i)*Conversion;
        end
        
        % Normalisation par le Rayon Min de la cellule
        resultat_norm=zeros(size(resultat,1),1);
        
        %Si cellule pinc?e, possible que plaque plus ?loign? que Rmin.
        % => uniquement pour garder la m?me dynamique pour toutes les cellules
        if max(resultat)>Rmin*Conversion;
            Rmin=max(resultat)/Conversion;
        end
        
        for i=1:size(resultat)
            resultat_norm(i)=resultat(i)/Rmin/Conversion*100;
        end
        
        
        %% Resultat
        resultat(i+1:end)=[]; %efface tout ce qui est ? zero
        %%
        figure('Name','Focal Adhesion Distribution','NumberTitle','off');
        subplot(2,2,1)
        imshow(new_Cell,[min(min(new_Cell)) max(max(new_Cell))])
        title('Cell intensity');
        hold on
        plot(Position_Plaques(:,2),Position_Plaques(:,1),'og')
        plot(Coor_Membrane(:,2),Coor_Membrane(:,1),'+r','LineWidth',0.1)
        hold off
        
        subplot(2,2,2)
        imshow(Cell_binaire);
        title('Binarised cell');
        hold on
        plot(Position_Plaques(:,2),Position_Plaques(:,1),'og')
        plot(Coor_Membrane(:,2),Coor_Membrane(:,1),'+r','LineWidth',0.1)
        hold off
        subplot(2,2,3)
        hist(resultat);
        xlabel(['Distance (?m)    N = ',num2str(size(resultat,1))])
        ylabel('count')
        titre='Histogram';
        
        subplot(2,2,4)
        hist(resultat_norm,10:10:90);
        xlabel(['Relative distance to the membrane  N = ',num2str(size(resultat,1))])
        ylabel('count')
        
        warning('OFF', 'MATLAB:MKDIR:DirectoryExists');
        mkdir(fullfile(path0,'Summary'));
        newPath=[path0,'/Summary'];
        % /Summary pour Mac
        cd(newPath)
        print ('-dtiff',titre);
        
        %% Enregistrements
        
        resultswrite=fopen(('ResultsDistrib.dat'),'w');
        fprintf(resultswrite,'FA Number\t');
        fprintf(resultswrite,'Area(?m?)\t');
        fprintf(resultswrite,'Mean Intensity\t');
        fprintf(resultswrite,'Membrane_Distance(?m)\t');
        fprintf(resultswrite,'Shape Factor\t');
        fprintf(resultswrite,'Length (?m)\t');
        fprintf(resultswrite,'Relative Distance\t');
        fprintf(resultswrite,'\n');
        for i=1:size(resultat)
            fprintf(resultswrite,[num2str(i),'\t']);
            fprintf(resultswrite,[num2str(Area(i)),'\t']);
            fprintf(resultswrite,[num2str(Intensity(i)),'\t']);
            fprintf(resultswrite,[num2str(resultat(i)),'\t']);
            fprintf(resultswrite,[num2str(Shape(i)),'\t']);
            fprintf(resultswrite,[num2str(Longueur(i)),'\t']);
            fprintf(resultswrite,[num2str(resultat_norm(i)),'\t']);
            fprintf(resultswrite,'\n');
        end
        fclose(resultswrite);
        
        Bilan_temp=[];
        Bilan_temp(1:size(Area,1),1)=ii;
        Bilan_temp=cat(2,Bilan_temp,Area);
        Bilan_temp=cat(2,Bilan_temp,Intensity);
        Bilan_temp=cat(2,Bilan_temp,resultat);
        Bilan_temp=cat(2,Bilan_temp,Shape);
        Bilan_temp=cat(2,Bilan_temp,Longueur);
        Bilan_temp=cat(2,Bilan_temp,resultat_norm);
        Bilan=cat(1,Bilan,Bilan_temp);
        Bilan_temp=[];
        
        
        %% Elimination d'objet
        if premier_passage==1 || Tous==0
            Choix = questdlg('FA Discrimination?', 'Discrimination','Oui', 'Non','Non');
            if strcmp(Choix,'Oui')
                Choix=1;
            else
                Choix=0;
            end
            Tous = questdlg('Do the same for all the images?', 'Discrimination','Oui', 'Non','Non');
            if strcmp(Tous,'Oui')
                Tous=1;
            else
                Tous=0;
            end
        end
        close;
        if Choix==1
            j=1;
            resultat_filtre=zeros(size(resultat,1),1);
            Intensity_filtre=zeros(size(Intensity,1),1);
            Area_filtre=zeros(size(Area,1),1);
            Position_Plaques_filtre=zeros(size(Position_Plaques,1),2);
            Shape_filtre=zeros(size(Shape,1),1);
            Longueur_filtre=zeros(size(Longueur,1),1);
            resultat_norm_filtre=zeros(size(resultat_norm,1),1);
            if Tous==0 || premier_passage==1
                Choix2 = questdlg('Discrimination by :','Choice', ...
                    'Area','Shape','Area');
                % Handle response
                switch Choix2
                    case 'Area'
                        Choix2=0;
                    case 'Shape'
                        Choix2=1;
                end
            end
            if Choix2==0
                if Tous==0 || premier_passage==1
                    premier_passage=2;
                    reponse = {'Min Size (?m)','Max Size (?m)',};
                    x = inputdlg(reponse,'Size Choice');
                    Critere_min=str2double(x{1});
                    Critere_max=str2double(x{2});
                end
                for i=1:size(resultat)
                    % Copie des plaques ayant une aire sup?rieur au critere
                    if Area(i)>Critere_min && Area(i)<Critere_max
                        resultat_filtre(j)=resultat(i);
                        Area_filtre(j)=Area(i);
                        Intensity_filtre(j)=Intensity(i);
                        Position_Plaques_filtre(j,1)=Position_Plaques(i,1);
                        Position_Plaques_filtre(j,2)=Position_Plaques(i,2);
                        Shape_filtre(j)=Shape(i);
                        Longueur_filtre(j)=Longueur(i);
                        resultat_norm_filtre(j)=resultat_norm(i);
                        j=j+1;
                    end
                end
            else
                if Tous==0 || premier_passage==1
                    premier_passage=2;
                    reponse = {'Shape Factor Min ','Shape Factor Max',};
                    x = inputdlg(reponse,'Size Choice');
                    Critere1=str2double(x{1});
                    Critere2=str2double(x{2});
                end
                for i=1:size(resultat)
                    if Shape(i)>=Critere1 && Shape(i)<=Critere2
                        resultat_filtre(j)=resultat(i);
                        Area_filtre(j)=Area(i);
                        Intensity_filtre(j)=Intensity(i);
                        Position_Plaques_filtre(j,1)=Position_Plaques(i,1);
                        Position_Plaques_filtre(j,2)=Position_Plaques(i,2);
                        Shape_filtre(j)=Shape(i);
                        Longueur_filtre(j)=Longueur(i);
                        resultat_norm_filtre(j)=resultat_norm(i);
                        j=j+1;
                    end
                end
            end
            
            % Efface tout ce qui est ? zero
            Area_filtre(Area_filtre==0)=[];
            taille=size(Area_filtre);
            
            % Copie uniquement des r?sultats r?pondant aux criteres
            resultat_filtre=resultat_filtre(1:taille);
            resultat_norm_filtre=resultat_norm_filtre(1:taille);
            Intensity_filtre=Intensity_filtre(1:taille);
            Shape_filtre=Shape_filtre(1:taille);
            Longueur_filtre=Longueur_filtre(1:taille);
            
            % Sauvegarde des donn?es pour concat?ner les r?sultats sur toutes
            % les cellules
            Bilan_temp=[];
            Bilan_temp(1:size(Area_filtre,1),1)=ii;
            Bilan_temp=cat(2,Bilan_temp,Area_filtre);
            Bilan_temp=cat(2,Bilan_temp,Intensity_filtre);
            Bilan_temp=cat(2,Bilan_temp,resultat_filtre);
            Bilan_temp=cat(2,Bilan_temp,Shape_filtre);
            Bilan_temp=cat(2,Bilan_temp,Longueur_filtre);
            Bilan_temp=cat(2,Bilan_temp,resultat_norm_filtre);
            Bilan_filtre=cat(1,Bilan_filtre,Bilan_temp);
            Bilan_temp=[];
            
            
            %% Affichage
            figure('Name','Focal Adhesion distribution (filtre)','NumberTitle','off');
            subplot(2,2,1)
            imshow(new_Cell)
            title('Cell intensity');
            hold on
            plot(Position_Plaques_filtre(:,2),Position_Plaques_filtre(:,1),'og')
            plot(Coor_Membrane(:,2),Coor_Membrane(:,1),'+r','LineWidth',0.1)
            hold off
            subplot(2,2,2)
            imshow(Cell_binaire)
            title('Binarised cell');
            hold on
            plot(Position_Plaques_filtre(:,2),Position_Plaques_filtre(:,1),'og')
            plot(Coor_Membrane(:,2),Coor_Membrane(:,1),'+r','LineWidth',0.1)
            hold off
            subplot(2,2,3)
            hist(resultat_filtre);
            xlabel(['Distance (?m)    N = ',num2str(size(resultat_filtre,1))])
            ylabel('count')
            titre2='Histogram_filtered';
            
            subplot(2,2,4)
            hist(resultat_norm_filtre,10:10:90);
            xlabel(['Relative distance to the membrane  N = ',num2str(size(resultat_norm_filtre,1))])
            ylabel('count')
            print ('-dtiff',titre2);
            close;
            
            
            resultswrite2=fopen(('ResultsDistrib_filtre.dat'),'w');
            fprintf(resultswrite2,'FA Number\t');
            fprintf(resultswrite2,'Area(?m?)\t');
            fprintf(resultswrite2,'Mean Intensity\t');
            fprintf(resultswrite2,'Membrane Distance(?m)\t');
            fprintf(resultswrite2,'Shape Factor\t');
            fprintf(resultswrite2,'Length (?m)\t');
            fprintf(resultswrite2,'Relative Distance\t');
            fprintf(resultswrite2,'\n');
            
            for i=1:size(resultat_filtre)
                fprintf(resultswrite2,[num2str(i),'\t']);
                fprintf(resultswrite2,[num2str(Area_filtre(i)),'\t']);
                fprintf(resultswrite2,[num2str(Intensity_filtre(i)),'\t']);
                fprintf(resultswrite2,[num2str(resultat_filtre(i)),'\t']);
                fprintf(resultswrite2,[num2str(Shape_filtre(i)),'\t']);
                fprintf(resultswrite2,[num2str(Longueur_filtre(i)),'\t']);
                fprintf(resultswrite2,[num2str(resultat_norm_filtre(i)),'\t']);
                fprintf(resultswrite2,'\n');
            end
            fclose(resultswrite2);
            
            
        end
    end
end


%% Copie du bilan sur toutes les cellules

% Donn?es brutes
cd(path);
resultswrite3=fopen(('Bilan.dat'),'w');
fprintf(resultswrite3,'Cell Number\t');
fprintf(resultswrite3,'FA Number\t');
fprintf(resultswrite3,'Area(?m?)\t');
fprintf(resultswrite3,'Mean Intensity\t');
fprintf(resultswrite3,'Membrane Distance(?m)\t');
fprintf(resultswrite3,'Shape Factor\t');
fprintf(resultswrite3,'Length (?m)\t');
fprintf(resultswrite3,'Relative Distance\t');
fprintf(resultswrite3,'\n');
for i=1:size(Bilan,1)
    fprintf(resultswrite3,[num2str(Bilan(i,1)),'\t']);
    fprintf(resultswrite3,[num2str(i),'\t']);
    fprintf(resultswrite3,[num2str(Bilan(i,2)),'\t']);
    fprintf(resultswrite3,[num2str(Bilan(i,3)),'\t']);
    fprintf(resultswrite3,[num2str(Bilan(i,4)),'\t']);
    fprintf(resultswrite3,[num2str(Bilan(i,5)),'\t']);
    fprintf(resultswrite3,[num2str(Bilan(i,6)),'\t']);
    fprintf(resultswrite3,[num2str(Bilan(i,7)),'\t']);
    fprintf(resultswrite3,'\n');
end
fclose(resultswrite3);

% Construction du tableau repr?sentant l'intensit? en fonction de la
% distance (echantillonage 1?m)

bin=1:ceil(max(Bilan(:,4)));
hist_intensity=zeros(length(bin),2);
hist_intensity_norm=zeros(100,2);
hist_intensity(:,2)=bin;
hist_intensity_norm(:,2)=1:100;
E=zeros(1,size(bin,2));
ENorm=zeros(1,100);


for i=1:size(bin,2);
    [row]=find(ceil(Bilan(:,4)) > (bin(i)-1) & ceil(Bilan(:,4))<= bin(i));
    hist_intensity(i)=mean(Bilan(row,3));
    CalcErr=Bilan(row,3);
    E(i)=std(CalcErr,1);
end

for i=1:100;
    [row]=find(ceil(Bilan(:,7)) > (i-1) & ceil(Bilan(:,7))<=i);
    hist_intensity_norm(i)=mean(Bilan(row,3));
    CalcErrNorm=Bilan(row,3);
    ENorm(i)=std(CalcErrNorm,1);
end



% Affichage

h2=figure('Name','Focal Adhesion Distribution','NumberTitle','off');
set(h2,'Units', 'Normalized', 'Position',[0.02 0.02 0.95 0.9])
subplot(3,2,1)
hist(Bilan(:,4));
title(['Number of FA = ',num2str(size(Bilan,1)),' , Number of cell = ',num2str(Num_cell)]);
xlabel('Distance (?m)');
ylabel('count');


subplot(3,2,2)
hist(Bilan(:,7));
xlabel('Relative distance to the membrane');
ylabel('count');




subplot(3,2,3);
h=errorbar(hist_intensity(:,2),hist_intensity(:,1),E,'+b');
hc = get(h, 'Children');
%set(hc(1),'color','+b') %// data
%set(hc(2),'color','r') %// error bars
axis([0 ceil(max(max(Bilan(:,4))))+1 0 2*max(Intensity)])
%bar((hist_intensity(:,1)));
xlabel('Distance (?m)');
ylabel('Intensity');

subplot(3,2,4)
h=errorbar(hist_intensity_norm(:,2),hist_intensity_norm(:,1),ENorm,'+b');
hc = get(h, 'Children');
%set(hc(1),'color','+b') %// data
%set(hc(2),'color','r') %// error bars
axis([0 100 0 2*max(Intensity)])
%bar((hist_intensity_norm(:,1)));
xlabel('Relative distance to the membrane');
ylabel('Intensity');



subplot(3,2,5)

% 30 barres entre le Max et le Min
ecart=(max(Area)-min(Area))/30;
hist(Bilan(:,2),0:ecart:max(Bilan(:,2)));
xlabel('Area (?m?)');
ylabel('count');

subplot(3,2,6)
ecart2=(max(Longueur)-min(Longueur))/30;
hist(Bilan(:,6),0:ecart2:max(Bilan(:,6)));
xlabel('Max length of the FA (?m)');
ylabel('count');
titre2='Histogram_bilan';
print ('-dtiff',titre2);
hgsave('Histogram_bilan.fig');

%% Copie histo_intensite
cd(path);
resultsbins=fopen(('Histogram_intensity.dat'),'w');
fprintf(resultsbins,'Bins (?m)\t');
fprintf(resultsbins,'Intensity\t');
fprintf(resultsbins,'STD\n');
for i=1:size(bin,2)
    fprintf(resultsbins,[num2str(hist_intensity(i,2)-1),'-',num2str(hist_intensity(i,2)),'\t']);
    fprintf(resultsbins,[num2str(hist_intensity(i,1)),'\t']);
    fprintf(resultsbins,[num2str(E(i)),'\n']);
end

fprintf(resultsbins,'\n\n\n\n');
fprintf(resultsbins,'Bins (%%) \t');
fprintf(resultsbins,'Intensity\t');
fprintf(resultsbins,'STD\n');

for i=1:100
    fprintf(resultsbins,[num2str(i-1),'-',num2str(i),'\t']);
    fprintf(resultsbins,[num2str(hist_intensity_norm(i,1)),'\t']);
    fprintf(resultsbins,[num2str(ENorm(i)),'\n']);
end

%% Filtre


if Choix==1
    bin_filtre=1:(ceil(max(Bilan_filtre(:,4)))+1);
    EFiltre=zeros(1,size(bin_filtre,2));
    EFiltreNorm=zeros(1,100);
    hist_intensity_filtre=zeros(length(bin_filtre),2);
    hist_intensity_filtre_norm=zeros(100,2);
    hist_intensity_filtre(:,2)=bin_filtre;
    hist_intensity_filtre_norm(:,2)=1:100;
    
    
    for i=1:size(bin_filtre,2);
        [row]=find(ceil(Bilan_filtre(:,4))>(bin_filtre(i)-1) &ceil(Bilan_filtre(:,4))<= bin_filtre(i));
        hist_intensity_filtre(i)=mean(Bilan_filtre(row,3));
        CalcErrFiltre=Bilan_filtre(row,3);
        EFiltre(i)=std(CalcErrFiltre,1);
    end
    
    for i=1:100;
        [row]=find(ceil(Bilan_filtre(:,7))>(i-1) & ceil(Bilan_filtre(:,7))<= i);
        hist_intensity_filtre_norm(i)=mean(Bilan_filtre(row,3));
        CalcErrFiltreNorm=Bilan_filtre(row,3);
        EFiltreNorm(i)=std(CalcErrFiltreNorm,1);
    end
    
    h1=figure('Name','Focal Adhesion distribution (filtre)','NumberTitle','off');
    set(h1,'Units', 'Normalized', 'Position',[0.02 0.02 0.95 0.9])
    subplot(3,2,1)
    hist(Bilan_filtre(:,4));
    title(['Number of FA = ',num2str(size(Bilan_filtre,1)),' , Number of cell = ',num2str(Num_cell)]);
    xlabel('Distance (?m)');
    ylabel('count')
    
    
    subplot(3,2,2)
    hist(Bilan_filtre(:,7),10:10:90);
    xlabel('Relative distance to the membrane');
    ylabel('count');
    
    
    subplot(3,2,3);
    h=errorbar(hist_intensity_filtre(:,2),hist_intensity_filtre(:,1),EFiltre,'+b');
    hc = get(h, 'Children');
    %set(hc(1),'color','+b') %// data
    %set(hc(2),'color','r') %// error bars
    axis([0 ceil(max(max(Bilan_filtre(:,4))))+1 0 2*max(Intensity)])
    %set(gca,'XTick',0:max(hist_intensity_filtre(:,2)));
    %bar((hist_intensity_filtre(:,1)));
    xlabel('Distance (?m)');
    ylabel('Intensity');
    
    subplot(3,2,4)
    h=errorbar(hist_intensity_filtre_norm(:,2),hist_intensity_filtre_norm(:,1),EFiltreNorm,'+b');
    hc = get(h, 'Children');
    %set(hc(1),'color','+b') %// data
    %set(hc(2),'color','r') %// error bars
    axis([0 100 0 2*max(Intensity)])
    %set(gca,'XTick',0:10:max(hist_intensity_filtre_norm(:,2)));
    %bar((hist_intensity_filtre_norm(:,1)));
    xlabel('Relative distance to the membrane');
    ylabel('Intensity');
    
    subplot(3,2,5)
    
    ecart=(max(Area_filtre)-min(Area_filtre))/30;
    hist(Area_filtre,0:ecart:max(Area_filtre));
    xlabel('Area (?m?)');
    ylabel('count');
    
    subplot(3,2,6)
    ecart2=(max(Longueur_filtre)-min(Longueur_filtre))/30;
    hist(Longueur_filtre,0:ecart2:max(Longueur_filtre));
    xlabel('Max length of the FA (?m)');
    ylabel('count');
    titre3='Histogram_bilan_filtered.tif';
    print ('-dtiff',titre3);
    hgsave('Histogram_bilan_filtered.fig');
    
    
    
    % Copie bilan avec filtrage
    cd(path);
    resultswrite4=fopen(('Bilan_filtre.dat'),'w');
    fprintf(resultswrite4,'Cell Number\t');
    fprintf(resultswrite4,'FA Number\t');
    fprintf(resultswrite4,'Area(?m?)\t');
    fprintf(resultswrite4,'Mean Intensity\t');
    fprintf(resultswrite4,'Membrane Distance(?m)\t');
    fprintf(resultswrite4,'Shape Factor\t');
    fprintf(resultswrite4,'Length (?m)\t');
    fprintf(resultswrite4,'Relative Distance\t');
    fprintf(resultswrite4,'\n');
    
    for i=1:size(Bilan_filtre,1)
        fprintf(resultswrite4,[num2str(Bilan_filtre(i,1)),'\t']);
        fprintf(resultswrite4,[num2str(i),'\t']);
        fprintf(resultswrite4,[num2str(Bilan_filtre(i,2)),'\t']);
        fprintf(resultswrite4,[num2str(Bilan_filtre(i,3)),'\t']);
        fprintf(resultswrite4,[num2str(Bilan_filtre(i,4)),'\t']);
        fprintf(resultswrite4,[num2str(Bilan_filtre(i,5)),'\t']);
        fprintf(resultswrite4,[num2str(Bilan_filtre(i,6)),'\t']);
        fprintf(resultswrite4,[num2str(Bilan_filtre(i,7)),'\t']);
        fprintf(resultswrite4,'\n');
    end
    cd(path);
    resultsbinsfiltre=fopen(('Histogram_intensity_filtre.dat'),'w');
    fprintf(resultsbinsfiltre,'Bins (?m)\t');
    fprintf(resultsbinsfiltre,'Intensity\t');
    fprintf(resultsbinsfiltre,'STD\n');
    for i=1:size(bin_filtre,2)
        fprintf(resultsbinsfiltre,[num2str(hist_intensity_filtre(i,2)-1),'-',num2str(hist_intensity_filtre(i,2)),'\t']);
        fprintf(resultsbinsfiltre,[num2str(hist_intensity_filtre(i,1)),'\t']);
        fprintf(resultsbinsfiltre,[num2str(EFiltre(i)),'\n']);
    end
    
    fprintf(resultsbinsfiltre,'\n\n\n\n');
    fprintf(resultsbinsfiltre,'Bins (%%) \t');
    fprintf(resultsbinsfiltre,'Intensity\t');
    fprintf(resultsbinsfiltre,'STD\n');
    
    for i=1:100
        fprintf(resultsbinsfiltre,[num2str(i-1),'-',num2str(i),'\t']);
        fprintf(resultsbinsfiltre,[num2str(hist_intensity_filtre_norm(i,1)),'\t']);
        fprintf(resultsbinsfiltre,[num2str(EFiltreNorm(i)),'\n']);
    end
    
    fclose(resultswrite4);
    fclose(resultsbinsfiltre);
    
    
end

fclose(resultsbins);
