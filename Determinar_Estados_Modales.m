%% Cargar el video y obtener una celda con cada uno de los frames
vidObj = VideoReader('incremento_FBG.avi');%Se importa el video para trabajar con el 
numFrames = 0;%se crea la variable que contará los frames
ImFrames = cell([],1) ;%se crea la celda que almacenará los frames

while hasFrame(vidObj) %ciclo WHILE para recorrer el objeto
    
    F = readFrame(vidObj);    % se leen los frames
    numFrames = numFrames + 1; %Incremento
    ImFrames{numFrames} = F(:,:,1) ; %se almacena solamente la componente R
                                     %del formato RGB de imagen
    %obteniendo el tamaño de cada frame
    
    [filas,columnas]=size(ImFrames{numFrames});
    %recortar la imagen
    rect=[ filas/2 columnas/2 filas columnas];
    ImFrames{numFrames}=imcrop(ImFrames{numFrames}, rect);
end


%% Procesamiento morfologico de la imagen (filtrado circular)
filtro_circ = fspecial('disk',3); %filtro de promediado circular

for i=1:1:numFrames
    ImFrames{i} = imbinarize(ImFrames{i}); %Hacer la imagen binaria 
    ImFrames{i} = imfilter(ImFrames{i},filtro_circ,'replicate');%filtrado circular
    ImFrames{i}=imclearborder(ImFrames{i}); %limpio los bordes    
end

%% Estados modales mediante similitud entre frames
DICE=zeros (1,numFrames-1); %se crea la variable 

for i=2:1:numFrames %se comienza en dos para no usar el primer frame usado anteriormente
    % Metodo dice para determinar la similitud de dos imagenes
    unionMask = ImFrames{1}|ImFrames{i};
    andMask = ImFrames{1}&ImFrames{i};
    DICE(i-1) = nnz(andMask) / nnz(unionMask);    
end
 DICE_inv = 1-DICE; % se invierte para trabajo posterior
%% Estados modales mediante Posicion e Intensidad de un ROI

%Definiendo la máscara que corresponde a la mancha de mayor tamano
reg=regionprops(ImFrames{1}); %Obtiene las propiedades del Frame 1
Areas=[reg.Area]; %Areas de las regiones del Frame 1
[AreaMax,Indice]=max(Areas);%Se obtiene la mayor área
Rec_min=reg(Indice).BoundingBox; %Bounding Box correspondiente a la mayor area 
Rec_min = round(Rec_min); %redondea el valor

diferencias = zeros (1,numFrames-1); %creo el arreglo de diferencias de intensidades
for i=2:1:numFrames %se comienza en dos para no usar el primer frame usado anteriormente
    %Area en el nuevo frame
    IMG=ImFrames{i};
    ROI=IMG(Rec_min(2):Rec_min(2)+Rec_min(4),Rec_min(1):Rec_min(1)+Rec_min(3));
    AreaROI=sum(ROI(:)); %Obtiene el área de la región de interés
    dif=abs(AreaMax-AreaROI);%Diferencias 
    diferencias(i-1)=dif;
end

%% Estados Modales generalizacion

 Kps = diferencias.*DICE_inv;


 %filtrado
 order = 5;
 framelen = 11;
 Kps = sgolayfilt(Kps,order,framelen);

 %% promediando para obtener la variazion de estados vs tiempo
 dif_mean = zeros (1,round(numFrames/vidObj.FrameRate));
 %se comienza en dos para no usar el primer frame usado anteriormente
 for i=1:1:round(numFrames/vidObj.FrameRate) 
     fin = 25*i;
     if i == 1
         dif_mean(i) = mean(Kps(1:fin));
     elseif i~=30
         dif_mean(i) = mean(Kps(ini:fin));
     else
         dif_mean(i) = mean(Kps(ini:fin-1));
     end
     ini = fin;
 end

 

 %% Figuras
 figure()
 x = [1:1:numFrames-1];
 plot(x,abs(Kps),'b',LineWidth=1.9)
 title('Factor de caracterizacion de los estados modales')
 grid on
 xlabel('Frames')
 ylabel('Factor de diferencia')



 figure()
 x = [1:1:30];
 plot(x,dif_mean,'r',LineWidth=1.9)
 title('Factor de caracterizacion de los estados modales')
 grid on
 xlabel('Tiempo en segundos')
 ylabel('Factor de diferencia')
 
 %% Usando las cinco areas mayores
 diferenciasAreas = zeros (5,numFrames-1);
 DICE=zeros (1,numFrames-1);

 %Definiendo la máscara que corresponde a la mancha de mayor tamano
 reg=regionprops(ImFrames{1}); %Obtiene las propiedades de las regiones
 Areas=[reg.Area]; %Crea un arreglo con las áreas de las regiones
 Areas_ord = sort(Areas,'descend');

 for a=1:1:5
     [Indice]=find(Areas==Areas_ord(a)) ;%Obtiene la region de mayor área
     Rec_min = reg(Indice).BoundingBox; %define el mínimo rectángulo en el que está contenida la región
     Rec_min = round(Rec_min); %redondea el valor

     for i=2:1:numFrames %se comienza en dos para no usar el primer frame usado anteriormente
         IMG=ImFrames{i};
         ROI=IMG(Rec_min(2):Rec_min(2)+Rec_min(4),Rec_min(1):Rec_min(1)+Rec_min(3));%Define la región de interés en el nuevo frame
         AreaROI=sum(ROI(:)); %Obtiene el área de la región de interés
         dif=(1/(Rec_min(3)*Rec_min(4)))*abs(Areas_ord(a)-AreaROI);
         diferenciasAreas(a,i-1)=dif;
         unionMask = ImFrames{1}|ImFrames{i};
         andMask = ImFrames{1}&ImFrames{i};
         DICE(1,i-1) = nnz(andMask) / nnz(unionMask);
     end
 end

 for i=1:1:(numFrames-1)
     diferenciasAreas(1,i) = mean(diferenciasAreas(:,i));
 end
 %multiplicando las dos metricas
 Kps = diferenciasAreas(1,:).*(1-DICE(1,:));
 %filtrado
 order = 3;
 framelen = 11;
 Kps = sgolayfilt(Kps,order,framelen);

 figure()
 x = [1:1:numFrames-1];
 plot(x,abs(Kps),'b',LineWidth=1.9)
 title('Factor de caracterizacion de los estados modales 5 areas')
 grid on
 xlabel('Frames')
 ylabel('Factor de diferencia')

 dif_mean = zeros (1,round(numFrames/vidObj.FrameRate));
 for i=1:1:round(numFrames/vidObj.FrameRate) %se comienza en dos para no usar el primer frame usado anteriormente
     fin = 25*i;
     if i == 1
         dif_mean(i) = mean(diferenciasAreas(1,1:fin));
     elseif i~=30
         dif_mean(i) = mean(diferenciasAreas(1,ini:fin));
     else
         dif_mean(i) = mean(diferenciasAreas(1,ini:fin-1));
     end
     ini = fin;
 end

 figure()
 x = [1:1:30];
 plot(x,dif_mean,'r',LineWidth=1.9);
 title('Factor de caracterizacion de los estados modales 5 areas')
 grid on
 xlabel('Tiempo en segundos')
 ylabel('Factor de diferencia')