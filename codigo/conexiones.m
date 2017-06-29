clc
clear all
close all
fprintf('\n                               "Flujo dinámico de potencia"')
fprintf('                               "javier50"\n')

%% CARGA DATOS A PARTIR DE UN ARCHIVO DE EXCEL
[LIBRO DIR]=uigetfile('*.xls','Archivo excel');% NOMBRE DEL ARCHIVO
wai= waitbar(0,'Por favor espere. Los datos se están cargando. Gracias');% LETRERO DE ESPERA
Nodos=xlsread(LIBRO,'Nodos');%Extrae Nodos de libro
waitbar(0.5,wai)% Letrero 50%
Lineas=xlsread(LIBRO,'Lineas');%Extrae Lineas de libro
waitbar(1,wai);% Letrero 100%
close(wai);%Cierra Letrero

%Tipos de Nodo
tipo1=3;
tipo2=2;
tipo3=0;

%Numero de columna de Nodo
colIndice=1;
colTipoNodo=2;
colPg=5;
colQg=6;
colPc=7;
colQc=8;

%Numero de columna de Linea
colDe=2;
colPara=3;
colSobreCarga=9;
colConectado=10;

%% Obtener Lineas conectadas
[nrLineas,nc]=size(Lineas);
listaLineas=[];

for k=1:nrLineas
   if Lineas(k,colConectado) == 1
       listaLineas=cat(1,listaLineas,Lineas(k,:));
   end
end
Lineas = listaLineas;


%% Encontrar nodos con conexion a otros nodos
[nrLineas,ncLineas]=size(Lineas);
[nrNodos,ncNodos]=size(Nodos);
nodosTemp=Nodos;
lineasTemp=Lineas;
nIslas=0;
hayNodosValidos=true;

%Termina ciclo hasta que todas las lineas hallan sido agregadas a una isla
while hayNodosValidos 
    nodosIslaTemp=[];
    lineasIslaTemp=[];
    agregoLinea=false;
    seAgregoLinea = true;
    nIslas=nIslas+1;
    
    %termina ciclo hasta que una isla se haya completado        
    while seAgregoLinea
        [nrNodosIsla,nc]=size(nodosIslaTemp);
        %itera todas las lineas en la isla para agregarla a la isla
        for i=1:nrLineas
            %Verifica si es el primer elemento
            agregarLinea=false;
            esPrimerElemento=false;
            indiceNodoAgregar=0;
            
            %agrega primer nodo si es el primero en la lista
            if nrNodosIsla==0
                for iNodos=1:nrNodos
                    %fprintf('\nindiceN: %d, tipoNodo:%d',Nodos(iNodos, colIndice), Nodos(iNodos, colTipoNodo));
                    if Nodos(iNodos, colTipoNodo) == tipo1 || Nodos(iNodos, colTipoNodo) == tipo2
                        nodosIslaTemp=Nodos(iNodos,:);
                        nodosTemp(iNodos,:)=[];
                        break;
                    end
                end
            end
            
            %Valida si la conexion coincide con algun nodo de la lista nodosIslaTemp
            [nrNodosIsla,nc]=size(nodosIslaTemp);
            %fprintf('\n[%d,%d]',size(nodosIslaTemp));
            for j=1:nrNodosIsla
            	%fprintf('\nindiceN:%d, indiceL:%d, de:%d, para:%d',nodosIslaTemp(j,colIndice), Lineas(i,colIndice), Lineas(i,colDe), Lineas(i,colPara));
                if nodosIslaTemp(j,colIndice)==Lineas(i,colDe) %si se encuentra indice se rompe el ciclo
                	agregarLinea=true;
                    indiceNodoAgregar=Lineas(i,colPara);%Posiblemente aun no se encuentre en en lista nodosIslaTemp
                    break;
                end

                if nodosIslaTemp(j,colIndice)==Lineas(i,colPara) %si se encuentra indice se rompe el ciclo
                	agregarLinea=true;
                    indiceNodoAgregar=Lineas(i,colDe);%Posiblemente aun no se encuentre en en lista nodosIslaTemp
                    break;
                end
            end
            
            %agrega Linea a lista de lineasIsla y la remueve de lineasTemp
            if agregarLinea
                [nrLineasTemp,ncLineasTemp]=size(lineasTemp);
                for iBorrarLinea=1:nrLineasTemp
                    if lineasTemp(iBorrarLinea,colIndice)== Lineas(i,colIndice)
                        lineasTemp(iBorrarLinea,:)=[];
                        break; 
                    end    
                end
                lineasIslaTemp=cat(1,lineasIslaTemp, Lineas(i,:));
                seAgregoLinea=true;
            else
                seAgregoLinea=false;
            end
            
            %Valida que aun no se encuentre en lista nodosIslaTemp
            [nrNodosIsla,nc]=size(nodosIslaTemp);
            agregarNodo=true;
            for j=1:nrNodosIsla
                if indiceNodoAgregar==nodosIslaTemp(j,colIndice) %si se encuentra indice se rompe el ciclo
                	agregarNodo=false;
                    break;
                end
            end
            
            %Agrega Nodo si aun no se encuentra en la lista de nodosIslaTemp
            if agregarNodo
                [nrNodos, nc]=size(Nodos);
                for iBorrarNodo=1:nrNodos
                    if indiceNodoAgregar == Nodos(iBorrarNodo,colIndice)
                        [nrNodosTemp,nc]=size(nodosTemp);
                        %busca nodo para eliminarlo
                        for jBorrarNodo=1:nrNodosTemp
                            if indiceNodoAgregar == nodosTemp(jBorrarNodo,colIndice)
                                nodosTemp(jBorrarNodo,:)=[];
                                break;
                            end
                        end

                        %Agrega Nodo a lista nodosIslaTemp y rompe ciclo
                        nodosIslaTemp=cat(1,nodosIslaTemp, Nodos(iBorrarNodo,:));
                        break;
                    end
                end
            end
        end
        
        
        Lineas=lineasTemp;
        [nrLineas,ncLineas]=size(Lineas);
    end%cierra while seAgregoLinea
    
    %Elimina primer columna de Lineas para poder entrar al metodo de
    %newthon rhapson    
    lineasIslaTemp(:,colIndice)=[];
    
    %% Imprime lista de Nodos de la Isla a iterar
    fprintf('\n\n------------------------------------------------------------------------------------------------------');
    fprintf('\n Nodos de Isla:%d\n', nIslas);
    [nrC, ncC]=size(nodosIslaTemp);
    fprintf('\n \tNodo \tTipo\t\tVolt\t\tAng \tPg\t \tQg\t \tPc\t \tQc\t\tBshk\t\tVmax\t\t\tVim\n');
    for iC=1:nrC
        for jC=1:ncC
            if jC==3 || jC==10 || jC==11
                fprintf('| \t%f\t',nodosIslaTemp(iC,jC));
            else
                fprintf('| \t%d\t',nodosIslaTemp(iC,jC));
            end
        end
        fprintf('|\n');
    end
    
    
    %% Imprime lista de Lineas de la Isla a iterar
    fprintf('\n Lineas de Isla:%d\n',nIslas);
    [nrL, ncL]=size(lineasIslaTemp);
    fprintf('\n \tDe\t \tPara\t \tr\t \t\t\tx\t\t\t\tBshl\t\tTap\t \tFi\tSobCarga\tConexion\n');
    for iL=1:nrL
        for jL=1:ncL
            if jL==3 || jL==4 || jL==5
                fprintf('| \t%f\t',lineasIslaTemp(iL,jL));
            else
                fprintf('| \t%d\t',lineasIslaTemp(iL,jL));
            end
        end
        fprintf('\n');
    end
    
    
    %% Valida que potencia isla pueda pasar al metodo de newtonRapson
    [nrC, ncC]=size(nodosIslaTemp);
    indiceNodoTip2=0;
    pgMayor=0;
    validacionIsla=false;    
    nodostipo2Lista = [];
    nodostipo3Lista = [];
    
    for iC=1:nrC
        switch nodosIslaTemp(iC,colTipoNodo)
            case tipo1
                validacionIsla = true;
                break;
            case tipo2
                nodostipo2Lista=cat(1,nodostipo2Lista, nodosIslaTemp(iC,:));
                
                %Obtiene el nodo con mayo Pg por si se necesita
                if nodosIslaTemp(iC,colPg) > pgMayor 
                    pgMayor = nodosIslaTemp(iC,colPg);
                    indiceNodoTip2=iC;
                end
                
                break;
            case tipo3
                nodostipo3Lista=cat(1,nodostipo3Lista, nodosIslaTemp(iC,:));
                break;
        end
        
        if validacionIsla
            break;
        end
    end
    
    %Valida las potencias del la isla si no hay nodo tipo 3
    if ~validacionIsla
        [nrN2, nc]=size(nodostipo2Lista);
        [nrN3, nc]=size(nodostipo3Lista);
        
        if nrN2 == 0
            fprintf('\nLa isla no tinene nodos Tipo 3 ni tipo 2');
        elseif nrN2 == 1
            pcTotalNodoT3=0;
            %Suma las Pg de los nodos tipo 3
            for iC=1: nrN3
                pcTotalNodoT3=pcTotalNodoT3+nodostipo3Lista(iC,colPc);
            end
            
            %Compara Pc de nodos tipo3 con Pg nodo tipo2
            if pcTotalNodoT3 >= nodostipo2Lista(1, colPg)
                nodosIslaTemp(indiceNodoTip2,colTipoNodo) = tipo1;
                validacionIsla = true;
            else
                fprintf('\nLa Pc de los nodos Tipo 3 es mayor a la Pg del nodo Tipo 2');
            end
                        
        elseif nrN2 > 1
            qcTotalNodoT3=0;
            qgTotalNodoT2=0;
            
            %Suma los Qc de los nodos tipo3
            for iC=1: nrN3
                qcTotalNodoT3=pcTotalNodoT3+nodostipo3Lista(iC,colQc);
            end
            
            for iC=1: nrN2
                qgTotalNodoT2=pcTotalNodoT3+nodostipo3Lista(iC,colQg);
            end
            
            %Compara Pc de nodos tipo3 con Pg nodo tipo2
            if qcTotalNodoT3 >= qgTotalNodoT2
                nodosIslaTemp(indiceNodoTip2,colTipoNodo) = tipo1;
                validacionIsla = true;
            else
                fprintf('\nLa Qc de los nodos Tipo 3 es mayor a la Qg del los nodos Tipo 2');
            end            
        end
    end
    
    
    %% Se pasa Isla a metodo de Newthon Rhapson
    
    if validacionIsla
        try
           newthonRhapson(nodosIslaTemp, lineasIslaTemp);
        catch exception
           fprintf('\nOcurrio un problema en el metodo Newton Rhapson');
        end        
    else
        fprintf('\nLa isla no cumple con lo necesario para entrar al metodo Newton Rhapson');
    end
        
    
    %% Actualiza Datos para la siguiente iteracion
    Nodos=nodosTemp;
    [nrNodos,ncNodos]=size(Nodos);
    
    
    %% Valida si hay nodos y lineas para formar otra Isla
    
    %Valida si hay nodos validos para formar otra Isla
    indicesNodoValido=[];
    hayNodosValidos=false;
    for iC=1:nrNodos
        if Nodos(iC,colTipoNodo) == tipo1 || Nodos(iC,colTipoNodo) == tipo2
            indicesNodoValido=cat(2,indicesNodoValido,Nodos(iC,colIndice));
        end
    end
       
    %Valida si hay Lineas validos para formar otra Isla
    [nrNodosV,nc]=size(indicesNodoValido);
    if nrNodosV>0
        for iL=1:nrLineas
            for jL=1:nrNodosV
                if indicesNodoValido(jL) == Lineas(iL,colDe) || indicesNodoValido(jL) == Lineas(iL,colPara)
                    hayNodosValidos=true;
                    break;
                end
            end
            
            if hayNodosValidos
                break;
            end
        end
    end
    
    
    %% Imprime nodos que no se pudieron conectar a una isla
    if ~hayNodosValidos
        fprintf('\n\n-------------------------------------------------------------------------------------------------------------\n');
        fprintf('\nNodos que no se pudieron conectar a una Isla');
        fprintf('\n \tNodo \tTipo\t\tVolt\t\tAng \tPg\t \tQg\t \tPc\t \tQc\t\tBshk\t\tVmax\t\t\tVim\n');
        for iC=1:nrNodos
            for jC=1:ncNodos
                if jC==3 || jC==10 || jC==11
                    fprintf('| \t%f\t',Nodos(iC,jC));
                else
                    fprintf('| \t%d\t',Nodos(iC,jC));
                end
            end
            fprintf('|\n');
        end
        
        fprintf('\n\nLineas que no se pudieron conectar a una Isla');
        fprintf('\n \tDe\t \tPara\t \tr\t \t\t\tx\t\t\t\tBshl\t\tTap\t \tFi\tSobCarga\tConexion\n');
        for iL=1:nrLineas
            for jL=2:ncLineas
                if jL==4 || jL==5 || jL==6
                    fprintf('| \t%f\t',Lineas(iL,jL));
                else
                    fprintf('| \t%d\t',Lineas(iL,jL));
                end
            end
            fprintf('\n');
        end    
    end
end