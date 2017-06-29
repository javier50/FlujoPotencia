clc
clear all
close all
fprintf('\n                               "Flujo dinámico de potencia"\n')
fprintf('                               "javier50"\n')
%% CARGA DATOS A PARTIR DE UN ARCHIVO DE EXCEL
% NOMBRE DEL ARCHIVO
[LIBRO DIR]=uigetfile('*.xls','Archivo excel');
% LETRERO DE ESPERA
wai= waitbar(0,'Por favor espere. Los datos se están cargando. Gracias');
Nodos=xlsread(LIBRO,'Nodos');
waitbar(0.5,wai)
%Datos líneas
Lineas=xlsread(LIBRO,'Lineas');
aLineas=Lineas;
waitbar(1,wai)
close(wai)
%% PARÁMETROS 
nombre_de_la_red ='IEEE 39 NODOS'; % Nombre de la red con la cual se está trabajando para nombre en pantalla
% Tolerancia en pu
tol=0.01;
% Número máximo de iteraciones permitido
itmax=200;
baseMVA=100;
graus_to_rad=pi/180;
rad_to_graus=180/pi;
tic;
[nb,columnas]=size(Nodos); % Número total de nodos
nuevb=nb;
[nr,columnas]=size(Lineas); % Número total de líneas
%% Obtiene valores individuales de la hoja de Excel 
for k=1:nb
    numext(k)=Nodos(k,1);
    tipo(k)=Nodos(k,2);
    v(k)=Nodos(k,3);
    ang(k)=Nodos(k,4)*graus_to_rad;
    pg(k)=Nodos(k,5)/baseMVA;
    qg(k)=Nodos(k,6)/baseMVA;
    pc(k)=(Nodos(k,7)/baseMVA);
    qc(k)=(Nodos(k,8)/baseMVA);
    bshk(k)=Nodos(k,9)/baseMVA;
    pnom(k)=pg(k)-(pc(k));
    qnom(k)=qg(k)-(qc(k));
    numint(Nodos(k,1))=k;
end
for l=1:nr
        de(l)=numint(Lineas(l,1));
        para(l)=numint(Lineas(l,2));
        r(l)=Lineas(l,3);
        x(l)=Lineas(l,4);
        bshl(l)=Lineas(l,5)/2.0;
        sobre(l)=Lineas(l,8);
        tap(l)=Lineas(l,6);
        if tap(l)==0
            tap(l)=1;
        end
        fi(l)=-Lineas(l,7)*graus_to_rad;
end
%% Construcción de la matriz de admitancia nodal
for k=1:nb
    Y(k,k)=i*bshk(k);
end
for l=1:nr
    k=de(l);
    m=para(l);
    y(l)=1/(r(l)+i*x(l));
    akk(l)=1/(tap(l)*tap(l));
    amm(l)=1.0;
    akm(l)=1/tap(l);
    Y(k,k)=Y(k,k)+akk(l)*y(l)+i*bshl(l);
    Y(m,m)=Y(m,m)+amm(l)*y(l)+i*bshl(l);
    Y(k,m)=Y(k,m)-akm(l)*y(l);
    Y(m,k)=Y(m,k)-akm(l)*y(l);
end
%% Conductancia y susceptancia
G=real(Y);
B=imag(Y);
%% MATRIZ YBUS FINAL
% disp('              --------------------------------------------------------------')
% fprintf('\n                               MATRIZ DE ADMITANCIA NODAL YBUS\n')
% disp('              --------------------------------------------------------------')
Ybus=Y;
Y;
%% MATRIZ ZBUS FINAL
% disp('              --------------------------------------------------------------')
% fprintf('\n                               MATRIZ DE IMPEDANCIA NODAL ZBUS\n')
% disp('              --------------------------------------------------------------')
% Construcción de la matriz Zbus
for z=1:nb
    if Nodos(z,2)==3
        Z=Ybus;
        Z(:,z)=[]; % Elimina la columna donde está el nodo Slack
        Z(z,:)=[]; % Elimina la fila donde está el nodo Slack
        Zbus=inv(Z);
    end
end
Zbus;
Ypolar_Mag=abs(Y);
Ypolar_ang=180*angle(Y)/pi;
%% Condiciones iniciales
for k=1:nb
    if tipo(k)~=3
        ang(k)=0.0;
        if tipo(k)<2
            v(k)=1.0;
        end
    end
end
%% MÉTODO ITERATIVO DE NEWTON-RHAPSON
iter=0;
maxDP=100;
maxDQ=100;
while (1)
    if((abs(maxDP)>tol)||(abs(maxDQ)>tol))
    for k=1:nb
        pcalc(k)=G(k,k)*v(k)*v(k);
        qcalc(k)=-B(k,k)*v(k)*v(k);
    end
    for l=1:nr
        k=de(l);
        m=para(l);
        ab=ang(k)-ang(m)+fi(l);
        gkm=akm(l)*real(y(l));
        bkm=akm(l)*imag(y(l));
        pcalc(k)=pcalc(k)+v(k)*v(m)*(-gkm*cos(ab)-bkm*sin(ab));
        pcalc(m)=pcalc(m)+v(k)*v(m)*(-gkm*cos(ab)+bkm*sin(ab));
        qcalc(k)=qcalc(k)+v(k)*v(m)*(-gkm*sin(ab)+bkm*cos(ab));
        qcalc(m)=qcalc(m)-v(k)*v(m)*(-gkm*sin(ab)-bkm*cos(ab));
    end
    DP=zeros(nb,1);
    DQ=zeros(nb,1);
    maxDP=0;
    maxDQ=0;
    busDP=0;
    busDQ=0;
    pcalc;
    qcalc;
    for k=1:nb
        if tipo(k)~=3
            DP(k)=pnom(k)-pcalc(k);
            if abs(DP(k))>abs(maxDP)
                maxDP=DP(k);
                busDP=numext(k);
            end
        end
        if tipo(k)<=1
            DQ(k)=qnom(k)-qcalc(k);
            if abs(DQ(k))>abs(maxDQ)
                maxDQ=DQ(k);
                busDQ=numext(k);
            end
        end
    end
    DP;
    DQ;
%% Construcción del Jacobiano
    H=spalloc(nb,nb,nb+2*nr);M=H;N=H;L=H;
    for k=1:nb
        H(k,k)=-qcalc(k)-v(k)*v(k)*B(k,k);
        N(k,k)=(pcalc(k)+v(k)*v(k)*G(k,k))/v(k);
        M(k,k)=pcalc(k)-v(k)*v(k)*G(k,k);
        L(k,k)=(qcalc(k)-v(k)*v(k)*B(k,k))/v(k);
        if tipo(k)==3
            H(k,k)=10^10;
        end
        if tipo(k)>=2
            L(k,k)=10^10;
        end
    end
    for l=1:nr
        k=de(l);
        m=para(l);
        ab=ang(k)-ang(m)+fi(l);
        H(k,m)=v(k)*v(m)*(G(k,m)*sin(ab)-B(k,m)*cos(ab));
        H(m,k)=v(k)*v(m)*(-G(k,m)*sin(ab)-B(k,m)*cos(ab));
        N(k,m)=v(k)*(G(k,m)*cos(ab)+B(k,m)*sin(ab));
        N(m,k)=v(m)*(G(k,m)*cos(ab)-B(k,m)*sin(ab));
        M(k,m)=-v(k)*v(m)*(G(k,m)*cos(ab)+B(k,m)*sin(ab));
        M(m,k)=-v(k)*v(m)*(G(k,m)*cos(ab)-B(k,m)*sin(ab));
        L(k,m)=v(k)*(G(k,m)*sin(ab)-B(k,m)*cos(ab));
        L(m,k)=-v(m)*(G(k,m)*sin(ab)+B(k,m)*cos(ab));
    end
%% Jacobiano
    J=[H N;M L];
    DS=[DP;DQ];
    DV=J\DS;
%% Actualiza variables
    for k=1:nb
        ang(k)=ang(k)+DV(k);
        v(k)=v(k)+DV(k+nb);
    end
    iter=iter+1;
    else
        break
    end
end 
%% Cálculo de flujo de potencia
for l=1:nr
    k=de(l);
    m=para(l);
    gkm=real(y(l));
    bkm=imag(y(l));
    ab=ang(k)-ang(m)+fi(l);
    vkm=v(k)*v(m);
    pkm(l)=akk(l)*v(k)*v(k)*gkm-akm(l)*vkm*(gkm*cos(ab)+bkm*sin(ab));
    pmk(l)=amm(l)*v(m)*v(m)*gkm-akm(l)*vkm*(gkm*cos(ab)-bkm*sin(ab));
    qkm(l)=-akk(l)*v(k)*v(k)*(bkm+bshl(l))+akm(l)*vkm*(bkm*cos(ab)-gkm*sin(ab));
    qmk(l)=-amm(l)*v(m)*v(m)*(bkm+bshl(l))+akm(l)*vkm*(bkm*cos(ab)+gkm*sin(ab));
    skm(l)=sqrt((pkm(l)^2)+(qkm(l)^2))*baseMVA;
    smk(l)=sqrt((pmk(l)^2)+(qmk(l)^2))*baseMVA;
    if skm(l)>=smk(l)
        saprox(l)=skm(l);
    else
        smk(l)>=skm(l);
        saprox(l)= smk(l);
    end
    pperdas(l)=pkm(l)+pmk(l);
    qperdas(l)=qkm(l)+qmk(l);
end
%% Muestra resultados en pantalla
fprintf(' Red: %s\n',nombre_de_la_red)
fprintf(' Converge en %d iteraciones\n\n',iter)
fprintf('                              CONDICIONES DE LA RED\n\n')
disp('----------------------------------------------------------------------------------------------')
fprintf('    Nodo       Tipo       Mag(p.u)    Ángulo(G)       P(MW)       Q(MVAR)        Qsh(MVAR)\n')
disp('----------------------------------------------------------------------------------------------')
for k=1:nb
    fprintf('%7d      %4d      %9.3f      %6.2f      %9.2f      %7.2f      %9.2f \n',numext(k),tipo(k),v(k),(ang(k)*rad_to_graus),baseMVA*(pcalc(k)+pc(k)),baseMVA*(qcalc(k)+qc(k)),baseMVA*bshk(k)*v(k)^2)
end
disp('----------------------------------------------------------------------------------------------')
fprintf('\n                                             RESULTADO DE FLUJOS DE POTENCIA\n\n')
disp('---------------------------------------------------------------------------------------------------------------------------------------')
fprintf('    Linea       De     Para      Pkm(MW)   Qkm(MVAR)     Pmk(MW)   Qmk(MVAR)     Ploss(MW) Qloss(MVAR) Skm(MVA)  Smk(MVA)    S(MVA)\n')
disp('---------------------------------------------------------------------------------------------------------------------------------------')
for l=1:nr
    fprintf('   %4d    %7d    %4d    %9.2f    %7.2f    %9.2f    %7.2f    %9.2f    %7.2f    %7.2f    %7.2f    %7.2f\n',l,de(l),para(l),baseMVA*pkm(l),baseMVA*qkm(l),baseMVA*pmk(l),baseMVA*qmk(l),baseMVA*pperdas(l),baseMVA*qperdas(l),skm(l),smk(l),saprox(l))
end
disp('---------------------------------------------------------------------------------------------------------------------------------------')
fprintf('\n               SOBRECARGA EN LAS LÍNEAS\n\n')
disp('  ----------------------------------------------------')
fprintf('        Linea         De       Para     Porcentaje  \n')
disp('  ----------------------------------------------------')
to=0;
for l=1:nr
    if saprox(l)>sobre(l)
        cargabl(l)=((saprox(l)/sobre(l))*100)-100;
        fprintf('       %4d      %7d      %4d      %4.5f\n',l,de(l),para(l),cargabl(l));
        to(l)=l;
    end
end
disp('  ----------------------------------------------------')
to(to==0)=[];
tosobre=length(to);
fprintf('        Un total de %d líneas tienen sobrecarga\n',tosobre);
disp('  ----------------------------------------------------')
fprintf('\n     EXCESO DE NIVEL DE TENSIÓN\n\n')
disp('     ---------------------------')
fprintf('        Nodo       Porcentaje\n')
disp('     ---------------------------')
to1=0;
for l=1:nb
    if v(l)>Nodos(l,10)
        ex(l)=((v(l)/Nodos(l,10))*100)-100;
        fprintf('       %4d        %4.5f\n',l,ex(l));
        to1(l)=l;
    end
end
disp('     ---------------------------')
to1(to1==0)=[];
exsobre=length(to1);
disp('-------------------------------------------------')
fprintf('\n  Un total de %d nodos tienen exceso de tensión\n',exsobre);
disp('-------------------------------------------------')
fprintf('\n        BAJO NIVEL DE TENSIÓN\n\n')
disp('     ---------------------------')
fprintf('        Nodo       Porcentaje\n')
disp('     ---------------------------')
ato1=0;
for l=1:nb 
    if v(l)<Nodos(l,11)
        ba(l)=((v(l)/Nodos(l,11))*100)-100;
        fprintf('       %4d        %4.5f   \n',l,abs(ba(l)));
        ato1(l)=l;
    end
end
disp('     ---------------------------')
ato1(ato1==0)=[];
basobre=length(ato1);
disp('-------------------------------------------------')
fprintf('\n  Un total de %d nodos tienen bajo nivel de tensión\n',basobre);
disp('-------------------------------------------------')