%% processamento do tempo

T=[0.14,0.05,0.1,0.17,0.14,0.36,0.2,0.2,0.24,0.23,0.13,0.18,0.07,0.23,0.23,0.11,0.09,0.2,0.13,.15,0.17,0.17,0.14];
% T sendo o tempo de reaçao de todos os alunos


mediaT=mean(T); % media do tempo de todos os alunos

varianceT=var(T); %variancia

desvT=std(T); %desvio padrão do tempo

desvMean=desvT/sqrt(23); %desvio padrão da média do tempo

TEMPO=cat(1,mediaT,varianceT,desvT,desvMean);

clear mediaT varianceT desvT desvMean

%% Massa do cobre

%mC1 a peça 1 que não possui marcação
%mC2 a peça 2 que possui marcaçao preta

mC1=1064.58;
mC2=1081.03;

m=mC1+mC2;

m=m/1000;

clear mC1 mC2

N=30;

%% mola 1
tMola1=[18.23,17.80,18.37]; %tempo de oscilação da mola 1

mediaMola1=mean(tMola1); %media da mola 1

varMola1=var(tMola1); %variancia da mola 1

desvMola1=std(tMola1); %desvio padrão da mola 1

desvMeanMola1=desvMola1/sqrt(3); %desvio padrão da media da mola 1

mola1=[mediaMola1,varMola1,desvMola1,desvMeanMola1]; %vetor com toda a informaçao da mola 1

clear mediaMola1 varMola1 desvMola1 desvMeanMola1
%% mola 2
tMola2=[18.14,17.87,18.30]; %tempo de oscilação da mola 1

mediaMola2=mean(tMola2); %media da mola 2

varMola2=var(tMola2); %variancia da mola 2

desvMola2=std(tMola2); %desvio padrão da mola 2

desvMeanMola2=desvMola2/sqrt(3); %desvio padrão da media da mola 2

mola2=[mediaMola2,varMola2,desvMola2,desvMeanMola2]; %vetor com toda a informaçao da mola 2

clear mediaMola2 varMola2 desvMola2 desvMeanMola2

%% mola serie
tMolaSerie=[25.48,26.3,26.95]; %tempo de oscilação da mola serie

mediaMolaSerie=mean(tMolaSerie); %media da mola serie

varMolaSerie=var(tMolaSerie); %variancia da mola serie

desvMolaSerie=std(tMolaSerie); %desvio padrão da mola serie

desvMeanMolaSerie=desvMolaSerie/sqrt(3); %desvio padrão da media da mola serie

molaSerie=[mediaMolaSerie,varMolaSerie,desvMolaSerie,desvMeanMolaSerie]; %vetor com toda a informaçao da mola serie

clear mediaMolaSerie varMolaSerie desvMolaSerie desvMeanMolaSerie


%% mola paralelo
tMolaParalelo=[12.91,12.8,13.21]; %tempo de oscilação da mola paralelo

mediaMolaParalelo=mean(tMolaParalelo); %media da mola paralelo

varMolaParalelo=var(tMolaParalelo); %variancia da mola paralelo

desvMolaParalelo=std(tMolaParalelo); %desvio padrão da mola paralelo

desvMeanMolaParalelo=desvMolaParalelo/sqrt(3); %desvio padrão da media da mola paralelo

molaParalelo=[mediaMolaParalelo,varMolaParalelo,desvMolaParalelo,desvMeanMolaParalelo]; %vetor com toda a informaçao da mola paralelo

clear mediaMolaParalelo varMolaParalelo desvMolaParalelo desvMeanMolaParalelo

%% 

molas=cat(1,mola1,mola2,molaSerie,molaParalelo); %concatenar todos os vetores em uma matriz

clear mola1 mola2 molaSerie molaParalelo

%% derivadas parciais para encontrar a formula da incerteza de k

syms M tn %definindo massa e tempo como simbolico

kk=4*pi*M*30^2/tn^2; %definindo a funçao k como simbolica

ddm=diff(kk,M); %derivando k em relaçao a massa

ddtn=diff(kk,tn); %derivando k em relaçao ao tempo

display(ddm); %mostrar ddm de uma maneira visualizavel para utilizar em seguida

display(ddtn); %mostrar ddtn de uma maneira visualizavel para utilizar em seguida


%% calculo de k e incerteza

k=NaN(4,2); %criar uma matriz vazia para receber o valor de k e a respectiva incerteza

for i=1:4
    k(i,1) = (4*(pi^2)*m*(N^2))/(molas(i,1)^2); %calcula o k para todas as médias de calculadas anteriormente
    k(i,2) = (((-(7200*pi*m)/molas(i,1)^3)*0.02)^2)+((((3600*pi)/molas(i,1)^2)*(0.01))^2)^(1/2); %calculo da incerteza baseado nas derivadas calculadas anteriormente
end

display(k); %mostrar o valor de k (coluna 1) e sua respectiva incerteza (coluna 2)
%Linha 1: mola1
%Linha 2: mola2
%Linha 3: mola serie
%Linha 4: mola paralelo

%% cálculo do erro normalizado

Va= k(1,1);
Vb= k(2,1);
Uva= k(1,2);
Uvb= k(2,2);

erroNormalizado= abs((Va-Vb)/((Uva^2+Uvb^2)^(1/2)));


%% tratamento do g

tempo10oscilacoes=[28.39,28.74,28.83,28.84];

mean10o=mean(tempo10oscilacoes);

var10o=var(tempo10oscilacoes);

desv10o=std(tempo10oscilacoes);

desvMean10o=desv10o/sqrt(3);


syms L t
g=L*(2*pi/t)^2;
diff(g,L)
diff(g,t)

L=2.1;


G=4*pi^2*L/((mean10o/10)^(2));
incertezaG=((((4*pi^2)/t^2)*0.01)^2) + (((-(8*L*pi^2)/t^3)*0.02)^2);

