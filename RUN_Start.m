clc
clear all
close all
addpath('Utils');

%OFDM Modulator
Ncarrier=1024;
NFFT=1024;
M=4;
Nbits=Ncarrier*10*100;
Fnull=0 ;
GI = 246/1000;
OSR=6;
binaryData=randi([0, 1], Nbits, 1);
[txSignal]=modul_ofdm2013(binaryData,M,Ncarrier,NFFT,GI,Fnull, OSR);

txSignal=txSignal./max(abs(txSignal));
%TX Imbalnce
if(0)
    Ki=0;
    Kq=0;
    pha=0;

    txSignal=IQ_Imbalance(txSignal, Ki, Kq, pha);
    txSignal=txSignal./max(abs(txSignal));
end

%---------------------------------------------------------------------------

%memoryless USRP RACOM PA
PA_coef=[0.845654964236630 + 1.49110374573613i;-2.76262037337481 + 1.88247060093166i;17.0635003318257 - 11.9190158113625i;-51.8235490530484 + 29.6408107759511i;81.3986654536564 - 40.1728308743692i;-64.9809973172151 + 27.2732197670390i;20.6182249621985 - 7.32252535080311i];  Kpa=7; Qpa=0;

rxSignal=PA_Model2(txSignal,PA_coef, Kpa, Qpa);  % Calculate PA model output
rxSignal=(rxSignal./max(abs(rxSignal))).';
%---------------------------------------------------------------------------
if(1)
%RX Imbalance
Ki=0;
Kq=0;
pha=0;

rxSignal=IQ_Imbalance(rxSignal, Ki, Kq, pha);
rxSignal=(rxSignal./max(abs(rxSignal)));

%First we model the PA (to have measurements on right-hand side
K=9;
Q=0;

txSignal=txSignal./max(abs(txSignal));
Y=DDR2_Matrix(txSignal,K,Q);

M=[real(Y) -imag(Y)];
x=real(rxSignal)./max(abs(real(rxSignal))); %normalize with respect to 1
b=inv((M')*M)*(M')*(x);
b2=b(1:(size(b,1)/2))+1i*b((size(b,1)/2)+1:end);
DPDoutput=DDR2_memory_polynomials(txSignal,K,Q,b2);

%---------------------------------------------------------------------------
%DPD
K=7;
Q=0;
Y=DDR2_Matrix(DPDoutput,K,Q);
x=txSignal;
b=inv((Y')*Y)*(Y')*(x);
DPDoutput=DDR2_memory_polynomials(txSignal,K,Q,b);
%---------------------------------------------------------------------------
end

DPDoutput_PA=PA_Model2(DPDoutput,PA_coef, Kpa, Qpa);
%---------------------------------------------------------------------------


if(1)
plot(abs(txSignal), abs(rxSignal),'.', abs(txSignal), abs(DPDoutput),'.', abs(txSignal), abs(DPDoutput_PA),'.')
%figure
%plot(abs(txSignal), angle(txSignal)-angle(rxSignal),'.', abs(txSignal), angle(txSignal)-angle(DPDoutput_PA.'),'.')
end


if(1)
figure
plot(10*log(abs(fftshift(fft(rxSignal)))));
hold on
plot(10*log(abs(fftshift(fft(DPDoutput_PA)))));

end





