% CODED BY: PUNEET DHEER (RF)
% Date: 01-August-2017
% MAIN_FUNCTION
% MAGNITUDE SQUARED COHERENCE (BAND_AVERAGED)
%%
function [cohh] = BA_msc(co1,co2,Ws,Lp,Fre_range)

Lw=1;
Z=Ws;

X1=co1;
X2=co2;

windows=ceil((length(X1)-Ws+1)/Lp); 
 
fft_window=hanning(Ws); %to reduce spectral leakage

for i=1:windows
    
%     nffts=2^nextpow2(length(A)); 
    A = fft(X1(Lw:Z).*fft_window');%fft of X
    B = fft(X2(Lw:Z).*fft_window');%fft of Y
    
%     A = fft(X1(Lw:Z));%fft of X
%     B = fft(X2(Lw:Z));%fft of Y
    
    A = A(1:floor(length(A)/2)+1); %selecting only first falf portion (1:length(x)/2) of the fft because the other 2nd half remains the same (reject the negative frequency)
    B = B(1:floor(length(B)/2)+1);
    

%%  AVERAGE across all frequency
%1
%     XY = (abs(mean(A.*conj(B))))^2;
%     XX = mean((A.*conj(A)));
%     YY = mean((B.*conj(B)));

%2cmtm
%     XY = A.*conj(B); %in complex form
%     XX = A.*conj(A); %already abs
%     YY = B.*conj(B); %already abs
%     Cxy = XY ./ sqrt((XX.*YY));
%     windowed_coh(i) = abs(mean(Cxy(Fre_range)));

%3    
    % AVERAGE across frequency range 
    XY = (abs(mean(A(Fre_range+1).*conj(B(Fre_range+1)))))^2;
    XX = mean((A(Fre_range+1).*conj(A(Fre_range+1))));
    YY = mean((B(Fre_range+1).*conj(B(Fre_range+1))));
    

%% OR
% Mathematically
%     XY = abs(A).*abs(B).*exp(1i*(angle(A)-angle(B)));
%     XY = (abs((XY))).^2;
% 
%     XX = abs(A).*abs(A);
%     XX = (XX);    
% 
%     YY = abs(B).*abs(B);
%     YY = (YY);
%%     
     
    coh = XY./(XX.*YY);

    windowed_coh(i) = (coh);
    
    Lw = Lw+Lp;
    Z = Z+Lp;      

%%
cohh = windowed_coh;

end