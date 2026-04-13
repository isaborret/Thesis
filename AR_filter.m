% AUTOREGRESSIVE, HIGH-PASS, AND LOW-PASS FILTER PROGRAM

% sig is a multivariate file with can channels
% can is the channel to be filtered
% p is the order of the AR filter

% AR filter separates drifts and fast rythms 
% by estimating the slow trend (fib, the low-pass output) 
% and subtracting it from the signal (fia, the high-pass output). 
% You then feed fia into the SPWVD to extract HF power, because the slow drift would otherwise distort the time-frequency representation.


function [fia,fib]=AR_filter(sig,can,p)

s=sig(:,can); % extracts column can from the input matrix. Since you call it with can = 1 and pass the RR series as a column vector, this just gives you the RR series itself.

s1=size(s,1);%d Gets the number of samples in the RR series

% Pre-allocate output vectors with zeros before filling them in the loops
usc = zeros(s1,1);
fib = zeros(s1,1);

% "andata" 
% First-order recursive low-pass filter 
% Each output sample is a weighted average of the previous output (p = 0.94 weight) and the current input (1-p = 0.06 weight).
usc(1)= s(1);
for i = 2 : s1
   usc(i)= p*usc(i-1)+(1-p)*s(i);
end 

% "ritorno"
% The same filter is applied again but running backwards through the signal. This is the zero-phase trick — running forward then backward cancels out any phase delay introduced by the filter
fib(s1)=usc(s1);
for i = (s1-1):-1:1
     fib(i)= p*fib(i+1)+(1-p)*usc(i);
end
   
%passa alto
% high pass filter
fia = s-fib+mean(s);
% Subtracts the slow trend from the original, leaving only fast fluctuations
% Adding mean(s) back preserves the DC level so the series stays centred around its original mean rather than around zero.


% % crea un vettore contenente segnale filtrato e da filtrare 
% % e li visualizza su un unico grafico
% figure(2);
% subplot(2,1,1);
% plot ([s fib]),title('originale+passa basso');
% zoom xon;
% subplot(2,1,2);
% plot ([s fia]),title('originale+passa alto');
% zoom xon;