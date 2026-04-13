%STIMA DELLO SPETTRO DI POTENZA CON IL METODO WEIGHTED COVARIANCE
%% This function estimates the power spectral density (PSD) of an EEG segment using the weighted covariance method.

function [totVar, pot, freq]=funct_CalcPSDmodified(signal,Bw,window,nfft,fc,plotFigure)
    x=signal';      % Transposes the input signal from a column vector to a row vector, which is the expected format for xcorr late
    totVar=var(x);  % Computes the total variance of the signal segment, which equals total power for a zero-mean signal
    
    m=round(1.273/Bw); %Computes the maximum lag m at which to truncate the autocorrelation, derived from the desired frequency resolution - wider bandwidth means shorter lag truncation
    if ~(m/2==fix(m/2)), m=m-1; end %ensures m is even, which is required for symmetric windowing

    x=x(:,1); % Ensures x is a single column
%     x = detrend(x);
    x=x-mean(x);
    [rx, lags]=xcorr(x,m,'biased'); %Estimates the biased autocorrelation of the signal up to lag m
    % = core of the weighted covariance method:
        % instead of computing the FFT directly on the signal, 
        % it first estimates how correlated the signal is with itself at different time lags, 
        % then transforms that into a power spectrum 
    nfft2=2*nfft;%Sets the FFT length to 2*512 = 1024 points => just makes the spectrum smoother


    eval(strcat('wi=',window,'(',int2str(2*m+1),');')); % Applies a Hamming window to the autocorrelation estimate
    rxw=rx.*wi;
    %% posso usare uno di questi metodi
     px=fft(rxw,nfft2); %stima veloce (in realtà sarebbe da usare la formula WC, ma il modulo è lo stesso)
%     px_1=fft(rxw,nfft2);
%     px_2=fft(rxw,nfft);
%     pt = fft(x,nfft2);
    px=(2/fc)*abs(px(1:nfft)); % Computes the FFT of the windowed autocorrelation
    
    %prendo il modulo, normalizzo per fnyquist, e prendo solo da zero a fnyquist %%errato
%     px1=abs(px_1(1:nfft)); %prendo il modulo, normalizzo per fnyquist, e prendo solo da zero a fnyquist
%     px2=abs(px_2);
%     pt1 = (1/length(x))*(abs(fft(x,nfft)).^2); 
%     pt2 = (1/length(x))*(abs(pt(1:nfft)).^2);
    f=(0:fc/(2*(nfft-1)):fc/2)'; %Constructs the frequency axis vector

    freq = f; % output: what the extraction script receives back and uses for band power integration
    pot = px; % output: what the extraction script receives back and uses for band power integration
      
    % LF_start = 0.04;
    % LF_end = 0.15;
    % HF_start = 0.15;
    % HF_end = 0.4;
    delta_start=0.5;
    delta_end=3;
    theta_start=3;
    theta_end=8;
    alpha_start=8;
    alpha_end=12;
    sigma_start=12;
    sigma_end= 16;
    beta_start=16;
    beta_end=25;


    
     %trova i valori di delta,theta,alpha,beta
        
     % posLF_freq = find(freq>LF_start & freq<=LF_end);
     % LF_freq = freq(posLF_freq);
     % LF_pot = pot(posLF_freq);
     % totLF_pot = mean(LF_pot) * (LF_end-LF_start);
     % 
     % posHF_freq = find(freq>HF_start & freq<=HF_end);
     % HF_freq = freq(posHF_freq);
     % HF_pot = pot(posHF_freq);
     % totHF_pot = mean(HF_pot) * (HF_end-HF_start);

     posdelta_freq = find(freq>delta_start & freq<=delta_end); % identifies which frequency bins fall within the band boundaries
     delta_freq = freq(posdelta_freq);
     delta_pot = pot(posdelta_freq);
     totdelta_pot = mean(delta_pot) * (delta_end-delta_start); % averages the PSD values across those bins
     % Multiplying by the bandwidth (delta_end - delta_start) converts from mean power density to total band power 
     % tot*_pot values are computed but never returned or used outside the function — they only feed into the optional plot

     postheta_freq = find(freq>theta_start & freq<=theta_end);
     theta_freq = freq(postheta_freq);
     theta_pot = pot(postheta_freq);
     tottheta_pot = mean(theta_pot) * (theta_end-theta_start);

     posalpha_freq = find(freq>alpha_start & freq<=alpha_end);
     alpha_freq = freq(posalpha_freq);
     alpha_pot = pot(posalpha_freq);
     totalpha_pot = mean(alpha_pot) * (alpha_end-alpha_start);

     posbeta_freq = find(freq>beta_start & freq<=beta_end);
     beta_freq = freq(posbeta_freq);
     beta_pot = pot(posbeta_freq);
     totbeta_pot = mean(beta_pot) * (beta_end-beta_start);

     possigma_freq = find(freq>sigma_start & freq<=sigma_end);
     sigma_freq = freq(possigma_freq);
     sigma_pot = pot(possigma_freq);
     totsigma_pot = mean(sigma_pot) * (sigma_end-sigma_start);
     
     if (plotFigure==1) % Only executes if you pass plotFigure=1
         figure("WindowState","maximized");
         h_spectrum=plot(freq,pot);

        %  hold on; h_LF=area(LF_freq, LF_pot,'FaceColor',[1.0000    0.6706    0.6902]);
        %  hold on; h_HF=area(HF_freq, HF_pot,'FaceColor',[0.3882    1.0000    0.8392]);
        %  set(gca, 'XTick', unique([LF_start LF_end HF_end, get(gca, 'XTick')]));
        %  xlabel('Frequency [Hz]'); ylabel('RRI Power'); %%set assi
        %  xlim([0 0.5]);
        %  title('PSD of RRI'); %%generico per lo script cond_subj
        %  lh = legend([h_spectrum h_LF h_HF], ...
        % 'Spettro', 'Potenza LF', 'Potenza HF');

         hold on; h_delta=area(delta_freq, delta_pot,'FaceColor',[1.0000    0.6706    0.6902]);
         hold on; h_theta=area(theta_freq, theta_pot,'FaceColor',[0.3882    1.0000    0.8392]);
         hold on; h_alpha=area(alpha_freq, alpha_pot,'FaceColor',[0.8078    0.6392    1.0000]);
         hold on; h_sigma=area(sigma_freq,sigma_pot, 'FaceColor',[0.7490    0.8510    1.0000]);
         hold on; h_beta=area(beta_freq, beta_pot,'FaceColor',[1     1     0]);
         
         set(gca, 'XTick', unique([delta_start delta_end theta_end alpha_end beta_end, get(gca, 'XTick')]));
         xlabel('Frequency [Hz]'); ylabel('EEG Power'); %%set assi
         title('PSD of EEG Signal'); %%generico per lo script cond_subj
         lh = legend([h_spectrum h_delta h_theta h_alpha h_beta,h_sigma], ...
        'Spettro', 'Potenza delta', 'Potenza theta','Potenza alpha', 'Potenza sigma', 'Potenza beta');

        set(gca,'FontSize',10)
     end
end