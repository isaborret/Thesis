%% GC Analysis: Conditional Granger Causality between beat-indexed EEG band
%  power and HF-HRV during sleep, adapted from MVGC toolbox demo (Barnett & Seth, 2014)
% specifically this demo: mvgc_demo_autocov.m

%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] A. B. Barrett, L. Barnett and A. K. Seth, "Multivariate Granger causality
% and generalized variance", _Phys. Rev. E_ 81(4), 2010.
%
% [3] L. Barnett and A. K. Seth, "Behaviour of Granger causality under
% filtering: Theoretical invariance and practical application", _J. Neurosci.
% Methods_ 201(2), 2011.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.

%% Paths
data_path = 'C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\data\derivatives\timeSeries';
output_path = 'C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\data\derivatives\GC';
addpath(genpath('C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\toolboxes\MVGC_toolbox'));

%% Parameters

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     =  10; %20;     % maximum model order for model order estimation
acmaxlags =  500 %1000;   % maximum autocovariance lags (empty for automatic calculation)
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test

electrodes = {'Cz', 'Pz'};
bands      = {'delta', 'theta', 'alpha', 'sigma', 'beta'};

% mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance') => ???????
% fs        = 200;    % sample rate (Hz)
% fres      = [];     % frequency resolution (empty for automatic calculation)
% seed      = 0;      % random seed (0 for unseeded)

%% LOAD DATA
load(fullfile(data_path, 'timeseries_results.mat'), 'results');  % loads struct array 'results'
n_subjects = length(results);

%% INITIALIZE OUTPUT
gc_results = struct();

%% Loop through subjects
for s = 1 %1:n_subjects
 
    subj = results(s);
    fprintf('\n=== Processing %s | %s | dream=%d ===\n', ...
        subj.subject_id, subj.sleep_stage, subj.dream);
 
    gc_results(s).subject_id  = subj.subject_id;
    gc_results(s).dream       = subj.dream;
    gc_results(s).sleep_stage = subj.sleep_stage;
    gc_results(s).n_beats     = subj.n_beats;
 
    % Loop over electrodes
    for e = 1:length(electrodes)
        elec = electrodes{e};
 
        % -----------------------------------------------------------------
        % Assemble the 6 x n_beats data matrix X
        % Row order: [delta; theta; alpha; sigma; beta; HFHRV]
        % Variable indices: 1=delta, 2=theta, 3=alpha, 4=sigma, 5=beta, 6=HRV
        % -----------------------------------------------------------------

        try
            X = [subj.(elec).PEEG;   % 5 x n_beats
                 subj.HFHRV];         % 1 x n_beats — gives 6 x n_beats total
        catch
            warning('%s: Could not assemble data matrix for %s. Skipping.', ...
                subj.subject_id, elec);
            gc_results(s).(elec) = [];
            continue
        end
         
        n_beats = size(X, 2);
        fprintf('  %s: n_beats = %d\n', elec, n_beats);
 
        % Require minimum number of observations for reliable VAR fitting
        if n_beats < 50
            warning('%s %s: fewer than 50 beats, skipping.', subj.subject_id, elec);
            gc_results(s).(elec) = [];
            continue
        end
 
        % Demean each series (required for VAR modelling)
        X = bsxfun(@minus, X, mean(X, 2));
 
        % MVGC expects [nvars x nobs x ntrials]; we have single trial
        X_mvgc = reshape(X, size(X,1), size(X,2), 1);
 
        % -----------------------------------------------------------------
        % Step 1: Model order estimation via AIC
        % -----------------------------------------------------------------
 
        [AIC, ~, moAIC, ~] = tsdata_to_infocrit(X_mvgc, momax, icregmode);
        %This function takes your data matrix X, fits VAR models of increasing order up to momax
        % for each order computes the AIC and BIC values.
            % returns four outputs: 
            % the full AIC curve across all orders, 
            % the full BIC curve, 
            % the order at which AIC is minimised (moAIC),
            % the order at which BIC is minimised (moBIC).
        
        % Plot information criteria
        % figure(1); clf;
        % plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
        % title('Model order estimation');
        % amo = size(AT,3); % actual model order
        % fprintf('\nbest model order (AIC) = %d\n',moAIC);
        % fprintf('best model order (BIC) = %d\n',moBIC);
        % fprintf('actual model order     = %d\n',amo);
        
        if strcmpi(morder, 'AIC')
            mord = moAIC;   % mord
        else
            mord = morder;  % if a fixed integer was supplied
        end
        fprintf('  %s: selected model order (AIC) = %d\n', elec, mord);
 
        % Guard against zero model order
        if mord < 1
            warning('%s %s: AIC selected model order 0, defaulting to 1.', ...
                subj.subject_id, elec);
            mord = 1;
        end
 
        % -----------------------------------------------------------------
        % Step 2: VAR model estimation
        % -----------------------------------------------------------------
 
        [A, SIG] = tsdata_to_var(X_mvgc, mord, regmode);
 
        if isbad(A)
            warning('%s %s: VAR estimation failed (ill-conditioned regression). Skipping.', ...
                subj.subject_id, elec);
            gc_results(s).(elec) = [];
            continue
        end
 
        % -----------------------------------------------------------------
        % Step 3: Autocovariance sequence
        % -----------------------------------------------------------------
 
        [G, info] = var_to_autocov(A, SIG, acmaxlags);
        %var_acinfo(info,true); % report results (and bail out on error)

        % Check for errors (non-stationarity, colinearity, etc.)
        % Prints a warning identifying exactly which subject and electrode failed, 
        % stores an empty result for that case, and uses continue to skip to the next iteration rather than crashing 
        % If there is only a warning (non-fatal), it logs the warning message into the output struct so you can inspect it later, but continues processing
        if info.error
            warning('%s %s: autocovariance computation failed: %s. Skipping.', ...
                subj.subject_id, elec, info.errmsg);
            gc_results(s).(elec) = [];
            continue
        end
 
        if info.warnings
            fprintf('  WARNING for %s %s: %s\n', subj.subject_id, elec, info.warnmsg);
            % Continue but flag in output
            gc_results(s).(elec).acov_warning = info.warnmsg;
        end
 
        % Log spectral radius as diagnostic: tells you how close the VAR is to being unstable
        % values approaching 1 are worth flagging when you review results
        gc_results(s).(elec).spectral_radius = info.rho;
        fprintf('  %s: spectral radius = %.4f\n', elec, info.rho);
 
        % -----------------------------------------------------------------
        % Step 4: Pairwise-conditional GC (time domain)
        %
        % autocov_to_pwcgc returns a 6x6 matrix F where F(i,j) is GC
        % from variable j to variable i, conditioned on all others.
        % Diagonal entries are NaN.
        % Variable 6 = HRV; variables 1-5 = EEG bands
        % -----------------------------------------------------------------
 
        F = autocov_to_pwcgc(G); % =core computation: takes the autocovariance sequence G and computes the full 6×6 matrix of pairwise-conditional GC values
        % Each entry F(i,j) answers the question: how much does knowing the past of variable j improve prediction of variable i, after conditioning on the past of all other four variables?
        % Diagonal is NaN because a variable cannot Granger-cause itself 

        if isbad(F, false)
            warning('%s %s: GC computation failed. Skipping.', ...
                subj.subject_id, elec);
            gc_results(s).(elec) = [];
            continue
        end

        % Significance test using theoretical null distribution, adjusting for multiple
        % hypotheses.
        % pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
        % sig  = significance(pval,alpha,mhtc);
        % % Plot time-domain causal graph, p-values and significance.
        % figure(2); clf;
        % sgtitlex('Pairwise-conditional Granger causality - time domain');
        % subplot(1,3,1);
        % plot_pw(F);
        % title('Pairwise-conditional GC');
        % subplot(1,3,2);
        % plot_pw(pval);
        % title('p-values');
        % subplot(1,3,3);
        % plot_pw(sig);
        % title(['Significant at p = ' num2str(alpha)])
 
        % -----------------------------------------------------------------
        % Step 5: Extract GC scalars of interest
        %   five values from row 6 (HRV as target, each EEG band as source, giving EEG-to-HRV GC) 
        %   and five values from column 6 (HRV as source, each EEG band as target, giving HRV-to-EEG GC)
        %   EEG -> HRV: F(6, band_idx) — HRV (row 6) predicted by EEG band
        %   HRV -> EEG: F(band_idx, 6) — EEG band (row) predicted by HRV (col 6)
        % -----------------------------------------------------------------
 
        for b = 1:length(bands)
            band = bands{b};
            band_idx = b;   % delta=1, theta=2, alpha=3, sigma=4, beta=5
 
            % GC from EEG band to HRV
            gc_results(s).(elec).(['EEG_to_HRV_' band]) = F(6, band_idx);
 
            % GC from HRV to EEG band
            gc_results(s).(elec).(['HRV_to_EEG_' band]) = F(band_idx, 6);
        end
 
        % Store full GC matrix for potential later inspection
        gc_results(s).(elec).F_full     = F;
        gc_results(s).(elec).morder     = mord;
        gc_results(s).(elec).n_beats    = n_beats;
 
        % -----------------------------------------------------------------
        % Step 6: Subject-level F-test (diagnostic only — not primary inference)
        % Primary inference is the group-level LMM in R
        % nobs and ntrials arguments: single trial, so ntrials=1
        % nvars-2 = conditioning set size = 4
        % -----------------------------------------------------------------
 
        nvars = 6;
        pval = mvgc_pval(F, mord, n_beats, 1, 1, 1, nvars-2, 'F');
        gc_results(s).(elec).pval_full = pval;
 
    end  % electrode loop
 
    fprintf('  Done.\n');
 
end  % subject loop
 
%% -------------------------------------------------------------------------
%  SAVE RESULTS
% --------------------------------------------------------------------------
 
save(fullfile(output_path, 'gc_results.mat'), 'gc_results');
fprintf('\nGC results saved to gc_results.mat\n');
 
%% -------------------------------------------------------------------------
%  EXPORT TO TABLE FOR R (optional convenience)
%  Creates a flat CSV with one row per subject per electrode
% --------------------------------------------------------------------------
 
fid = fopen(fullfile(output_path, 'gc_results_for_R.csv'), 'w');
header = 'subject_id,dream,sleep_stage,n_beats,electrode,morder,spectral_radius';
for b = 1:length(bands)
    header = [header ',' sprintf('EEG_to_HRV_%s', bands{b})];
end
for b = 1:length(bands)
    header = [header ',' sprintf('HRV_to_EEG_%s', bands{b})];
end
fprintf(fid, '%s\n', header);
 
for s = 1:length(gc_results)
    for e = 1:length(electrodes)
        elec = electrodes{e};
        if isempty(gc_results(s).(elec))
            continue
        end
        r = gc_results(s).(elec);
        row = sprintf('%s,%d,%s,%d,%s,%d,%.6f', ...
            gc_results(s).subject_id, ...
            gc_results(s).dream, ...
            gc_results(s).sleep_stage, ...
            gc_results(s).n_beats, ...
            elec, ...
            r.morder, ...
            r.spectral_radius);
        for b = 1:length(bands)
            row = [row ',' sprintf('%.8f', r.(['EEG_to_HRV_' bands{b}]))];
        end
        for b = 1:length(bands)
            row = [row ',' sprintf('%.8f', r.(['HRV_to_EEG_' bands{b}]))];
        end
        fprintf(fid, '%s\n', row);
    end
end
 
fclose(fid);
fprintf('Flat CSV exported to gc_results_for_R.csv\n');
 
%% Old part of demo
% %% Model order estimation (<mvgc_schema.html#3 |A2|>)
% 
% 
% % Select model order.
% if     strcmpi(morder,'actual')
%     morder = amo;
%     fprintf('\nusing actual model order = %d\n',morder);
% elseif strcmpi(morder,'AIC')
%     morder = moAIC;
%     fprintf('\nusing AIC best model order = %d\n',morder);
% elseif strcmpi(morder,'BIC')
%     morder = moBIC;
%     fprintf('\nusing BIC best model order = %d\n',morder);
% else
%     fprintf('\nusing specified model order = %d\n',morder);
% end
% 
% %% VAR model estimation (<mvgc_schema.html#3 |A2|>)
% % Estimate VAR model of selected order from data.
% ptic('\n*** tsdata_to_var... ');
% [A,SIG] = tsdata_to_var(X,morder,regmode);
% ptoc;
% % Check for failed regression
% assert(~isbad(A),'VAR estimation failed');
% % NOTE: at this point we have a model and are finished with the data! - all
% % subsequent calculations work from the estimated VAR parameters A and SIG.
% %% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)
% % The autocovariance sequence drives many Granger causality calculations (see
% % next section). Now we calculate the autocovariance sequence G according to the
% % VAR model, to as many lags as it takes to decay to below the numerical
% % tolerance level, or to acmaxlags lags if specified (i.e. non-empty).
% ptic('*** var_to_autocov... ');
% [G,info] = var_to_autocov(A,SIG,acmaxlags);
% ptoc;
% % The above routine does a LOT of error checking and issues useful diagnostics.
% % If there are problems with your data (e.g. non-stationarity, colinearity,
% % etc.) there's a good chance it'll show up at this point - and the diagnostics
% % may supply useful information as to what went wrong. It is thus essential to
% % report and check for errors here.
% var_acinfo(info,true); % report results (and bail out on error)
% %% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)
% % Calculate time-domain pairwise-conditional causalities - this just requires
% % the autocovariance sequence.
% ptic('*** autocov_to_pwcgc... ');
% F = autocov_to_pwcgc(G);
% ptoc;
% % Check for failed GC calculation
% assert(~isbad(F,false),'GC calculation failed');
% % Significance test using theoretical null distribution, adjusting for multiple
% % hypotheses.
% pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
% sig  = significance(pval,alpha,mhtc);
% % Plot time-domain causal graph, p-values and significance.
% figure(2); clf;
% sgtitlex('Pairwise-conditional Granger causality - time domain');
% subplot(1,3,1);
% plot_pw(F);
% title('Pairwise-conditional GC');
% subplot(1,3,2);
% plot_pw(pval);
% title('p-values');
% subplot(1,3,3);
% plot_pw(sig);
% title(['Significant at p = ' num2str(alpha)])
% % For good measure we calculate Seth's causal density (cd) measure - the mean
% % pairwise-conditional causality. We don't have a theoretical sampling
% % distribution for this.
% cd = mean(F(~isnan(F)));
% fprintf('\ncausal density = %f\n',cd);
