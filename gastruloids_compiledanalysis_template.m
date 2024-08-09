% representative code to compile and average axial profile measurments across conditions
% written by H McNamara
% N.B.: template modified and adapted for specific experiments
% v1 April 2022
% updated Aug 2024

cd 'normdata' % navigate to appropriate directory


% get list of data files
files = ls('*.mat');
[nf,~] = size(files);

% setup arrays to store profiles
GFPprofArray = cell(nf,1);
DsRedprofArray = cell(nf,1);
iRFPprofArray = cell(nf,1);
FracprofArray = cell(nf,1);
LArray = cell(nf,1);
% and to store conditions
% here, averaging over differet dox induction times, and different
% fixation times
% can adjust labels depending on experiment
tfix = zeros(nf,1) ; % store fixation time of each profile
tdox = zeros(nf,1); % store dox induciton time of each profile
% for each file
for f = 1:nf
    f
    % read filename
% EXAMPLE: fn = 't120_dox072_stack1_maxmin_oid1_data'; 
    fn = files(f,:);
    % read conditions
    % assumes convention that fixation time is after 't' and dox time is
    % after 'dox'
    tidx = strfind(fn,'t'); tidx = tidx(1);
    didx = strfind(fn,'dox');
    tfix(f) = str2num(fn(tidx+4:tidx+6));
    tdox(f) = str2num(fn(didx+3:didx+5));

    load(fn);
    GFPprofArray{f} = GFPprof;
    DsRedprofArray{f} = DsRedprof;
    iRFPprofArray{f} = iRFPprof;
    FracprofArray{f} = FracLabelProf;
    
    LArray{f} = L;
end
save('compiled_data.mat', 'GFPprofArray', 'DsRedprofArray', 'iRFPprofArray', 'FracprofArray', 'LArray', 'tfix', 'tdox', 'nf');

% average profiles across conditions
%  ID each unique condition
tfix_unique = sort(unique(tfix));
tdox_unique = sort(unique(tdox));

% set up arrays to store length-normalized traces
GFPprofArray_Lnorm = {};
DsRedprofArray_Lnorm = {};
iRFPprofArray_Lnorm = {};
FracprofArray_Lnorm = {};
 
% set up arrays to store conditions
GFPprofArray_Avg = {};
DsRedprofArray_Avg = {};
FracprofArray_Avg = {};
iRFPprofArray_Avg = {};
testArray_Avg = {}; % to test that averaging done correctly

GFPprofArray_Std = {};
DsRedprofArray_Std = {};
FracprofArray_Std = {};
iRFPprofArray_Std = {};

% also set up arrays for version where length is renormalized to unit A-P
% axis

GFPprofArray_Avg_Lnorm = {};
DsRedprofArray_Avg_Lnorm = {};
FracprofArray_Avg_Lnorm = {};
iRFPprofArray_Avg_Lnorm = {};

GFPprofArray_Std_Lnorm = {};
DsRedprofArray_Std_Lnorm = {};
FracprofArray_Std_Lnorm = {};
iRFPprofArray_Std_Lnorm = {};
LNorm = [0:0.01:1]; %interpolate to percentiles of length

oidnums = [];
tfix_comp = [];
tdox_comp = [];


% average profiles
cidx = 0;
traceidx = 0;
for a = 1:length(tfix_unique)
    for b = 1:length(tdox_unique)
        % find all gastruloids corresponding to this condition
        idxA = tfix == tfix_unique(a);
        idxB = tdox == tdox_unique(b);
        idxC = and(idxA,idxB);

        if(sum(idxC)>0) % check that we have data for this condition
            cidx = cidx + 1; % increment counter

            oidnums(cidx) = sum(idxC); % track oids per condition
            tfix_comp(cidx) = tfix_unique(a);
            tdox_comp(cidx) = tdox_unique(b);
            
            % find file indices which have this conition
            fileidx = find(idxC==1);
            % setup averaging arrays
            gfp_sum = [];
            dsred_sum = [];
            irfp_sum = [];
            frac_sum = [];
            counts = [];
            l_comp = [];

            % and for normalized data
            gfp_sum_lnorm = zeros(size(LNorm));
            dsred_sum_lnorm = zeros(size(LNorm));
            irfp_sum_lnorm = zeros(size(LNorm));
            frac_sum_lnorm = zeros(size(LNorm));

            for c = 1:sum(idxC)
                oididx = fileidx(c); % get index for gastruloid
                % pull data
                gfp = GFPprofArray{oididx};
                dsred = DsRedprofArray{oididx};
                frac = gfp ./ (gfp + dsred);
                irfp = iRFPprofArray{oididx};
                l = LArray{oididx};

                gfp_interp = interp1(l/max(l),gfp,LNorm);
                dsred_interp = interp1(l/max(l),dsred,LNorm);
                irfp_interp = interp1(l/max(l),irfp,LNorm);
                frac_interp = gfp_interp./(gfp_interp + dsred_interp);

                gfp_sum_lnorm = gfp_sum_lnorm + gfp_interp;
                dsred_sum_lnorm = dsred_sum_lnorm + dsred_interp;
                irfp_sum_lnorm = irfp_sum_lnorm + irfp_interp;
                frac_sum_lnorm = frac_sum_lnorm + frac_interp;


                sz_current = length(gfp);
                sz_sum = length(gfp_sum);

                if sz_current>sz_sum % if longer than longest so far
                    % then expand arrays
                    tmp = zeros(sz_current,1);
                    tmp((sz_current-sz_sum+1):end) = gfp_sum;
                    gfp_sum = tmp;
                    tmp = zeros(sz_current,1);
                    tmp((sz_current-sz_sum+1):end) = dsred_sum;
                    dsred_sum= tmp;
                    tmp = zeros(sz_current,1);
                    tmp((sz_current-sz_sum+1):end) = frac_sum;
                    frac_sum= tmp;
                    tmp = zeros(sz_current,1);
                    tmp((sz_current-sz_sum+1):end) = irfp_sum;
                    irfp_sum= tmp;
                    tmp = zeros(sz_current,1);
                    tmp((sz_current-sz_sum+1):end) = counts;
                    counts= tmp;
                    % and update axis label
                    l_comp = l(end:-1:1); % reverse orientation b/c measured from posterior pole now

                    % update sz_sum variable in case lengthened
                    sz_sum = max(sz_sum, sz_current);
                
                end
                
                % add to accumulators
                % indexing to align to posterior most postiion for averaging
                gfp_sum((sz_sum-sz_current+1):end) = gfp_sum((sz_sum-sz_current+1):end) + gfp';
                dsred_sum((sz_sum-sz_current+1):end) = dsred_sum((sz_sum-sz_current+1):end) + dsred';
                frac_sum((sz_sum-sz_current+1):end) = frac_sum((sz_sum-sz_current+1):end) + frac';
                irfp_sum((sz_sum-sz_current+1):end) = irfp_sum((sz_sum-sz_current+1):end) + irfp';
                counts((sz_sum-sz_current+1):end) = counts((sz_sum-sz_current+1):end) + 1;

                % save individual traces
                traceidx = traceidx + 1;

                GFPprofArray_Lnorm{traceidx} = gfp_interp;
                DsRedprofArray_Lnorm{traceidx}  = dsred_interp;
                iRFPprofArray_Lnorm{traceidx}  = irfp_interp;
                FracprofArray_Lnorm{traceidx}  = frac_interp;
 
                
            end
            

            
            % calculate averages
            gfp_avg = gfp_sum ./ counts;
            dsred_avg = dsred_sum ./ counts;
            frac_avg = frac_sum ./ counts;
            irfp_avg = irfp_sum ./ counts;

            gfp_avg_lnorm = gfp_sum_lnorm/sum(idxC);
            dsred_avg_lnorm = dsred_sum_lnorm/sum(idxC);
            frac_avg_lnorm= frac_sum_lnorm/sum(idxC);
            irfp_avg_lnorm = irfp_sum_lnorm/sum(idxC);


            % now, loop again to calculate variances
            gfp_var = zeros(size(gfp_sum));
            dsred_var = zeros(size(gfp_sum));
            frac_var = zeros(size(gfp_sum));
            irfp_var = zeros(size(gfp_sum));

            gfp_var_lnorm = zeros(size(LNorm));
            dsred_var_lnorm = zeros(size(LNorm));
            irfp_var_lnorm = zeros(size(LNorm));
            frac_var_lnorm = zeros(size(LNorm));

            for c = 1:sum(idxC)
                oididx = fileidx(c); % get index for gastruloid
                % pull data again
                gfp = GFPprofArray{oididx};
                dsred = DsRedprofArray{oididx};
                frac = gfp ./ (gfp + dsred);
                irfp = iRFPprofArray{oididx};
                l = LArray{oididx};

                %interpolate
                gfp_interp = interp1(l/max(l),gfp,LNorm);
                dsred_interp = interp1(l/max(l),dsred,LNorm);
                irfp_interp = interp1(l/max(l),irfp,LNorm);
                frac_interp = gfp_interp./(gfp_interp + dsred_interp);

                sz_current = length(gfp);
                sz_sum = length(gfp_sum);

                dgfp = (gfp' - gfp_avg((sz_sum-sz_current+1):end)).^2; 
                ddsred = (dsred' - dsred_avg((sz_sum-sz_current+1):end)).^2; 
                dfrac = (frac' - frac_avg((sz_sum-sz_current+1):end)).^2;
                dirfp = (irfp' - irfp_avg((sz_sum-sz_current+1):end)).^2; 

                % and for interpolated oids
                dgfp_interp = (gfp_interp - gfp_avg_lnorm).^2;
                ddsred_interp = (dsred_interp - dsred_avg_lnorm).^2;
                dirfp_interp = (irfp_interp - irfp_avg_lnorm).^2;
                dfrac_interp = (frac_interp - frac_avg_lnorm).^2;

                gfp_var((sz_sum-sz_current+1):end) = gfp_var((sz_sum-sz_current+1):end) + dgfp;
                dsred_var((sz_sum-sz_current+1):end) = dsred_var((sz_sum-sz_current+1):end) + ddsred;
                frac_var((sz_sum-sz_current+1):end) = frac_var((sz_sum-sz_current+1):end) + dfrac;
                irfp_var((sz_sum-sz_current+1):end) = irfp_var((sz_sum-sz_current+1):end) + dirfp;

                gfp_var_lnorm = gfp_var_lnorm + dgfp_interp;
                dsred_var_lnorm = dsred_var_lnorm + ddsred_interp;
                irfp_var_lnorm = irfp_var_lnorm + dirfp_interp;
                frac_var_lnorm = frac_var_lnorm + dfrac_interp;
            end

          
            % trim to only save axial positions with some minimum number of
            % gatruloids
            mincounts = 2;
            keepidx = counts >= mincounts;
            % trim accordingly
            gfp_avg = gfp_avg(keepidx);
            gfp_var = gfp_var(keepidx);
            dsred_avg = dsred_avg(keepidx);
            dsred_var = dsred_var(keepidx);
            frac_avg = frac_avg(keepidx);
            frac_var = frac_var(keepidx);
            irfp_avg= irfp_avg(keepidx);
            irfp_var = irfp_var(keepidx);
            counts = counts(keepidx);
            l_comp = l_comp(keepidx);

            % convert to std deviations
            gfp_std = sqrt(gfp_var./(counts - 1));
            dsred_std = sqrt(dsred_var./(counts - 1));
            frac_std = sqrt(frac_var./(counts - 1));
            irfp_std = sqrt(irfp_var./(counts-1));

            gfp_std_lnorm = sqrt(gfp_var_lnorm/(sum(idxC)-1));
            dsred_std_lnorm = sqrt(dsred_var_lnorm/(sum(idxC)-1));
            irfp_std_lnorm = sqrt(irfp_var_lnorm/(sum(idxC)-1));
            frac_std_lnorm = sqrt(frac_var_lnorm/(sum(idxC)-1));

            % and save
            GFPprofArray_Avg{cidx} = gfp_avg;
            DsRedprofArray_Avg{cidx} = dsred_avg;
            FracprofArray_Avg{cidx} = frac_avg;
            iRFPprofArray_Avg{cidx} = irfp_avg;
            LArray_comp{cidx} = l_comp;

            GFPprofArray_Std{cidx} = gfp_std;
            DsRedprofArray_Std{cidx} = dsred_std;
            FracprofArray_Std{cidx} = frac_std;
            iRFPprofArray_Std{cidx} = irfp_std;

            GFPprofArray_Avg_Lnorm{cidx} = gfp_avg_lnorm;
            DsRedprofArray_Avg_Lnorm{cidx} = dsred_avg_lnorm;
            FracprofArray_Avg_Lnorm{cidx} = frac_avg_lnorm;
            iRFPprofArray_Avg_Lnorm{cidx} = irfp_avg_lnorm;

            GFPprofArray_Std_Lnorm{cidx} = gfp_std_lnorm;
            DsRedprofArray_Std_Lnorm{cidx} = dsred_std_lnorm;
            FracprofArray_Std_Lnorm{cidx} = frac_std_lnorm;
            iRFPprofArray_Std_Lnorm{cidx} = irfp_std_lnorm;
        end
    end
end

% save output

save('averaged_data.mat', 'GFPprofArray_Avg', 'DsRedprofArray_Avg', 'FracprofArray_Avg', 'iRFPprofArray_Avg', 'GFPprofArray_Std', 'DsRedprofArray_Std', 'FracprofArray_Std', 'iRFPprofArray_Std', 'LArray_comp', 'tfix_comp', 'tdox_comp', 'oidnums', 'GFPprofArray_Std_Lnorm', 'GFPprofArray_Std_Lnorm', 'DsRedprofArray_Avg_Lnorm', 'DsRedprofArray_Std_Lnorm', 'FracprofArray_Avg_Lnorm', 'FracprofArray_Std_Lnorm', 'iRFPprofArray_Avg_Lnorm', 'iRFPprofArray_Std_Lnorm');
            

% plot axial profiles?
% using auxiliary function boundedline.m to plate std devition as shaded
% region about mean

cmap = colormap();

% all profiles for differnet times
close all; figure(); hold on;

boundedline(LNorm, double(FracprofArray_Avg_Lnorm{1}), double(FracprofArray_Std_Lnorm{1}), 'cmap', cmap(1,:));
boundedline(LNorm, double(FracprofArray_Avg_Lnorm{2}), double(FracprofArray_Std_Lnorm{2}), 'cmap', cmap(2,:));
boundedline(LNorm, double(FracprofArray_Avg_Lnorm{3}), double(FracprofArray_Std_Lnorm{3}), 'cmap', cmap(3,:));
boundedline(LNorm, double(FracprofArray_Avg_Lnorm{4}), double(FracprofArray_Std_Lnorm{4}), 'cmap', cmap(4,:));

% etc...
