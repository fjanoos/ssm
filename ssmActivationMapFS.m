function fs_map = ssmActivationMapFS( sptmap, rf )
% create the stimulus/condition specific activation map (cond_act) in feature space
    
% $Id: ssmActivationMapFS.m v0.01 2012-05-30 18:35:12 fj $

%%%% propaganda
eval([ 'env.logo.' mfilename ' = cafe_logo( mfilename) ;' ])

    sptmap.xSPM.thresDesc = 'FWE';
    sptmap.xSPM.u = rf/50;
    [S,cond_act] = spm_getSPM(sptmap.xSPM);

    try
        topoFDR = spm_get_defaults('stats.topoFDR');
    catch
        topoFDR = true;
    end

    tol = eps*10;
    %-Extract data from cond_act
    %----------------------------------------------------------------------
    S     = cond_act.S;
    VOX   = cond_act.VOX;
    DIM   = cond_act.DIM;
    n     = cond_act.n;
    STAT  = cond_act.STAT;
    df    = cond_act.df;
    u     = cond_act.u;
    M     = cond_act.M;
    k     = cond_act.k;
    
    try, QPs = cond_act.Ps; end
    try, QPp = cond_act.Pp; end
    try, QPc = cond_act.Pc; end
   
    if STAT~='P'
        R     = full(cond_act.R);
        FWHM  = full(cond_act.FWHM);
    end
    try
        units = cond_act.units;
    catch
        units = {'mm' 'mm' 'mm'};
    end
    units{1}  = [units{1} ' '];
    units{2}  = [units{2} ' '];

    DIM       = DIM > 1;              % non-empty dimensions
    D         = sum(DIM);             % highest dimension
    VOX       = VOX(DIM);             % scaling

    if STAT ~= 'P'
        FWHM  = FWHM(DIM);            % Full width at max/2
        FWmm  = FWHM.*VOX;            % FWHM {units}
        V2R   = 1/prod(FWHM);         % voxels to resels
        k     = k*V2R;                % extent threshold in resels
        R     = R(1:(D + 1));         % eliminate null resel counts
        try, QPs = sort(QPs(:)); end  % Needed for voxel   FDR
        try, QPp = sort(QPp(:)); end  % Needed for peak    FDR
        try, QPc = sort(QPc(:)); end  % Needed for cluster FDR
    end

    %-Number and separation for maxima 
    Num    = 3;
    Dis    = 8;


     %------------------------------------------------------------------
    Pz           = spm_P(1,0,u,df,STAT,1,n,S);
    Pu           = spm_P(1,0,u,df,STAT,R,n,S);
    [P Pn Ec Ek] = spm_P(1,k,u,df,STAT,R,n,S);
    
    %-Workaround for negative thresholds
    %----------------------------------------------------------------------
    minz          = abs(min(min(cond_act.Z)));
    zscores       = 1 + minz + cond_act.Z;
    [N Z XYZ A L] = spm_max(zscores,cond_act.XYZ);
    Z             = Z - minz - 1;    
        
        
    %-- 420
    c       = max(A);                                  %-Number of clusters
    try
        NONSTAT = spm_get_defaults('stats.rft.nonstat');
    catch
        NONSTAT = 0;
    end
    if STAT ~= 'P'
        if NONSTAT
            K     = zeros(c,1);
            for i = 1:c
                
                %-Get LKC for voxels in i-th region
                %----------------------------------------------------------
                LKC  = spm_get_data(cond_act.VRpv,L{i});
                
                %-Compute average of valid LKC measures for i-th region
                %----------------------------------------------------------
                valid = ~isnan(LKC);
                if any(valid)
                    LKC = sum(LKC(valid)) / sum(valid);
                else
                    LKC = V2R; % fall back to whole-brain resel density
                end
                
                %-Intrinsic volume (with surface correction)
                %----------------------------------------------------------
                IV   = spm_resels([1 1 1],L{i},'V');
                IV   = IV*[1/2 2/3 2/3 1]';
                K(i) = IV*LKC;
                
            end
            K   = K(A);
        else
            K   = N*V2R;
        end
    end

    %-Convert maxima locations from voxels to mm
    %----------------------------------------------------------------------
    XYZmm = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];
    Pc    = spm_P(c,k,u,df,STAT,R,n,S);  %line 476
    
    
    % --493
    fs_map.dat = {Pc,c};            %-Table data
    TabLin     = 1;                 %-Table data line
    
    %- 497
    %----------------------------------------------------------------------
    HlistXYZ = [];
    while numel(find(isfinite(Z)))
       
        %- 516
        %------------------------------------------------------------------
        [U,i]   = max(Z);           % largest maxima
        j       = find(A == A(i));  % maxima in cluster


        %-Compute peak-level activity for this cluster
        %------------------------------------------------------------------
        if STAT ~= 'P'
            
            % p-values (FWE)
            %--------------------------------------------------------------
            Pz      = spm_P(1,0,   U,df,STAT,1,n,S);  % uncorrected p value
            Pu      = spm_P(1,0,   U,df,STAT,R,n,S);  % FWE-corrected {based on Z}
            [Pk Pn] = spm_P(1,K(i),u,df,STAT,R,n,S);  % [un]corrected {based on K}
            
            % q-values (FDR)
            %--------------------------------------------------------------
            if topoFDR
                Qc  = spm_P_clusterFDR(K(i),df,STAT,R,n,u,QPc); % based on K
                Qp  = spm_P_peakFDR(U,df,STAT,R,n,u,QPp);       % based on Z
                Qu  = [];
            else
                Qu  = spm_P_FDR(U,df,STAT,n,QPs);     % voxel FDR-corrected
                Qc  = [];
                Qp  = [];
            end

            % Equivalent Z-variate
            %--------------------------------------------------------------
            if Pz < tol
                Ze  = Inf;
            else
                Ze  = spm_invNcdf(1 - Pz);
            end
        else
            Pz      = [];
            Pu      = [];
            Qu      = [];
            Pk      = [];
            Pn      = [];
            Qc      = [];
            Qp      = [];
            ws      = warning('off','outOfRangeNormal');
            Ze      = spm_invNcdf(U);
            warning(ws);
        end

        if topoFDR
        [fs_map.dat{TabLin,3:12}] = deal(Pk,Qc,N(i),Pn,Pu,Qp,U,Ze,Pz,XYZmm(:,i));
        else
        [fs_map.dat{TabLin,3:12}] = deal(Pk,Qc,N(i),Pn,Pu,Qu,U,Ze,Pz,XYZmm(:,i));
        end
        TabLin = TabLin + 1;

        
        %- 629
        %------------------------------------------------------------------
        [l q] = sort(-Z(j));                % sort on Z value
        D     = i;
        for i = 1:length(q)
            d = j(q(i));
            if min(sqrt(sum((XYZmm(:,D)-XYZmm(:,d)*ones(1,size(D,2))).^2)))>Dis;

                if length(D) < Num                    

                    % voxel-level activity values 
                    %------------------------------------------------------
                    if STAT ~= 'P'
                        Pz    = spm_P(1,0,Z(d),df,STAT,1,n,S);
                        Pu    = spm_P(1,0,Z(d),df,STAT,R,n,S);
                        if topoFDR
                            Qp = spm_P_peakFDR(Z(d),df,STAT,R,n,u,QPp);
                            Qu = [];
                        else
                            Qu = spm_P_FDR(Z(d),df,STAT,n,QPs);
                            Qp = [];
                        end
                        if Pz < tol
                            Ze = Inf;
                        else
                            Ze = spm_invNcdf(1 - Pz); 
                        end
                    else
                        Pz    = [];
                        Pu    = [];
                        Qu    = [];
                        Qp    = [];
                        ws      = warning('off','outOfRangeNormal');
                        Ze    = spm_invNcdf(Z(d));
                        warning(ws);
                    end


                    % specifically modified line to use hMIPax
                    %------------------------------------------------------
                    tXYZmm = XYZmm(DIM,d);
                    D     = [D d];
                    if topoFDR
                    [fs_map.dat{TabLin,7:12}] = ...
                        deal(Pu,Qp,Z(d),Ze,Pz,XYZmm(:,d));
                    else
                    [fs_map.dat{TabLin,7:12}] = ...
                        deal(Pu,Qu,Z(d),Ze,Pz,XYZmm(:,d));
                    end
                    TabLin = TabLin+1;
                end
            end
        end %for i
        Z(j) = NaN;     % Set local maxima to NaN
    end                 % while numel
