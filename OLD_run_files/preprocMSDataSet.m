function pp = preprocMSDataSet( VECT, step) ;
% 

% %     %%%% propaganda
% %     Mout					= IAM( mfilename) ;
% %     ILM.date					= Mout.date ;
% %     ILM.version					= Mout.version ;
% %     ILM.mfilename				= Mout.name ;
% %     clear Mout

    %%%% home
    [ status output]				= unix([ 'echo $HOME']) ;
    pp.s.home					= '/home/fjanoos' ;

    %%%% our main data path
    pp.s.path					= '/data/kodaly/r/ms/' ;
    
    %%%% task
    pp.s.task					= { 'SIRPv0' } ;
    
    %%%% our OUTPUT path
% $$$     pp.s.pathSPM				= [ 'SPM_' pp.s.task '_block_paramT_091001' ] ;
% %     pp.s.pathSPM				= ILM.mfilename ;
    
    %%%% mask
    pp.s.Mask					= '/home/pisti/programs/matlab/SPM_related/pisti/Mask_ROIs_EPI_brain_20_255_79x95x68.nii' ;

    %%%% path of this mfile
    [a,b,c]					= fileparts( mfilename('fullpath')) ;
    pp.s.pathMfile				= a ;
    
    %%%% template job file
    if step == 1
        pp.s.jobfile				= fullfile( pp.s.pathMfile, 'job_fMRI_BC_20070518_SPM_20111221_9cond.mat') ;		% HRF job
    elseif step == 2 
        pp.s.jobfile				= fullfile( pp.s.pathMfile, 'job_fMRI_BC_20070518_SPM_20111221_9cond_FIR.mat') ;	% FIR job
    elseif step == 3 
        pp.s.jobfile				= fullfile( pp.s.pathMfile, 'job_fMRI_BC_20070518_SPM_20111221_9cond_FIR_motion.mat') ;	% ICA job
    end
    
    
    %%%%
    pp.s.fnameFilter				= '^we.*' ;	%%%% '^sw.*'
    pp.s.MPfnameFilter				= '^rp.*txt' ;	%%%% '^sw.*'

    
    %%%% determines N cycles of loop
    mySize					= size( VECT.ExtrSubj, 2) ;

    %%%% The Subject Loop ------------------------------------------------------------------------------------------------------
    for counter_subj = 1 : mySize
	
	mySubj					= VECT.ExtrSubj{ counter_subj} ;
	%%IAM( ILM.mfilename, ILM.version, ILM.date, 'messg', [ 'mySubj is called :\r\t\t\t\t' char(39) mySubj char(39) ]) ;

	%%%% loading job.mat batchfile
	clear matlabbatch
	load( pp.s.jobfile)
	pp.job					= matlabbatch ;
	
 	%%%% removing temporarily 'sess' 
	pp.job{1}.spm.stats.fmri_spec		= rmfield( pp.job{1}.spm.stats.fmri_spec, 'sess') ;
	
	%%%% updating matlabbatch path and maskfile
	switch step 
	    case 1
		preproc_dir = 'preproc_hrf';
	    case 2
		preproc_dir = 'preproc_fir';
	    case 3
		preproc_dir = 'ica';
	end
	myCurrDir				=  fullfile(fullfile(fullfile( pp.s.path, [ 'fMRI_' mySubj ]), 'ssm'), preproc_dir) ;
	pp.job{1}.spm.stats.fmri_spec.dir	= { myCurrDir } ;
	pp.job{1}.spm.stats.fmri_spec.mask	= { pp.s.Mask } ;
	
	%%%% extracts myV
	myVECT.ExtrSubj				= { mySubj} ;
	myV					= cafe_central( myVECT) ;
	
	%%%% generalizes current subject name
	pp.subject				= eval([ 'myV.' mySubj ]) ;
	
	%%%% goes now through entire loop of possible seriess funct for given exam
	counter_series				= 0 ;
	
	%%%% The Series Loop
	for i = 1 : size( pp.subject.db.funct, 2)
	    
	    %%%% checking if given functs# corresponds with myTask
	    if isempty( strmatch( pp.subject.db.funct( i).task , pp.s.task )) == 0  &  isempty( find( pp.subject.db.funct( i).str == '?')) == 1

		myTask				= pp.s.task{ strmatch( pp.subject.db.funct( i).task , pp.s.task )} ;
		counter_series			= counter_series + 1 ;
		myAleph				= pp.subject.db.funct( i).taskAleph ;
		
		%%%% reloading one series per loop
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series)	= matlabbatch{1}.spm.stats.fmri_spec.sess(1) ;
		
		%%%% loading list of series-specific data files : CAVE : assumes same length of list of files as in loaded jobfile !!!
		[ myFnames, sts]		= spm_select( 'List', fullfile( pp.s.path , [ 'fMRI_' mySubj ] , pp.subject.db.funct( i).series) , pp.s.fnameFilter , '1' ) ;
		
		for j = 1 : size( myFnames, 1)
		    
		    pp.job{1}.spm.stats.fmri_spec.sess( counter_series).scans{j,1}	= ...
			[ fullfile( pp.s.path , [ 'fMRI_' mySubj ] , pp.subject.db.funct( i).series, myFnames( j, :)) ',1'] ;
		    
		end
		
		%%%% regressors : onset vectors
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(1).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load1.encode_time' ] ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(1).onset		= eval([ 'rot90( pp.subject.' myTask '.' myAleph '.load1.encode_time, 3) ;']) ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(2).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load3.encode_time' ] ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(2).onset		= eval([ 'rot90( pp.subject.' myTask '.' myAleph '.load3.encode_time, 3) ;']) ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(3).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load5.encode_time' ] ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(3).onset		= eval([ 'rot90( pp.subject.' myTask '.' myAleph '.load5.encode_time, 3) ;']) ;

		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(4).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load1.target_time' ] ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(4).onset		= eval([ 'rot90( pp.subject.' myTask '.' myAleph '.load1.target_time, 3) ;']) ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(5).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load3.target_time' ] ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(5).onset		= eval([ 'rot90( pp.subject.' myTask '.' myAleph '.load3.target_time, 3) ;']) ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(6).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load5.target_time' ] ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(6).onset		= eval([ 'rot90( pp.subject.' myTask '.' myAleph '.load5.target_time, 3) ;']) ;
		
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(7).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load1.foil_time' ] ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(7).onset		= eval([ 'rot90( pp.subject.' myTask '.' myAleph '.load1.foil_time, 3) ;']) ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(8).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load3.foil_time' ] ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(8).onset		= eval([ 'rot90( pp.subject.' myTask '.' myAleph '.load3.foil_time, 3) ;']) ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(9).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load5.foil_time' ] ;
		pp.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(9).onset		= eval([ 'rot90( pp.subject.' myTask '.' myAleph '.load5.foil_time, 3) ;']) ;
		
		
		%%%% PISTI : 2012-04-02 - adjusting motion_parameter fnames - needed by GIFT for component sorting
		if step == 3
		
% $$$ 		    unix(['cp -axv ../preproc_fir/job_0.mat ./.']) ;
% $$$ 		    load job_0.mat
% $$$ 		    pp.job				= matlabbatch ;
		
		    %%%% example : /misc/kodaly/r/data/ms/fMRI_CM_20070511/e4781s4a1/rp_e4781s4a1i_001.txt
		
		    [ myMPfname, sts]		= spm_select( 'List', fullfile( pp.s.path , [ 'fMRI_' mySubj ] , pp.subject.db.funct( i).series) , pp.s.MPfnameFilter , '1' ) ;
		
		    %%%% inserting MP fname into job structure - we assume the order of PM files correspond with imported job/SPM file !  no checks currently !
% $$$ 		    matlabbatch{1}.spm.stats.fmri_spec.sess( counter_series).multi_reg	= myMPfname ;
		    pp.job{1}.spm.stats.fmri_spec.sess( counter_series).multi_reg	= { fullfile( pp.s.path , [ 'fMRI_' mySubj ] , pp.subject.db.funct( i).series, myMPfname) } ;
		
		end
	    end
	end
	
	clear matlabbatch
	
	%%%% do job only if series dirs were found
	if counter_series > 0
	    
	    %%%% doing 1st level ...................................................................
	    
	    myJ					= 'job_0' ;
	    
	    % % 	    %%%% creating directories
	    % % 	    IAM( ILM.mfilename, ILM.version, ILM.date, 'messg', [ 'creating directory :\n\t\t\t\t' myCurrDir ]) ;
	    % % 	    [ status, result ]			= unix([ 'mkdir -pv ' myCurrDir ]) ;
	    % % 	    cd( myCurrDir) ;
	    
	    %%%% FIRDAUS : creating subdir
	    if ~exist( myCurrDir)
		mkdir( fullfile( fullfile( pp.s.path, [ 'fMRI_' mySubj ]), 'ssm'), preproc_dir) ;
	    end
	    cd( myCurrDir) ;
	    
	    unix([ 'touch timestamp_' myJ '_begin' ]) ;
	    
	    %%%% saving new job file in pwd
	    matlabbatch				= pp.job ;
	    save job_0.mat matlabbatch
	    
	    %%%% do the SPECIFY job !
	    spm_jobman( 'run', matlabbatch) ;
	    
	    unix([ 'touch timestamp_' myJ '_end' ]) ;
	    clear matlabbatch

	    %%%% FIRDAUS : STOP
	    if step == 2 | step == 3
		return
	    end

	    
	    %%%% doing 'estimate' level ............................................................
	    
	    myJ					= 'job_1' ;
	    unix([ 'touch timestamp_' myJ '_begin' ]) ;
	    
	    %%%% creating and saving those simple entries
	    myJob1{1}.spm.stats.fmri_est.method.Classical	= 1 ;
	    myJob1{1}.spm.stats.fmri_est.spmmat	= { fullfile( myCurrDir, 'SPM.mat') } ;
	    save job_1.mat myJob1
	    
	    %%%% do the ESTIMATE job !
	    spm_jobman( 'run', myJob1) ;
	    
	    unix([ 'touch timestamp_' myJ '_end' ]) ;
	    
	    
	    %%%% doing 'contrasts' level ...........................................................
	    
	    myJ					= 'job_2' ;
	    unix([ 'touch timestamp_' myJ '_begin' ]) ;
	    
	    %%%% assessing number of regressors for contrast manager
	    mySeriesN				= size( pp.job{1}.spm.stats.fmri_spec.sess,2) ;
	    myCondN				= size( pp.job{1}.spm.stats.fmri_spec.sess(1).cond, 2) ;
	    myParamN				= size( pp.job{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod, 2) ;	%%%% parametric regressor
	    myRegrN				= 1 + myParamN ;							%%%% onset regressor + parametric regressor(s)
	    myColN				= mySeriesN * myCondN * myRegrN + mySeriesN ;
	    
	    load (fullfile( myCurrDir, 'SPM.mat'))
	    cname				= [] ;
	    c					= [] ;
	    
	    %%%% looping through all tasks (here just one actually = prosody_v5)
	    for myTaskN = 1 : size( pp.s.task,1)

		%%%% ripping task and series name
		[ myFPa myFPb myFPc ]		= fileparts( pp.job{1}.spm.stats.fmri_spec.sess(1).cond(i).name) ;
		[ myFPaa myFPbb myFPcc ]	= fileparts( myFPb) ;
		
		%%%% for SIRPv0 only : ENCODE
		cname{ end+1}			= [ myFPbb ' [encode ALL]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 1 1 1 0 0 0 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
		
		%%%% for SIRPv0 only : TARGET
		cname{ end+1}			= [ myFPbb ' [target ALL]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 1 1 1 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
		
		%%%% for SIRPv0 only : FOIL
		cname{ end+1}			= [ myFPbb ' [foil ALL]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 0 0 0 1 1 1 ], 1, mySeriesN) == 1), myColN) ;
		
		%%%% for ALL
		cname{ end+1}			= [ myFPbb ' [ALL]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 1 1 1 1 1 1 1 1 1 ], 1, mySeriesN) == 1), myColN) ;
        
% $$$ 		cname{ end+1}			= [ myFPbb ' [encode load 1]' ] ;
% $$$ 		c{ end+1}			= [repmat([1 0 0 0 0 0 0 0 0], [1,4]), 0  0 0 0];
% $$$         
% $$$  		cname{ end+1}			= [ myFPbb ' [encode load 3]' ] ;
% $$$ 		c{ end+1}			= [repmat([0 1 0 0 0 0 0 0 0], [1,4]), 0  0 0 0];;
% $$$ 		
% $$$ 		cname{ end+1}			= [ myFPbb ' [encode load 5]' ] ;
% $$$ 		c{ end+1}			= [repmat([0 0 1 0 0 0 0 0 0], [1,4]), 0  0 0 0];;
% $$$ 		
% $$$ 		cname{ end+1}			= [ myFPbb ' [target load 1]' ] ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 1 0 0 0 0 0], [1,4]), 0  0 0 0];;
% $$$ 		
% $$$ 		cname{ end+1}			= [ myFPbb ' [target load 3]' ] ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 1 0 0 0 0], [1,4]), 0  0 0 0];;
% $$$ 		
% $$$ 		cname{ end+1}			= [ myFPbb ' [target load 5]' ] ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 0 1 0 0 0], [1,4]), 0  0 0 0];;
% $$$ 		
% $$$ 		cname{ end+1}			= [ myFPbb ' [foil load 1]' ] ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 0 0 1 0 0], [1,4]), 0  0 0 0];;
% $$$ 		
% $$$ 		cname{ end+1}			= [ myFPbb ' [foil load 3]' ] ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 0 0 0 1 0], [1,4]), 0  0 0 0];;
% $$$ 		
% $$$ 		cname{ end+1}			= [ myFPbb ' [foil load 5]' ] ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 0 0 0 0 1], [1,4]), 0  0 0 0];;

		%%%% PISTI : 2012-04-02		
    		cname{ end+1}			= [ myFPbb ' [encode load 1]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 1 0 0 0 0 0 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
% $$$ 		c{ end+1}			= [repmat([1 0 0 0 0 0 0 0 0], [1,mySeriesN]), 0  0 0 0];;
        
 		cname{ end+1}			= [ myFPbb ' [encode load 3]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 1 0 0 0 0 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
% $$$ 		c{ end+1}			= [repmat([0 1 0 0 0 0 0 0 0], [1,mySeriesN]), 0  0 0 0];;
		
		cname{ end+1}			= [ myFPbb ' [encode load 5]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 1 0 0 0 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
% $$$ 		c{ end+1}			= [repmat([0 0 1 0 0 0 0 0 0], [1,mySeriesN]), 0  0 0 0];;
		
		cname{ end+1}			= [ myFPbb ' [target load 1]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 1 0 0 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 1 0 0 0 0 0], [1,mySeriesN]), 0  0 0 0];;
		
		cname{ end+1}			= [ myFPbb ' [target load 3]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 0 1 0 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 1 0 0 0 0], [1,mySeriesN]), 0  0 0 0];;
		
		cname{ end+1}			= [ myFPbb ' [target load 5]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 0 0 1 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 0 1 0 0 0], [1,mySeriesN]), 0  0 0 0];;
		
		cname{ end+1}			= [ myFPbb ' [foil load 1]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 0 0 0 1 0 0 ], 1, mySeriesN) == 1), myColN) ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 0 0 1 0 0], [1,mySeriesN]), 0  0 0 0];;
		
		cname{ end+1}			= [ myFPbb ' [foil load 3]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 0 0 0 0 1 0 ], 1, mySeriesN) == 1), myColN) ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 0 0 0 1 0], [1,mySeriesN]), 0  0 0 0];;
		
		cname{ end+1}			= [ myFPbb ' [foil load 5]' ] ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 0 0 0 0 0 1 ], 1, mySeriesN) == 1), myColN) ;
% $$$ 		c{ end+1}			= [repmat([0 0 0 0 0 0 0 0 1], [1,mySeriesN]), 0  0 0 0];;
        
		%%%% looping through all conditions
		for i = 1 : myCondN
		    
		    cname{ end+1}		= pp.job{1}.spm.stats.fmri_spec.sess(1).cond(i).name ;
		    c{ end+1}			= IAM_spm_pad_with_zeros( i, myColN ) ;
		    
		    %%%% if number of series == 2
		    if mySeriesN == 2
			
			cname{ end}		= [ cname{ end} ' [A+B]' ] ;
			c{ end}			= IAM_spm_pad_with_zeros( [ i i+myCondN ] , myColN ) ;
			
		    end
		end
		
	    end

	    %%%% inserting the contrasts
	    for i = 1 : size( cname, 2)

		if i == 1
		    SPM.xCon			= spm_FcUtil( 'Set', cname{i}, 'T', 'c', rot90(c{i}, 3), SPM.xX.xKXs) ;
		else
		    SPM.xCon(end + 1)		= spm_FcUtil( 'Set', cname{i}, 'T', 'c', rot90(c{i}, 3), SPM.xX.xKXs) ;		%%%% PISTI : here it broke for subject SD
		end
	    end

        %%%% run the contrast stuff...
        spm_contrasts(SPM) ;

	    %%%% linking the mean image to pwd
	    unix([ 'ln -sv ../../e*/wmean*nii 1.nii']) ;
	    
	    unix([ 'touch timestamp_' myJ '_end' ]) ;
	
	end
	
	clear mySubj myCurrDir myAleph myV
	
    end
    
    cd( pp.s.home) ;
    spm( 'FigName', 'done') ;
    spm( 'Pointer', 'Arrow') ;
