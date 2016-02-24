function myB = myB_20111219_SIRPv0_test_v2a( VECT) ;
% 
% supposedly works with spm8
%
% need to load VECT manually by doing :
%
%	VECT = cafe_central( {'t=dur_v0,freq_v0'}) ;
%
%
% $Id: spm_batch_dur_v0_1st_level_SPM_091001a.mat v0.01 2009-10-01 21:22:50 pisti $
% $Id: myB_091014_mismatch_dur_freq_v0_PAUSE_paramT_v0a v0.02 2009-10-15 00:09:45 pisti $
% $Id: myB_091201_prosody_v5_4conditions_v0a.m v0.03 2009-12-01 10:53:45 pisti $
% $Id: myB_091201_prosody_v5_4conditions_v0a.m v0.04 2009-12-01 18:03:43 pisti $
% $Id: myB_091202_prosody_v5_4conditionsAll_v0b.m v0.04 2009-12-02 22:33:02 pisti $
% $Id: myB_20111219_SIRPv0_test_v0a.m v0.01 2011-12-19 21:40:25 pisti $
% $Id: myB_20111219_SIRPv0_test_v1a.m v0.02 2011-12-21 14:52:04 pisti $ added TARGET and FOIL
% $Id: myB_20111219_SIRPv0_test_v2a.m v0.03 2011-12-21 23:05:22 pisti $ changed order of contrasts , corrected for smaller then expected number of series 
cd /home/pisti/munka/USA/Boston/BWH/ms/analysis/

    %%%% propaganda
    Mout					= IAM( mfilename) ;
    ILM.date					= Mout.date ;
    ILM.version					= Mout.version ;
    ILM.mfilename				= Mout.name ;
    clear Mout

    %%%% home
    [ status output]				= unix([ 'echo $HOME']) ;
    myB.s.home					= '/home/pisti' ;

    %%%% our main data path
    myB.s.path					= '/data/kodaly/r/ms/' ;
    
    %%%% task
% $$$     myB.s.task					= 'dur_v0' ;
% $$$     myB.s.task					= { 'freq_v0' ; 'dur_v0' } ;
% $$$     myB.s.task					= { 'prosody_v5' } ;
    myB.s.task					= { 'SIRPv0' } ;
    
    %%%% our OUTPUT path
% $$$     myB.s.pathSPM				= [ 'SPM_' myB.s.task '_block_paramT_091001' ] ;
    myB.s.pathSPM				= ILM.mfilename ;
    
    %%%% mask
% $$$     myB.s.Mask					= '/home/pisti/programs/matlab/SPM_related/pisti/Mask_ROIs_EPI_brain_20_255_79x95x69_40o57o26.img' ;
% $$$     myB.s.Mask					= '/home/pisti/programs/matlab/SPM_related/pisti/Mask_ROIs_EPI_brain_20_255_79x95x69_40o57o26.nii' ;
    myB.s.Mask					= '/home/pisti/programs/matlab/SPM_related/pisti/Mask_ROIs_EPI_brain_20_255_79x95x68.nii' ;

    %%%% path of this mfile
    [a,b,c]					= fileparts( mfilename('fullpath')) ;
    myB.s.pathMfile				= a ;
    
    %%%% template job file
% $$$     myB.s.jobfile				= '/data/kodaly/r/dickey/jobs/job_mat_fMRI_BL_040686_SPM_dur_v0_AB_block_paramT_091001.mat' ;
% $$$     myB.s.jobfile				= '/data/kodaly/r/dickey/jobs/job_mat_fMRI_BL_040686_SPM_mismatch_dur_frequ_v0_PAUSE_block_paramT_091014.mat' ;
% $$$     myB.s.jobfile				= '/data/kodaly/r/dickey/jobs/job_mat_fMRI_WM_070104_SPM_prosody_v5_4conditions_v0a_091201.mat' ;
% $$$     myB.s.jobfile				= fullfile( myB.s.pathMfile, 'job_fMRI_BC_20070518_SPM_20111214_b.mat') ;
    myB.s.jobfile				= fullfile( myB.s.pathMfile, 'job_fMRI_BC_20070518_SPM_20111221_9cond.mat') ;
    
    
    %%%%
    myB.s.fnameFilter				= '^we.*' ;	%%%% '^sw.*'

    % ls /data/kodaly/r/dickey/fMRI_??_??????/logfiles/STIM_fMRI_??_??????_sentence_prosody_v5_*_B.mat
    
    %%%% determines N cycles of loop
    mySize					= size( VECT.ExtrSubj, 2) ;

    %%%% The Subject Loop ------------------------------------------------------------------------------------------------------
    for counter_subj = 1 : mySize
	
	mySubj					= VECT.ExtrSubj{ counter_subj} ;
	IAM( ILM.mfilename, ILM.version, ILM.date, 'messg', [ 'mySubj is called :\r\t\t\t\t' char(39) mySubj char(39) ]) ;

	%%%% loading job.mat batchfile
	clear matlabbatch
	load( myB.s.jobfile)
	myB.job					= matlabbatch ;
	
 	%%%% removing temporarily 'sess' 
	myB.job{1}.spm.stats.fmri_spec		= rmfield( myB.job{1}.spm.stats.fmri_spec, 'sess') ;
	
	%%%% updating matlabbatch path and maskfile
	myCurrDir				= fullfile( myB.s.path, [ 'fMRI_' mySubj ], 'ssm_preproc', myB.s.pathSPM) ;
	myB.job{1}.spm.stats.fmri_spec.dir	= { myCurrDir } ;
	myB.job{1}.spm.stats.fmri_spec.mask	= { myB.s.Mask } ;
	
	%%%% extracts myV
	myVECT.ExtrSubj				= { mySubj} ;
	myV					= cafe_central( myVECT) ;
	
	%%%% generalizes current subject name
	myB.subject				= eval([ 'myV.' mySubj ]) ;
	
	%%%% goes now through entire loop of possible seriess funct for given exam
	counter_series				= 0 ;
	
	%%%% The Series Loop
	for i = 1 : size( myB.subject.db.funct, 2)
	    
	    %%%% checking if given functs# corresponds with myTask
% $$$ 	    if isempty( strmatch( myB.subject.db.funct( i).task , myB.s.task )) == 0
	    if isempty( strmatch( myB.subject.db.funct( i).task , myB.s.task )) == 0  &  isempty( find( myB.subject.db.funct( i).str == '?')) == 1

		myTask				= myB.s.task{ strmatch( myB.subject.db.funct( i).task , myB.s.task )} ;
		counter_series			= counter_series + 1 ;
		myAleph				= myB.subject.db.funct( i).taskAleph ;
		
		%%%% reloading one series per loop
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series)	= matlabbatch{1}.spm.stats.fmri_spec.sess(1) ;
		
		%%%% loading list of series-specific data files : CAVE : assumes same length of list of files as in loaded jobfile !!!
		[ myFnames, sts]		= spm_select( 'List', fullfile( myB.s.path , [ 'fMRI_' mySubj ] , myB.subject.db.funct( i).series) , myB.s.fnameFilter , '1' ) ;
		
		for j = 1 : size( myFnames, 1)
		    
		    myB.job{1}.spm.stats.fmri_spec.sess( counter_series).scans{j,1}	= ...
			[ fullfile( myB.s.path , [ 'fMRI_' mySubj ] , myB.subject.db.funct( i).series, myFnames( j, :)) ',1'] ;
		    
		end
		
		%%%% regressors : onset vectors
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(1).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load1.encode_time' ] ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(1).onset		= eval([ 'rot90( myB.subject.' myTask '.' myAleph '.load1.encode_time, 3) ;']) ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(2).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load3.encode_time' ] ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(2).onset		= eval([ 'rot90( myB.subject.' myTask '.' myAleph '.load3.encode_time, 3) ;']) ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(3).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load5.encode_time' ] ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(3).onset		= eval([ 'rot90( myB.subject.' myTask '.' myAleph '.load5.encode_time, 3) ;']) ;

		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(4).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load1.target_time' ] ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(4).onset		= eval([ 'rot90( myB.subject.' myTask '.' myAleph '.load1.target_time, 3) ;']) ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(5).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load3.target_time' ] ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(5).onset		= eval([ 'rot90( myB.subject.' myTask '.' myAleph '.load3.target_time, 3) ;']) ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(6).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load5.target_time' ] ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(6).onset		= eval([ 'rot90( myB.subject.' myTask '.' myAleph '.load5.target_time, 3) ;']) ;
		
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(7).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load1.foil_time' ] ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(7).onset		= eval([ 'rot90( myB.subject.' myTask '.' myAleph '.load1.foil_time, 3) ;']) ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(8).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load3.foil_time' ] ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(8).onset		= eval([ 'rot90( myB.subject.' myTask '.' myAleph '.load3.foil_time, 3) ;']) ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(9).name		= [ 'myV.' mySubj '.' myTask '.' myAleph '.load5.foil_time' ] ;
		myB.job{1}.spm.stats.fmri_spec.sess( counter_series).cond(9).onset		= eval([ 'rot90( myB.subject.' myTask '.' myAleph '.load5.foil_time, 3) ;']) ;
	    end
	end
	
	clear matlabbatch
	
	%%%% do job only if series dirs were found
	if counter_series > 0
	    
	    %%%% doing 1st level ...................................................................
	    
	    myJ					= 'job_0' ;
	    
	    %%%% creating directories
	    IAM( ILM.mfilename, ILM.version, ILM.date, 'messg', [ 'creating directory :\n\t\t\t\t' myCurrDir ]) ;
	    [ status, result ]			= unix([ 'mkdir -pv ' myCurrDir ]) ;
	    cd( myCurrDir) ;
	    
	    unix([ 'touch timestamp_' myJ '_begin' ]) ;
	    
	    %%%% saving new job file in pwd
	    matlabbatch				= myB.job ;
	    save job_0.mat matlabbatch
	    
	    %%%% do the job !
	    spm_jobman( 'run', matlabbatch) ;
	    
	    unix([ 'touch timestamp_' myJ '_end' ]) ;
	    clear matlabbatch
	    
	    
	    %%%% doing 'estimate' level ............................................................
	    
        
        
	    myJ					= 'job_1' ;
	    unix([ 'touch timestamp_' myJ '_begin' ]) ;
	    
	    %%%% creating and saving those simple entries
	    myJob1{1}.spm.stats.fmri_est.method.Classical	= 1 ;
	    myJob1{1}.spm.stats.fmri_est.spmmat	= { fullfile( myCurrDir, 'SPM.mat') } ;
	    save job_1.mat myJob1
	    
	    %%%% do the job !
	    %spm_jobman( 'run', myJob1) ;
	    
	    unix([ 'touch timestamp_' myJ '_end' ]) ;
	    
	    
	    %%%% doing 'contrasts' level ...........................................................
	    
	    myJ					= 'job_2' ;
	    unix([ 'touch timestamp_' myJ '_begin' ]) ;
	    
	    %%%% assessing number of regressors for contrast manager
	    mySeriesN				= size( myB.job{1}.spm.stats.fmri_spec.sess,2) ;
	    myCondN				= size( myB.job{1}.spm.stats.fmri_spec.sess(1).cond, 2) ;
	    myParamN				= size( myB.job{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod, 2) ;	%%%% parametric regressor
	    myRegrN				= 1 + myParamN ;							%%%% onset regressor + parametric regressor(s)
	    myColN				= mySeriesN * myCondN * myRegrN + mySeriesN ;
	    
	    load SPM.mat
	    cname				= [] ;
	    c					= [] ;
	    
	    %%%% looping through all tasks (here just one actually = prosody_v5)
	    for myTaskN = 1 : size( myB.s.task,1)

		%%%% ripping task and series name
		[ myFPa myFPb myFPc ]		= fileparts( myB.job{1}.spm.stats.fmri_spec.sess(1).cond(i).name) ;
		[ myFPaa myFPbb myFPcc ]	= fileparts( myFPb) ;
		
		%%%% for SIRPv0 only : ENCODE
		cname{ end+1}			= [ myFPbb ' [encode ALL]' ] ;
% $$$ 		c{ end+1}			= IAM_spm_pad_with_zeros( [ 1 : mySeriesN * myCondN * myRegrN ] , myColN ) ;
% $$$ 		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 1 1 1 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 1 1 1 0 0 0 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
		
		%%%% for SIRPv0 only : TARGET
		cname{ end+1}			= [ myFPbb ' [target ALL]' ] ;
% $$$ 		c{ end+1}			= IAM_spm_pad_with_zeros( [ 1 : mySeriesN * myCondN * myRegrN ] , myColN ) ;
% $$$ 		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 1 1 1 ], 1, mySeriesN) == 1), myColN) ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 1 1 1 0 0 0 ], 1, mySeriesN) == 1), myColN) ;
		
		%%%% for SIRPv0 only : FOIL
		cname{ end+1}			= [ myFPbb ' [foil ALL]' ] ;
% $$$ 		c{ end+1}			= IAM_spm_pad_with_zeros( [ 1 : mySeriesN * myCondN * myRegrN ] , myColN ) ;
% $$$ 		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 1 1 1 ], 1, mySeriesN) == 1), myColN) ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 0 0 0 0 0 0 1 1 1 ], 1, mySeriesN) == 1), myColN) ;
		
		%%%% for ALL
		cname{ end+1}			= [ myFPbb ' [ALL]' ] ;
% $$$ 		c{ end+1}			= IAM_spm_pad_with_zeros( [ 1 : mySeriesN * myCondN * myRegrN ] , myColN ) ;
		c{ end+1}			= IAM_spm_pad_with_zeros( find( repmat( [ 1 1 1 1 1 1 1 1 1 ], 1, mySeriesN) == 1), myColN) ;
        
        cname{ end+1}			= [ myFPbb ' [encode load 1]' ] ;
		c{ end+1}		  = [repmat([1 0 0 0 0 0 0 0 0], [1,4]), 0  0 0 0];
        
 		cname{ end+1}			= [ myFPbb ' [encode load 3]' ] ;
		c{ end+1}		  = [repmat([0 1 0 0 0 0 0 0 0], [1,4]), 0  0 0 0];;
        
        cname{ end+1}			= [ myFPbb ' [encode load 5]' ] ;
		c{ end+1}		  = [repmat([0 0 1 0 0 0 0 0 0], [1,4]), 0  0 0 0];;

        cname{ end+1}			= [ myFPbb ' [target load 1]' ] ;
		c{ end+1}		  = [repmat([0 0 0 1 0 0 0 0 0], [1,4]), 0  0 0 0];;
              
        cname{ end+1}			= [ myFPbb ' [target load 3]' ] ;
		c{ end+1}		  = [repmat([0 0 0 0 1 0 0 0 0], [1,4]), 0  0 0 0];;
        
        cname{ end+1}	  = [ myFPbb ' [target load 5]' ] ;
		c{ end+1}		  = [repmat([0 0 0 0 0 1 0 0 0], [1,4]), 0  0 0 0];;
        
        cname{ end+1}			= [ myFPbb ' [foil load 1]' ] ;
		c{ end+1}		  = [repmat([0 0 0 0 0 0 1 0 0], [1,4]), 0  0 0 0];;
              
        cname{ end+1}			= [ myFPbb ' [foil load 3]' ] ;
		c{ end+1}		  = [repmat([0 0 0 0 0 0 0 1 0], [1,4]), 0  0 0 0];;
        
        cname{ end+1}	  = [ myFPbb ' [foil load 5]' ] ;
		c{ end+1}		  = [repmat([0 0 0 0 0 0 0 0 1], [1,4]), 0  0 0 0];;
        
				
		%%%% looping through all conditions
		for i = 1 : myCondN
		    
		    cname{ end+1}		= myB.job{1}.spm.stats.fmri_spec.sess(1).cond(i).name ;
% $$$ 		    [ myFPa myFPb myFPc ]	= fileparts( myCondStr) ;
% $$$ 		    [ myFPaa myFPbb myFPcc ]	= fileparts( myFPb) ;
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
		    SPM.xCon(end + 1)		= spm_FcUtil( 'Set', cname{i}, 'T', 'c', rot90(c{i}, 3), SPM.xX.xKXs) ;
		end
	    end

% $$$ %			    %%%% checking for task type
% $$$ %			    if isempty( strfind( myB.job{1}.spm.stats.fmri_spec.sess(j).cond(i).name , myB.s.task{ myTaskN})) == 0
% $$$ %				
% $$$ %				myIndex( end+1)		= i ;
% $$$ %	i
% $$$ %	myIndex
% $$$ %	
% $$$ %	% $$$ 			%%%% determines task name
% $$$ %	% $$$ 			myTaskName		= myB.job{1}.spm.stats.fmri_spec.sess(i).cond(1).pmod(1).name ;
% $$$ %				myTaskName		= myB.job{1}.spm.stats.fmri_spec.sess(i).cond(1).name ;
% $$$ %				
% $$$ %			    end
% $$$ %			end
% $$$ %			
% $$$ %			if length( myIndex) > 0
% $$$ %			    
% $$$ %			    myTaskName			= [ myTaskName ' [A+B]' ] 
% $$$ %			    cname{ end+1}		= myTaskName 
% $$$ %			    c{ end+1}			= IAM_spm_pad_with_zeros( [ myRegrN * myIndex], myColN ) 
% $$$ %	
% $$$ %			end
% $$$ %		    end
% $$$ %	
% $$$ %		    
% $$$ %		    %%%% inserting the contrasts
% $$$ %		    for i = 1 : size( cname, 2)
% $$$ %			if i == 1
% $$$ %			    SPM.xCon			= spm_FcUtil( 'Set', cname{i}, 'T', 'c', rot90(c{i}, 3), SPM.xX.xKXs) ;
% $$$ %			else
% $$$ %			    SPM.xCon(end + 1)		= spm_FcUtil( 'Set', cname{i}, 'T', 'c', rot90(c{i}, 3), SPM.xX.xKXs) ;
% $$$ %			end
% $$$ %		    end
% $$$ 
% $$$ % $$$ 		%%%% for onset regressor
% $$$ % $$$ 		cname{ end+1}		= myB.job{1}.spm.stats.fmri_spec.sess(i).cond(1).pmod(1).name
% $$$ 
% $$$ % $$$ 		    
% $$$ % $$$ 		    %%%% in case two series were run -- and not threee or more...!
% $$$ % $$$ 		    if mySeriesN == 2
% $$$ % $$$ 			cname{ end}		= [ cname{ end} ' [A+B]' ] ;
% $$$ % $$$ 		    end
% $$$ % $$$ 		    
% $$$ % $$$ 		    %%%% indexes only work for up to max mySeriesN == 2
% $$$ % $$$ 		    counter_index		= counter_index + 1 ;
% $$$ % $$$ 		    c{ end+1}			= IAM_spm_pad_with_zeros( [ counter_index ( counter_index + ( (mySeriesN-1) * ( myCondN * myRegrN ))) ], myColN ) ;
% $$$ % $$$ 		    
% $$$ % $$$ 		    for k = 1 : myParamN
% $$$ % $$$ 			
% $$$ % $$$ 			
% $$$ % $$$ 			%%%% for parametric regressor(s)
% $$$ % $$$ 			cname{ end+1}		= myB.job{1}.spm.stats.fmri_spec.sess(i).cond(j).pmod(k).name ;
% $$$ % $$$ 			
% $$$ % $$$ 			%%%% in case two series were run -- and not threee or more...!
% $$$ % $$$ 			if mySeriesN == 2
% $$$ % $$$ 			    cname{ end}		= [ cname{ end} ' [A+B]' ] ;
% $$$ % $$$ 			end
% $$$ % $$$ 			
% $$$ % $$$ 			%%%% indexes only work for up to max mySeriesN == 2 
% $$$ % $$$ 			counter_index		= counter_index + 1 ;
% $$$ % $$$ 			c{ end+1}		= IAM_spm_pad_with_zeros( [ counter_index ( counter_index + ( (mySeriesN-1) * ( myCondN * myRegrN ))) ], myColN ) ;
% $$$ % $$$ 			
% $$$ % $$$ 		    end
% $$$ % $$$ 		end
% $$$ % $$$ 	    end
% $$$ % $$$ 		
% $$$ % $$$ 	    %%%% inserting the stuff into SPM variable
% $$$ % $$$ 	    for i = 1 : size( cname, 2)
% $$$ % $$$ 		
% $$$ % $$$ 		if i == 1
% $$$ % $$$ 		    SPM.xCon			= spm_FcUtil( 'Set', cname{i}, 'T', 'c', rot90(c{i}, 3), SPM.xX.xKXs) ;
% $$$ % $$$ 		else
% $$$ % $$$ 		    SPM.xCon(end + 1)		= spm_FcUtil( 'Set', cname{i}, 'T', 'c', rot90(c{i}, 3), SPM.xX.xKXs) ;
% $$$ % $$$ 		end
% $$$ % $$$ 		
% $$$ % $$$ 	end
	
	%%%% run the contrast stuff...
	spm_contrasts(SPM) ;

	    %%%% linking the mean image to pwd
	    unix([ 'ln -sv ../../e*/wmean*nii 1.nii']) ;
	    
	    unix([ 'touch timestamp_' myJ '_end' ]) ;
	
	end
	
	clear mySubj myCurrDir myAleph myV
	
    end
    
    cd( myB.s.home) ;
    spm( 'FigName', 'done') ;
    spm( 'Pointer', 'Arrow') ;
