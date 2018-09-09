% fMRI META-ANALYSIS WITH ADJUSTMENT FOR VARIABLE DROPOUT
% Jo Cutler - J.Cutler@sussex.ac.uk
%
% This script runs an adaptation of the meta-analysis technique AES:SDM,
% created by Joaquim Radua, Anton Albajes-Eizagirre and others - www.sdmproject.com/people/
% 
% INSTRUCTIONS FOR USE:
% The following script can be used in two ways:
%       - with provided data as an example
%       - on your own data
% Further instructions are given for each use separately below.
% 
%% DEPENDENCIES: You need SPM downloaded and on the matlab path to run this script. 
% To download SPM visit www.fil.ion.ucl.ac.uk/spm/
% Follow the instructions online then add the location of SPM on your
% computer here if it's not already on your path:
try
   spm('defaults','fmri'); % tries to load spm fmri defaults to check whether spm on path
catch % if can't find spm functions, add spm to path
    try
        addpath % ADD YOUR PATH TO SPM HERE e.g. addpath C:\Programs\spm12
    catch
        spmdirectory = uigetdir([],'Select the spm folder'); % if spm path not added above, opens
        addpath(spmdirectory); % adds spm directory to path
    end
    spm('defaults','fmri'); % loads spm fmri defaults
end

%% RUNNING AN EXAMPLE ANALYSIS WITH THE DATA PROVIDED
% Download and unzip the data folder into your preferred location
% Run this script and select the directory containing all of the unzipped
% files in the window that appears.

% There are 18 datafiles for this example which were randomly allocated
% between two groups of 9 for the mixed effects analysis.

% We are incredibly grateful to the authors of these datasets for uploading
% their data to NeuroVault - https://neurovault.org/ for reuse such as in
% this project. More information on the data can be found in document 'NeuroVault_studies'.

%% USING THE METHOD ON YOUR OWN META-ANALSYS
% This script requires the first steps of a meta-analysis in AES:SDM to be run BEFORE
% running the adjusted analysis here.
% Please visit www.sdmproject.com, download the software and follow the
% steps described up to at least the preprocessing stage.
% The instructions for this will explain what types of data can be included
% (peaks and maps) and demonstrate how to organise the data.
% You may want to run all the stages of the AES:SDM method to get the
% unadjusted results for your meta-analysis. These can be compared to the
% results from this script to check that it has worked correctly (look for
% the files with "All" in the name from this output to compare with
% AES:SDM).

% IMPORTANT - AES:SDM preprocessing outputs the files as .nii.gz but this
% script requires .nii input. Simply select all of the files beginning
% "pp_" from the AES:SDM output and unzip them into the same folder. (These
% will be pp_(study)_var.nii.gz and pp_(study).nii.gz files.) Also unzip
% the file sdm_mask.nii.gz so it becomes sdm_mask.nii.

% When you have done this, run this script and select the directory
% containing all of the unzipped files as well as the sdm_table in the
% window that appears.

%% CHANGE LOG

% 07/09/18 (JC): Initial script created in MATLAB 2017, using SPM12 on a
% Windows 64-bit PC running Windows 7.

%%

clear all % clears the workspace
directory = uigetdir([],'Select the meta-analysis folder');
cd(directory); % changes to the above directory

answer = questdlg('What type of analysis would you like to run?','Analysis type','Single group mean (random effects)','Two groups mean and comparison (mixed effects)','Single group mean (random effects)');
switch answer
    case 'Single group mean (random effects)'
        numbergroups = 1;
    case 'Two groups mean and comparison (mixed effects)'
        numbergroups = 2;
end
clear answer

names = dir ('pp_*_var.nii'); % identifies the names of the files to use
[numstud, ~] = size(names); % determines how many studies are in the meta-analysis

if numstud == 0 % checks whether any files have been found and if not, returns an error
    error('Cannot find any of the relevant files, please check directory (above) is correct and that you have followed the instructions for running the AES:SDM analysis and unzipping the relevant files.'); % error if can't find map or coordinates
else
end

if exist('sdm_mask.nii', 'file') == 2
else
    error('Cannot find the file "sdm_mask.nii", please check that you have followed the instructions for running the AES:SDM analysis and unzipped "sdm_mask.nii.gz" into the directory.'); % error if can't find sdm_mask.nii
end

if numbergroups == 1 % if specified a single group mean random effects analysis
    
    numberstudies = num2str(numstud); % number of studies as a string
    quest = ['Identified ', numberstudies, ' studies to be included in the analysis. Is this correct?'];
    answer = questdlg(quest,'Check number studies','Yes - continue','No - exit','Yes - continue'); % checks whether this is the correct number of studies
    
    switch answer
        case 'No - exit'
            error(['To correct the number of studies, check there is an unzipped version of pp_NAME_var.nii and a pp_NAME.nii file for each study in the folder ', directory,'. Then try running the script again.']); % error if said wrong number of studies
        case 'Yes - continue'
            disp(['Beginning analysis with ',numberstudies,' studies.']) % message if said right number of studies
    end
    clear answer
    
elseif numbergroups == 2 % if running a 2 group mixed effects analysis extra checks to do:
    
    try
        sdm_table = tdfread('sdm_table.txt'); % checks can read the sdm table which is used to define the groups
    catch
        error('Error reading sdm_table - please check there are the same number of columns and delimiters in each row (check for extra tabs after the final column).')
    end
    
    colnames = fieldnames(sdm_table); % gathers information about the names of colunms in the table
    btntext1 = colnames{end-1,1};
    btntext2 = colnames{end,1};
    
    answer = questdlg('Which column defines the groups of studies?', 'Grouping variable?', btntext1, btntext2, 'Other column', 'Other column');
    switch answer % gets the user to state which column has the coding of groups in
        case btntext1
            colnum = (length(colnames))-1; % assumes most likely columns are those at the end
        case btntext2
            colnum = length(colnames);
        case 'Other column' % option to type in the column title if many columns
            coltext = inputdlg('Please type the name of the column defining the groups exactly as it appears in the sdm table:');
            coltext = coltext{1};
            colnum = find(~cellfun(@isempty,(strfind(colnames, coltext))));
            while length(colnum) ~= 1
                coltext = inputdlg('Please type the name of the column defining the groups exactly as it appears in the sdm table:');
                coltext = coltext{1};
                colnum = find(~cellfun(@isempty,(strfind(colnames, coltext))));
            end
    end
    clear answer
    
    sdm_tab_cell = struct2cell(sdm_table); % converts sdm_table to cell array
    checknumstud = length(sdm_tab_cell{colnum});
    checknumberstudies = num2str(checknumstud); % checks the number of studies according to the sdm_table
    
    if numstud ~= checknumstud % checks this number matches the number of files
        error(['Found ',numberstudies,' image files but only detected ', checknumberstudies, ' rows in the sdm_table. Please ensure these match.'])
    else
        for k=1:numstud
            WithinVarName = names(k).name;
            Name = WithinVarName(4:end-8);
            position = find(~cellfun(@isempty,(strfind(cellstr(sdm_tab_cell{1,1}), Name)))); % for each study a file has been found for, checks can find in the sdm_table so the group can be identified
            if length(position) ~= 1
                error(['Name of image file ', Name, ' not found in sdm_table.']);
            end
        end
        
        groupingvals = unique(sdm_tab_cell{colnum}); % identifies the values used to define the groups - for AES:SDM these should be 0 and 1
        group0ind = find(sdm_tab_cell{colnum} == groupingvals(1));
        group0numstud = length(group0ind);
        group0numberstudies = num2str(group0numstud);
        group1ind = find(sdm_tab_cell{colnum} == groupingvals(2));
        group1numstud = length(group1ind);
        group1numberstudies = num2str(group1numstud); % counts the number of studies in each group to check with the user below
        quest = ['Identified ', checknumberstudies, ' studies to be included in the analysis. ', group0numberstudies, ' in one group and ', group1numberstudies, ' in the other. Is this correct?'];
        answer = questdlg(quest,'Check number studies','Yes - continue','No - exit','Yes - continue'); % checks whether this is the correct number of studies
        switch answer
            case 'No - exit'
                error(['To correct the number of studies, check there is an unzipped version of pp_NAME_var.nii and a pp_NAME.nii file for each study in the folder ', directory,'. Then try running the script again.']); % error if said wrong number of studies
            case 'Yes - continue'
                disp(['Beginning analysis with ',checknumberstudies,' studies. (', group0numberstudies, ' in one group and ', group1numberstudies, ' in the other.)']) % message if said right number of studies
        end
        clear answer
    end
    
end

group0count = 1; % starts index for group 0 if 2 group mixed effects analysis
group1count = 1; % starts index for group 1 if 2 group mixed effects analysis

for k=1:numstud % loop to run for each study
    
    WithinVarName = names(k).name; % sets file name for the within-study variance map
    WithinVarFile = [directory,'\',WithinVarName,',1']; % sets file path for the within-study variance map
    EstimatesName = [WithinVarName(1:end-8),'.nii']; % sets file name for the study effect-size estimates map
    EstimatesFile = [directory,'\',EstimatesName,',1']; % sets file path for the study effect-size estimates map
    Name = WithinVarName(4:end-8); % sets the name of the study
    if numbergroups == 2
        position = find(~cellfun(@isempty,(strfind(cellstr(sdm_tab_cell{1,1}), Name))));
        group = sdm_tab_cell{colnum,1}(position);
    end
    BinarisedName = ['Binarised_',EstimatesName]; % sets the file name for the binarised map of where there is data (created below)
    BinarisedFile = [directory,'\',BinarisedName,',1']; % sets the file path for the binarised map of where there is data (created below)
    FE_WtName = ['FE_Wt_',Name]; % sets file name for the within-study inverse variance map (created below)
    FE_WtFile = [directory,'\',FE_WtName,'.nii,1']; % sets file path for within-study inverse variance map (created below)
    FEWt_x_EstName = ['FEWt_x_Estimates_',Name]; % sets file name for the study weighted effect-size estimates map (created below)
    FEWt_x_EstFile = [directory,'\',FEWt_x_EstName,'.nii,1']; % set file path for within-study weighted effect-size estimates
    FE_Adj_WtName = ['FE_Adj_Wt_',Name]; % sets file name for the adjusted within study inverse variance map (FE_Wt * Bin - created below)
    FE_Adj_WtFile = [directory,'\',FE_Adj_WtName,'.nii,1']; % set file path for within-study adjusted inverse variance
    FE_All_Wt_SqName = ['FE_All_Wt_Sq_',Name]; % sets file name for fixed effects inverse variance squared for all data (created below)
    FE_All_Wt_SqFile = [directory,'\',FE_All_Wt_SqName,'.nii,1'];  % set file path for fixed effects inverse variance squared for original analysis (with all data; created below)
    FE_Adj_Wt_SqName = ['FE_Adj_Wt_Sq_',Name]; % sets file name for fixed effects inverse variance squared for adjusted analysis (created below)
    FE_Adj_Wt_SqFile = [directory,'\',FE_Adj_Wt_SqName,'.nii,1']; % set file path for fixed effects inverse variance squared for adjusted analysis (created below)
    FEWt_x_EstSqName = ['FEWt_x_EstimatesSq_',Name];  % sets file name for effect size map squared multipled by the fixed effects inverse variance (created below)
    FEWt_x_EstSqFile = [directory,'\',FEWt_x_EstSqName,'.nii,1'];  % sets file path for effect-size estimates map squared multipled by the fixed effects inverse variance (created below)
    BinarisedList{k,1}=BinarisedFile; % create a list of binarised map file paths
    
    if numbergroups == 1
        FE_WtList{k,1}=FE_WtFile; % create a list of within-study inverse variance file paths
        FEWt_x_EstList{k,1}=FEWt_x_EstFile; % create a list of within-study weighted effect-size estimates file paths
        FE_Adj_WtList{k,1}=FE_Adj_WtFile; % create a list of within-study adjusted inverse variance file paths
        FE_All_Wt_SqList{k,1}=FE_All_Wt_SqFile; % create a list of fixed effects inverse variance squared (original) file paths
        FE_Adj_Wt_SqList{k,1}=FE_Adj_Wt_SqFile; % create a list of fixed effects inverse variance squared (adjusted) file paths
        FEWt_x_EstSqList{k,1}=FEWt_x_EstSqFile; % create a list of effect-size estimates map squared multipled by the fixed effects inverse variance file paths
    elseif numbergroups == 2 && group == 0 % same as above but only for studies in group 0 if 2 group mixed effects analysis
        BinarisedListGroup0{group0count,1}=BinarisedFile;
        FE_WtListGroup0{group0count,1}=FE_WtFile;
        FEWt_x_EstListGroup0{group0count,1}=FEWt_x_EstFile;
        FE_Adj_WtListGroup0{group0count,1}=FE_Adj_WtFile;
        FE_All_Wt_SqListGroup0{group0count,1}=FE_All_Wt_SqFile;
        FE_Adj_Wt_SqListGroup0{group0count,1}=FE_Adj_Wt_SqFile;
        FEWt_x_EstSqListGroup0{group0count,1}=FEWt_x_EstSqFile;
        group0count = group0count + 1;
    elseif numbergroups == 2 && group == 1 % same as above but only for studies in group 1 if 2 group mixed effects analysis
        BinarisedListGroup1{group1count,1}=BinarisedFile; 
        FE_WtListGroup1{group1count,1}=FE_WtFile;
        FEWt_x_EstListGroup1{group1count,1}=FEWt_x_EstFile;
        FE_Adj_WtListGroup1{group1count,1}=FE_Adj_WtFile;
        FE_All_Wt_SqListGroup1{group1count,1}=FE_All_Wt_SqFile;
        FE_Adj_Wt_SqListGroup1{group1count,1}=FE_Adj_Wt_SqFile;
        FEWt_x_EstSqListGroup1{group1count,1}=FEWt_x_EstSqFile;
        group1count = group1count + 1;
    end
    
    GzName = [directory,'\',Name,'.nii.gz']; % sets file name for raw data map
    if exist(GzName, 'file') == 2 % checks whether this file exists - if it does, a map is available - if not, will check for coordinates
        InputFiles = {EstimatesFile};
        run_imcalc(InputFiles, BinarisedName, directory, '(i1>0.00001 | i1<-0.00001)') % create map of where data exists if a map
    else
        TextName = [directory,'\',Name,'.*.txt']; % if no map, check for coordinates text file
        if ~isempty(dir(TextName)) == 1 % if there is a coordinates file,
            copyfile ('sdm_mask.nii', BinarisedName) % set the binarised file to the whole mask (coordinates not adjusted as zeros are meaningful)
        else
            error(['Cannot determine whether file ',Name,' is an image or text file']); % error if can't find map or coordinates
        end
    end
    
    InputFiles = {WithinVarFile};
    run_imcalc(InputFiles, FE_WtName, directory, '1./i1') % calculate the inverse of the within-study variance map
    
    InputFiles = {
        EstimatesFile
        FE_WtFile
        };
    run_imcalc(InputFiles, FEWt_x_EstName, directory, 'i1.*i2') % weight (multiply) the effect-size estimates by the inverse of the variance
    
    InputFiles = {
        BinarisedFile
        FE_WtFile
        };
    run_imcalc(InputFiles, FE_Adj_WtName, directory, 'i1.*i2') % adjust the inverse variance so values are only present where there is data
    
    InputFiles = {
        FE_WtFile
        };
    run_imcalc(InputFiles, FE_All_Wt_SqName, directory, 'i1.^2') % square the fixed effects (FE) weighting (inverse variance)
    
    InputFiles = {
        FE_All_Wt_SqFile
        BinarisedFile
        };
    run_imcalc(InputFiles, FE_Adj_Wt_SqName, directory, 'i1.*i2') % adjust the squared FE weighting so only where data present
    
    InputFiles = {
        EstimatesFile
        FE_WtFile
        };
    run_imcalc(InputFiles, FEWt_x_EstSqName, directory, '(i1.^2).*i2') % multiply the effect-size estimates by the FE weighting
    
end

sumcalc = 'i1 + i'; % creates a string to use for summing maps
index = 2;

for i=2:numstud % creates string to sum for the correct number of studies
    next = num2str(index);
    if i == numstud
        sumcalc = [sumcalc, next];
    else
        sumcalc = [sumcalc, next, ' + i'];
    end
    index = index + 1;
end

if numbergroups == 1 % if doing a fixed effects analysis on one group, runs the relevant calculations
    
    InputFiles = FEWt_x_EstList;
    run_imcalc(InputFiles, 'Sum_FE_Wt_x_Estimates', directory, sumcalc) % sum the FE effect-size estimates
    
    InputFiles = FE_WtList;
    run_imcalc(InputFiles, 'Sum_FE_All_Wt', directory, sumcalc) % sum the FE weightings (original)
    
    InputFiles = FE_Adj_WtList;
    run_imcalc(InputFiles, 'Sum_FE_Adj_Wt', directory, sumcalc) % sum the FE weightings (adjusted)
    
    Sum_FE_Wt_x_EstimatesFile = [directory,'\','Sum_FE_Wt_x_Estimates.nii,1']; % sets file path for the sum of FE effect-size estimates maps
    Sum_FE_All_WtFile = [directory,'\','Sum_FE_All_Wt.nii,1']; % sets file path for the sum of FE weightings (original)
    Sum_FE_Adj_WtFile = [directory,'\','Sum_FE_Adj_Wt.nii,1']; % sets file path for the sum of FE weightings (adjusted)
    
    InputFiles = {
        Sum_FE_Wt_x_EstimatesFile
        Sum_FE_All_WtFile
        };
    run_imcalc(InputFiles, 'FE_All_EffectSize', directory, 'i1./i2') % calculate the FE meta-analysed effect size map (original)
    
    InputFiles = {
        Sum_FE_Wt_x_EstimatesFile
        Sum_FE_Adj_WtFile
        };
    run_imcalc(InputFiles, 'FE_Adj_EffectSize', directory, 'i1./i2') % calcuate the FE meta-analysed effect size map (adjusted)
    
    InputFiles = {
        Sum_FE_All_WtFile
        };
    run_imcalc(InputFiles, 'FE_All_Variance', directory, '1./i1') % calculate the FE meta-analysed variance map (original)
    
    InputFiles = {
        Sum_FE_Adj_WtFile
        };
    run_imcalc(InputFiles, 'FE_Adj_Variance', directory, '1./i1') % calculate the FE meta-analysed variance map (adjusted)
    
    FE_All_VarianceFile = [directory,'\','FE_All_Variance.nii,1']; % sets file path for the map of FE variance (original)
    FE_Adj_VarianceFile = [directory,'\','FE_Adj_Variance.nii,1']; % sets file path for the map of FE variance (adjusted)
    
    InputFiles = {
        FE_All_VarianceFile
        };
    run_imcalc(InputFiles, 'FE_All_StandardError', directory, 'sqrt(i1)') % calculate the FE meta-analysed standard error map (original)
    
    InputFiles = {
        FE_Adj_VarianceFile
        };
    run_imcalc(InputFiles, 'FE_Adj_StandardError', directory, 'sqrt(i1)') % calculate the FE meta-analysed standard error map (adjusted)
    
    FE_All_StandardErrorFile = [directory,'\','FE_All_StandardError.nii,1']; % sets file path for the map of FE standard error (original)
    FE_Adj_StandardErrorFile = [directory,'\','FE_Adj_StandardError.nii,1']; % sets file path for the map of FE standard error (adjusted)
    FE_All_EffectSizeFile = [directory,'\','FE_All_EffectSize.nii,1']; % sets file path for the map of FE effect sizes (original)
    FE_Adj_EffectSizeFile = [directory,'\','FE_Adj_EffectSize.nii,1']; % sets file path for the map of FE effect sizes (adjusted)
    
    InputFiles = {
        FE_All_EffectSizeFile
        FE_All_StandardErrorFile
        };
    run_imcalc(InputFiles, 'FE_All_Z', directory, 'i1./i2') % calculate the FE meta-analysed SDM-Z score map (original)
    
    InputFiles = {
        FE_Adj_EffectSizeFile
        FE_Adj_StandardErrorFile
        };
    run_imcalc(InputFiles, 'FE_Adj_Z', directory, 'i1./i2') % calculate the FE meta-analysed SDM-Z score map (adjusted)
    
    InputFiles = FE_All_Wt_SqList;
    run_imcalc(InputFiles, 'Sum_FE_All_Wt_Sq', directory, sumcalc); % sum the squared FE weightings (original)
    
    InputFiles = FE_Adj_Wt_SqList;
    run_imcalc(InputFiles, 'Sum_FE_Adj_Wt_Sq', directory, sumcalc); % sum the squared FE weightings (adjusted)
    
    InputFiles = FEWt_x_EstSqList;
    run_imcalc(InputFiles, 'Sum_FEWt_x_EstSq', directory, sumcalc) % sum the squared effect-size estimates weighted with FE inverse variance
    
    Sum_FE_All_WtFile = [directory,'\Sum_FE_All_Wt.nii,1']; % sets file path for the sum of weightings (original)
    Sum_FE_Adj_WtFile = [directory,'\Sum_FE_Adj_Wt.nii,1']; % sets file path for the sum of weightings (adjusted)
    Sum_FEWt_x_EstSqFile = [directory,'\Sum_FEWt_x_EstSq.nii,1']; % sets file path for the sum of weightings multiplied by squared effect-size estimates (adjusted)
    Sum_FE_All_Wt_SqFile = [directory,'\Sum_FE_All_Wt_Sq.nii,1']; % sets file path for the sum of weightings squared (adjusted)
    Sum_FE_Adj_Wt_SqFile = [directory,'\Sum_FE_Adj_Wt_Sq.nii,1']; % sets file path for the sum of weightings squared (adjusted)
    
    InputFiles = {
        Sum_FEWt_x_EstSqFile
        Sum_FE_Wt_x_EstimatesFile
        Sum_FE_All_WtFile
        };
    run_imcalc(InputFiles, 'RE_All_Q', directory, 'i1-((i2.^2)./i3)') % calculate map of Q for the original random effects
    
    InputFiles = {
        Sum_FEWt_x_EstSqFile
        Sum_FE_Wt_x_EstimatesFile
        Sum_FE_Adj_WtFile
        };
    run_imcalc(InputFiles, 'RE_Adj_Q', directory, 'i1-((i2.^2)./i3)') % calculate map of Q for the adjusted random effects
    
    InputFiles = {
        Sum_FE_All_WtFile
        Sum_FE_All_Wt_SqFile
        };
    run_imcalc(InputFiles, 'RE_All_C', directory, 'i1-i2./i1') % calculate map of C for the original random effects
    
    InputFiles = {
        Sum_FE_Adj_WtFile
        Sum_FE_Adj_Wt_SqFile
        };
    run_imcalc(InputFiles, 'RE_Adj_C', directory, 'i1-i2./i1') % calculate map of C for the adjusted random effects
    
    REAllDoFcalc = ['i1.*',num2str(numstud-1)];
    
    InputFiles = {[directory,'\sdm_mask.nii']};
    run_imcalc(InputFiles, 'RE_All_DoF', directory, REAllDoFcalc) % calculate the degrees of freedom for each voxel in the original analysis
    
    REAdjDoFcalc = ['(',sumcalc,')-1'];
    
    InputFiles = BinarisedList;
    run_imcalc(InputFiles, 'RE_Adj_DoF', directory, REAdjDoFcalc) % calculate the degrees of freedom for each voxel in the adjusted analysis
    
    RE_All_QFile = [directory,'\RE_All_Q.nii,1']; % sets file path for the map of Q (original)
    RE_All_DoFFile = [directory,'\RE_All_DoF.nii,1']; % sets file path for the degrees of freedom map (original)
    RE_Adj_QFile = [directory,'\RE_Adj_Q.nii,1']; % sets file path for the map of Q (adjusted)
    RE_Adj_DoFFile = [directory,'\RE_Adj_DoF.nii,1']; % sets file path for the degrees of freedom map (adjusted)
    
    numcalc = 'max((i1-i2),0)';
    
    InputFiles = {
        RE_All_QFile
        RE_All_DoFFile
        };
    run_imcalc(InputFiles, 'RE_All_Numerator', directory, numcalc) % calculate map of the numerator for the original random effects
    
    InputFiles = {
        RE_Adj_QFile
        RE_Adj_DoFFile
        };
    run_imcalc(InputFiles, 'RE_Adj_Numerator', directory, numcalc) % calculate map of the numerator for the adjusted random effects
    
    RE_All_CFile = [directory,'\RE_All_C.nii,1']; % sets file path for the map of C (original)
    RE_All_NumeratorFile = [directory,'\RE_All_Numerator.nii,1']; % sets file path for the map of the numerator (original)
    RE_Adj_CFile = [directory,'\RE_Adj_C.nii,1']; % sets file path for the map of C (adjusted)
    RE_Adj_NumeratorFile = [directory,'\RE_Adj_Numerator.nii,1']; % sets file path for the map of the numerator (adjusted)
    
    InputFiles = {
        RE_All_NumeratorFile
        RE_All_CFile
        };
    run_imcalc(InputFiles, 'RE_All_TauSq', directory, 'i1./i2') % calculate map of tau squared for the original random effects
    
    RE_All_TauSqFile = [directory,'\RE_All_TauSq.nii,1']; % sets file path for the map of tau squared (original)
    
    InputFiles = {
        RE_Adj_NumeratorFile
        RE_Adj_CFile
        };
    run_imcalc(InputFiles, 'RE_Adj_TauSq', directory, 'i1./i2') % calculate map of tau squared for the adjusted random effects
    
    RE_Adj_TauSqFile = [directory,'\RE_Adj_TauSq.nii,1']; % sets file path for the map of tau squared (adjusted)
    
    %%
    
    for k = 1:numstud
        
        WithinVarName = names(k).name; % sets file name for the within-study variance map
        WithinVarFile = [directory,'\',WithinVarName,',1']; % sets file path for the within-study variance map
        EstimatesName = [WithinVarName(1:end-8),'.nii']; % sets file name for the study effect size map
        EstimatesFile = [directory,'\',EstimatesName,',1']; % sets file path for the study effect size map
        Name = WithinVarName(4:end-8); % sets the name of the study
        BinarisedName = ['Binarised_',EstimatesName]; % sets  file name for the binarised map of where there is data (created below)
        BinarisedFile = [directory,'\',BinarisedName,',1']; % sets file path for the binarised map of where there is data (created below)
        RE_All_Wt_x_EstName = ['RE_All_Wt_x_Est_',Name]; % sets file name for the random effects weighting multipled by the effect-size estimate (original)
        RE_All_WtName = ['RE_All_Wt_',Name]; % sets file name for the random effects weighting (original)
        RE_All_WtFile = [directory,'\',RE_All_WtName,'.nii,1']; % sets file path for the random effects weighting (original)
        RE_Adj_Wt_x_EstName = ['RE_Adj_Wt_x_Est_',Name]; % sets file name for the random effects weighting multipled by the effect-size estimate (adjusted)
        RE_Adj_WtName = ['RE_Adj_Wt_',Name]; % sets file name for the random effects weighting (adjusted)
        RE_Adj_WtFile = [directory,'\',RE_Adj_WtName,'.nii,1']; % sets file path for the random effects weighting (adjusted)
        
        InputFiles = {
            RE_All_TauSqFile
            WithinVarFile
            };
        run_imcalc(InputFiles, RE_All_WtName, directory, '1./(i1+i2)') % calculates total inverse variance (RE weighting, original analysis) for each study by 1 / (between + within study variances)
        
        InputFiles = {
            RE_Adj_TauSqFile
            WithinVarFile
            BinarisedFile
            };
        run_imcalc(InputFiles, RE_Adj_WtName, directory, '(1./(i1+i2)).*i3') % calculates total inverse variance (RE weighting, adjusted analysis) for each study by 1 / (between + within study variances)
        
        InputFiles = {
            EstimatesFile
            RE_All_WtFile
            };
        run_imcalc(InputFiles, RE_All_Wt_x_EstName, directory, 'i1.*i2') % multiples the weighting (original) by the effect size estimate
        
        InputFiles = {
            EstimatesFile
            RE_Adj_WtFile
            };
        run_imcalc(InputFiles, RE_Adj_Wt_x_EstName, directory, 'i1.*i2') % multiples the weighting (adjusted) by the effect size estimate
        
        RE_All_Wt_x_EstFile = [directory,'\',RE_All_Wt_x_EstName,'.nii,1']; % sets the file paths and creates lists of these for the files just created
        RE_All_Wt_x_EstList{k,1}=RE_All_Wt_x_EstFile;
        RE_All_WtFile = [directory,'\',RE_All_WtName,'.nii,1'];
        RE_All_WtList{k,1}=RE_All_WtFile;
        
        RE_Adj_Wt_x_EstFile = [directory,'\',RE_Adj_Wt_x_EstName,'.nii,1'];
        RE_Adj_Wt_x_EstList{k,1}=RE_Adj_Wt_x_EstFile;
        RE_Adj_WtFile = [directory,'\',RE_Adj_WtName,'.nii,1'];
        RE_Adj_WtList{k,1}=RE_Adj_WtFile;
        
    end
    
    %%
    
    InputFiles = RE_All_Wt_x_EstList;
    run_imcalc(InputFiles, 'Sum_RE_All_Estimates', directory, sumcalc) % sums the effect size estimates (original analysis)
    
    InputFiles = RE_Adj_Wt_x_EstList;
    run_imcalc(InputFiles, 'Sum_RE_Adj_Estimates', directory, sumcalc) % sums the effect size estimates (adjusted analysis)
    
    InputFiles = RE_All_WtList;
    run_imcalc(InputFiles, 'Sum_RE_All_Wt', directory, sumcalc);  % sums the weightings (original analysis)
    
    InputFiles = RE_Adj_WtList;
    run_imcalc(InputFiles, 'Sum_RE_Adj_Wt', directory, sumcalc); % sums the weightings (adjusted analysis)
    
    Sum_RE_All_EstimatesFile = [directory,'\','Sum_RE_All_Estimates.nii,1']; % sets the file paths for the files just created
    Sum_RE_All_WtFile = [directory,'\','Sum_RE_All_Wt.nii,1'];
    Sum_RE_Adj_EstimatesFile = [directory,'\','Sum_RE_Adj_Estimates.nii,1'];
    Sum_RE_Adj_WtFile = [directory,'\','Sum_RE_Adj_Wt.nii,1'];
    
    InputFiles = {
        Sum_RE_All_EstimatesFile
        Sum_RE_All_WtFile
        };
    run_imcalc(InputFiles, 'RE_All_EffectSize', directory, 'i1./i2') % divides sum of the estimates by the sum of the weightings to give the random effects effect size image for original analysis
    
    InputFiles = {
        Sum_RE_Adj_EstimatesFile
        Sum_RE_Adj_WtFile
        };
    run_imcalc(InputFiles, 'RE_Adj_EffectSize', directory, 'i1./i2') % divides sum of the estimates by the sum of the weightings to give the random effects effect size image for adjusted analysis
    
    InputFiles = {Sum_RE_All_WtFile};
    run_imcalc(InputFiles, 'RE_All_Variance', directory, '1./i1') % calculates the variance for the original analysis (1 / sum of weightings)
    
    InputFiles = {Sum_RE_Adj_WtFile};
    run_imcalc(InputFiles, 'RE_Adj_Variance', directory, '1./i1') % calculates the variance for the adjusted analysis (1 / sum of weightings)
    
    RE_All_VarianceFile = [directory,'\','RE_All_Variance.nii,1']; % sets file paths
    RE_Adj_VarianceFile = [directory,'\','RE_Adj_Variance.nii,1'];
    
    RE_var_neg_vals = 1;
    for x = 1:2 % checks whether any voxels in the variance maps have negative values which may happen in the adjusted analysis
        tocheck = {'RE_All_Variance.nii', 'RE_Adj_Variance.nii'};
        varfile = tocheck{x};
        readin = spm_vol(varfile);
        getfile = spm_read_vols(readin);
        minimum = min(min(min(getfile)));
        if minimum < 0 % any negative values need to be removed before conversion to standard error (square root of variance)
            copyname = [varfile(1:end-4),'_NegValues.nii'];
            copyfile(varfile, copyname);
            NegValuesFile = [directory,'\',copyname,',1']; % creates a copy of the original variance map for the user to investigate which voxels the issues were in
            InputFiles = {NegValuesFile};
            run_imcalc(InputFiles, varfile(1:end-4), directory, 'i1.*(i1>0)')
            RE_var_neg_vals_log{RE_var_neg_vals,1} = varfile; % keeps a log of any files which have been edited to report at the end
            RE_var_neg_vals = RE_var_neg_vals + 1;
        else
        end
    end
    
    InputFiles = {RE_All_VarianceFile};
    run_imcalc(InputFiles, 'RE_All_StandardError', directory, 'sqrt(i1)') % converts variance to standard error for the original analysis
    
    InputFiles = {RE_Adj_VarianceFile};
    run_imcalc(InputFiles, 'RE_Adj_StandardError', directory, 'sqrt(i1)') % converts variance to standard error for the adjusted analysis
    
    RE_All_StandardErrorFile = [directory,'\','RE_All_StandardError.nii,1']; % sets file paths
    RE_All_EffectSizeFile = [directory,'\','RE_All_EffectSize.nii,1'];
    RE_Adj_StandardErrorFile = [directory,'\','RE_Adj_StandardError.nii,1'];
    RE_Adj_EffectSizeFile = [directory,'\','RE_Adj_EffectSize.nii,1'];
    
    InputFiles = {
        RE_All_EffectSizeFile
        RE_All_StandardErrorFile
        };
    run_imcalc(InputFiles, 'RE_All_Z', directory, 'i1./i2') % calculates the map of Z scores for the random effects, original analysis
    
    InputFiles = {
        RE_Adj_EffectSizeFile
        RE_Adj_StandardErrorFile
        };
    run_imcalc(InputFiles, 'RE_Adj_Z', directory, 'i1./i2') % calculates the map of Z scores for the random effects, adjusted analysis
    
    disp('Random effects analysis for single group mean complete.')
    
    if exist('RE_var_neg_vals_log') == 1 % if any variance files were changed due to negative values, reports them here
        disp('IMPORTANT - the following variance files contained negative values so were masked to exclude these before standard error could be calculated. See files *_Variance_*_NegValues.nii for the location of these voxels.')
        celldisp(RE_var_neg_vals_log)
    else
    end
    
    %% MIXED EFFECTS
    
elseif numbergroups == 2
    
    sumcalcGroup0 = 'i1 + i'; % creates a string to use for summing maps in group 0
    index = 2;
    
    for i=2:group0numstud % creates string to sum for the correct number of studies
        next = num2str(index);
        if i == group0numstud
            sumcalcGroup0 = [sumcalcGroup0, next];
        else
            sumcalcGroup0 = [sumcalcGroup0, next, ' + i'];
        end
        index = index + 1;
    end
    
    sumcalcGroup1 = 'i1 + i'; % creates a string to use for summing maps in group 1
    index = 2;
    
    for i=2:group1numstud % creates string to sum for the correct number of studies
        next = num2str(index);
        if i == group1numstud
            sumcalcGroup1 = [sumcalcGroup1, next];
        else
            sumcalcGroup1 = [sumcalcGroup1, next, ' + i'];
        end
        index = index + 1;
    end
    
    InputFiles = FEWt_x_EstListGroup0;
    run_imcalc(InputFiles, 'Sum_FE_Wt_x_Estimates_Group0', directory, sumcalcGroup0) % sum the FE effect-size estimates for group 0
    
    InputFiles = FEWt_x_EstListGroup1;
    run_imcalc(InputFiles, 'Sum_FE_Wt_x_Estimates_Group1', directory, sumcalcGroup1) % sum the FE effect-size estimates for group 1
    
    InputFiles = FE_WtListGroup0;
    run_imcalc(InputFiles, 'Sum_FE_All_Wt_Group0', directory, sumcalcGroup0) % sum the FE weightings (original) for group 0
    
    InputFiles = FE_WtListGroup1;
    run_imcalc(InputFiles, 'Sum_FE_All_Wt_Group1', directory, sumcalcGroup1) % sum the FE weightings (original) for group 1
    
    InputFiles = FE_Adj_WtListGroup0;
    run_imcalc(InputFiles, 'Sum_FE_Adj_Wt_Group0', directory, sumcalcGroup0) % sum the FE weightings (adjusted) for group 0
    
    InputFiles = FE_Adj_WtListGroup1;
    run_imcalc(InputFiles, 'Sum_FE_Adj_Wt_Group1', directory, sumcalcGroup1) % sum the FE weightings (adjusted) for group 1
    
    Sum_FE_Wt_x_EstimatesFileGroup0 = [directory,'\','Sum_FE_Wt_x_Estimates_Group0.nii,1']; % sets file path for the sum of FE effect-size estimates maps for group 0
    Sum_FE_Wt_x_EstimatesFileGroup1 = [directory,'\','Sum_FE_Wt_x_Estimates_Group1.nii,1']; % sets file path for the sum of FE effect-size estimates maps for group 1
    Sum_FE_All_WtFileGroup0 = [directory,'\Sum_FE_All_Wt_Group0.nii,1']; % sets file path for the sum of weightings (original) for group 0
    Sum_FE_Adj_WtFileGroup0 = [directory,'\Sum_FE_Adj_Wt_Group0.nii,1']; % sets file path for the sum of weightings (adjusted) for group 0
    Sum_FEWt_x_EstSqFileGroup0 = [directory,'\Sum_FEWt_x_EstSq_Group0.nii,1']; % sets file path for the sum of weightings multiplied by squared effect-size estimates (adjusted) for group 0
    Sum_FE_All_Wt_SqFileGroup0 = [directory,'\Sum_FE_All_Wt_Sq_Group0.nii,1']; % sets file path for the sum of weightings squared (adjusted) for group 0
    Sum_FE_Adj_Wt_SqFileGroup0 = [directory,'\Sum_FE_Adj_Wt_Sq_Group0.nii,1']; % sets file path for the sum of weightings squared (adjusted) for group 0
    
    Sum_FE_All_WtFileGroup1 = [directory,'\Sum_FE_All_Wt_Group1.nii,1']; % sets file path for the sum of weightings (original) for group 1
    Sum_FE_Adj_WtFileGroup1 = [directory,'\Sum_FE_Adj_Wt_Group1.nii,1']; % sets file path for the sum of weightings (adjusted) for group 1
    Sum_FEWt_x_EstSqFileGroup1 = [directory,'\Sum_FEWt_x_EstSq_Group1.nii,1']; % sets file path for the sum of weightings multiplied by squared effect-size estimates (adjusted) for group 1
    Sum_FE_All_Wt_SqFileGroup1 = [directory,'\Sum_FE_All_Wt_Sq_Group1.nii,1']; % sets file path for the sum of weightings squared (adjusted) for group 1
    Sum_FE_Adj_Wt_SqFileGroup1 = [directory,'\Sum_FE_Adj_Wt_Sq_Group1.nii,1']; % sets file path for the sum of weightings squared (adjusted) for group 1  
    
    InputFiles = FE_All_Wt_SqListGroup0;
    run_imcalc(InputFiles, 'Sum_FE_All_Wt_Sq_Group0', directory, sumcalcGroup0); % sum the squared FE weightings (original) for group 0
    
    InputFiles = FE_All_Wt_SqListGroup1;
    run_imcalc(InputFiles, 'Sum_FE_All_Wt_Sq_Group1', directory, sumcalcGroup1); % sum the squared FE weightings (original) for group 1
    
    InputFiles = FE_Adj_Wt_SqListGroup0;
    run_imcalc(InputFiles, 'Sum_FE_Adj_Wt_Sq_Group0', directory, sumcalcGroup0); % sum the squared FE weightings (adjusted) for group 0
    
    InputFiles = FE_Adj_Wt_SqListGroup1;
    run_imcalc(InputFiles, 'Sum_FE_Adj_Wt_Sq_Group1', directory, sumcalcGroup1); % sum the squared FE weightings (adjusted) for group 1
    
    InputFiles = FEWt_x_EstSqListGroup0;
    run_imcalc(InputFiles, 'Sum_FEWt_x_EstSq_Group0', directory, sumcalcGroup0) % sum the squared effect-size estimates weighted with FE inverse variance for group 0
    
    InputFiles = FEWt_x_EstSqListGroup1;
    run_imcalc(InputFiles, 'Sum_FEWt_x_EstSq_Group1', directory, sumcalcGroup1) % sum the squared effect-size estimates weighted with FE inverse variance for group 1
    
    InputFiles = {
        Sum_FEWt_x_EstSqFileGroup0
        Sum_FE_Wt_x_EstimatesFileGroup0
        Sum_FE_All_WtFileGroup0
        };
    run_imcalc(InputFiles, 'ME_All_Q_Group0', directory, 'i1-((i2.^2)./i3)') % calculate map of Q for the original mixed effects for group 0
    
    InputFiles = {
        Sum_FEWt_x_EstSqFileGroup1
        Sum_FE_Wt_x_EstimatesFileGroup1
        Sum_FE_All_WtFileGroup1
        };
    run_imcalc(InputFiles, 'ME_All_Q_Group1', directory, 'i1-((i2.^2)./i3)') % calculate map of Q for the original mixed effects for group 1
    
    InputFiles = {
        Sum_FEWt_x_EstSqFileGroup0
        Sum_FE_Wt_x_EstimatesFileGroup0
        Sum_FE_Adj_WtFileGroup0
        };
    run_imcalc(InputFiles, 'ME_Adj_Q_Group0', directory, 'i1-((i2.^2)./i3)') % calculate map of Q for the adjusted mixed effects for group 0
    
    InputFiles = {
        Sum_FEWt_x_EstSqFileGroup1
        Sum_FE_Wt_x_EstimatesFileGroup1
        Sum_FE_Adj_WtFileGroup1
        };
    run_imcalc(InputFiles, 'ME_Adj_Q_Group1', directory, 'i1-((i2.^2)./i3)') % calculate map of Q for the adjusted mixed effects for group 1
    
    InputFiles = {
        Sum_FE_All_WtFileGroup0
        Sum_FE_All_Wt_SqFileGroup0
        };
    run_imcalc(InputFiles, 'ME_All_C_Group0', directory, 'i1-i2./i1') % calculate map of C for the original mixed effects for group 0
    
    InputFiles = {
        Sum_FE_All_WtFileGroup1
        Sum_FE_All_Wt_SqFileGroup1
        };
    run_imcalc(InputFiles, 'ME_All_C_Group1', directory, 'i1-i2./i1') % calculate map of C for the original mixed effects for group 1
    
    InputFiles = {
        Sum_FE_Adj_WtFileGroup0
        Sum_FE_Adj_Wt_SqFileGroup0
        };
    run_imcalc(InputFiles, 'ME_Adj_C_Group0', directory, 'i1-i2./i1') % calculate map of C for the adjusted mixed effects for group 0
    
    InputFiles = {
        Sum_FE_Adj_WtFileGroup1
        Sum_FE_Adj_Wt_SqFileGroup1
        };
    run_imcalc(InputFiles, 'ME_Adj_C_Group1', directory, 'i1-i2./i1') % calculate map of C for the adjusted mixed effects for group 1
    
    MEAllDoFcalc = ['i1.*',num2str(numstud-2)];
    
    InputFiles = {[directory,'\sdm_mask.nii']};
    run_imcalc(InputFiles, 'ME_All_DoF', directory, MEAllDoFcalc) % calculate the degrees of freedom for each voxel in the original analysis
    
    MEAdjDoFcalc = ['(',sumcalc,')-2'];
    
    InputFiles = BinarisedList;
    run_imcalc(InputFiles, 'ME_Adj_DoF', directory, MEAdjDoFcalc) % calculate the degrees of freedom for each voxel in the adjusted analysis
    
    ME_All_QFileGroup0 = [directory,'\ME_All_Q_Group0.nii,1']; % sets file path for the map of Q (original) for group 0
    ME_All_QFileGroup1 = [directory,'\ME_All_Q_Group1.nii,1']; % sets file path for the map of Q (original) for group 1
    ME_All_QFile = [directory,'\ME_All_Q.nii,1']; % sets file path for the map of Q (original) for groups 0 and 1 together
    ME_All_DoFFile = [directory,'\ME_All_DoF.nii,1']; % sets file path for the degrees of freedom map (original)
    ME_All_CFileGroup0 = [directory,'\ME_All_C_Group0.nii,1']; % sets file path for the map of C (original) for group 0
    ME_All_CFileGroup1 = [directory,'\ME_All_C_Group1.nii,1']; % sets file path for the map of C (original) for group 1
    ME_All_CFile = [directory,'\ME_All_C.nii,1']; % sets file path for the map of C (original) for groups 0 and 1 together
    ME_Adj_QFileGroup0 = [directory,'\ME_Adj_Q_Group0.nii,1']; % sets file path for the map of Q (adjusted) for group 0
    ME_Adj_QFileGroup1 = [directory,'\ME_Adj_Q_Group1.nii,1']; % sets file path for the map of Q (adjusted) for group 1
    ME_Adj_QFile = [directory,'\ME_Adj_Q.nii,1']; % sets file path for the map of Q (adjusted) for groups 0 and 1 together
    ME_Adj_DoFFile = [directory,'\ME_Adj_DoF.nii,1']; % sets file path for the degrees of freedom map (adjusted)
    ME_Adj_CFileGroup0 = [directory,'\ME_Adj_C_Group0.nii,1']; % sets file path for the map of C (adjusted) for group 0
    ME_Adj_CFileGroup1 = [directory,'\ME_Adj_C_Group1.nii,1']; % sets file path for the map of C (adjusted) for group 1
    ME_Adj_CFile = [directory,'\ME_Adj_C.nii,1']; % sets file path for the map of C (adjusted) for groups 0 and 1 together
    
    InputFiles = {
        ME_All_QFileGroup0
        ME_All_QFileGroup1
        };
    run_imcalc(InputFiles, 'ME_All_Q', directory, 'i1+i2') % calculate map of Q for the original mixed effects analysis for group 0 and 1
    
    InputFiles = {
        ME_Adj_QFileGroup0
        ME_Adj_QFileGroup1
        };
    run_imcalc(InputFiles, 'ME_Adj_Q', directory, 'i1+i2') % calculate map of Q for the adjusted mixed effects analysis for group 0 and 1
    
    InputFiles = {
        ME_All_CFileGroup0
        ME_All_CFileGroup1
        };
    run_imcalc(InputFiles, 'ME_All_C', directory, 'i1+i2') % calculate map of C for the original mixed effects analysis for group 0 and 1
    
    InputFiles = {
        ME_Adj_CFileGroup0
        ME_Adj_CFileGroup1
        };
    run_imcalc(InputFiles, 'ME_Adj_C', directory, 'i1+i2') % calculate map of C for the adjusted mixed effects analysis for group 0 and 1
    
    numcalc = 'max((i1-i2),0)';
    
    InputFiles = {
        ME_All_QFile
        ME_All_DoFFile
        };
    run_imcalc(InputFiles, 'ME_All_Numerator', directory, numcalc) % calculate map of the numerator based on the original analysis
    
    InputFiles = {
        ME_Adj_QFile
        ME_Adj_DoFFile
        };
    run_imcalc(InputFiles, 'ME_Adj_Numerator', directory, numcalc) % calculate map of the numerator based on the adjusted analysis
    
    ME_All_NumeratorFile = [directory,'\ME_All_Numerator.nii,1']; % sets file path for the map of the numerator (original)
    ME_Adj_NumeratorFile = [directory,'\ME_Adj_Numerator.nii,1']; % sets file path for the map of the numerator (adjusted)
    
    InputFiles = {
        ME_All_NumeratorFile
        ME_All_CFile
        };
    run_imcalc(InputFiles, 'ME_All_TauSq', directory, 'i1./i2') % calculate map of mixed effects tau squared for the original analysis
    
    ME_All_TauSqFile = [directory,'\ME_All_TauSq.nii,1']; % sets file path for the map of tau squared (original)
    
    InputFiles = {
        ME_Adj_NumeratorFile
        ME_Adj_CFile
        };
    run_imcalc(InputFiles, 'ME_Adj_TauSq', directory, 'i1./i2') % calculate map of mixed effects tau squared for the adjusted analysis
    
    ME_Adj_TauSqFile = [directory,'\ME_Adj_TauSq.nii,1']; % sets file path for the map of tau squared (adjusted)
    
    %%
    
    group0count = 1; % starts index for group 0 if 2 group mixed effects analysis
    group1count = 1; % starts index for group 1 if 2 group mixed effects analysis
    
    for k = 1:numstud
        
        WithinVarName = names(k).name; % sets file name for the within-study variance map
        WithinVarFile = [directory,'\',WithinVarName,',1']; % sets file path for the within-study variance map
        EstimatesName = [WithinVarName(1:end-8),'.nii']; % sets file name for the study effect size map
        EstimatesFile = [directory,'\',EstimatesName,',1']; % sets file path for the study effect size map
        Name = WithinVarName(4:end-8); % sets the name of the study
        position = find(~cellfun(@isempty,(strfind(cellstr(sdm_tab_cell{1,1}), Name)))); % finds the position of the study in the sdm_table
        group = sdm_tab_cell{colnum,1}(position); % finds which group the study is in
        BinarisedName = ['Binarised_',EstimatesName]; % sets  file name for the binarised map of where there is data (created below)
        BinarisedFile = [directory,'\',BinarisedName,',1']; % sets file path for the binarised map of where there is data (created below)
        ME_All_Wt_x_EstName = ['ME_All_Wt_x_Est_',Name]; % sets file name for the mixed effects weighting multipled by the effect-size estimate (original)
        ME_All_WtName = ['ME_All_Wt_',Name]; % sets file name for the mixed effects weighting (original)
        ME_All_WtFile = [directory,'\',ME_All_WtName,'.nii,1']; % sets file path for the mixed effects weighting (original)
        ME_Adj_Wt_x_EstName = ['ME_Adj_Wt_x_Est_',Name]; % sets file name for the mixed effects weighting multipled by the effect-size estimate (adjusted)
        ME_Adj_WtName = ['ME_Adj_Wt_',Name]; % sets file name for the mixed effects weighting (adjusted)
        ME_Adj_WtFile = [directory,'\',ME_Adj_WtName,'.nii,1']; % sets file path for the mixed effects weighting (adjusted)
        
        InputFiles = {
            ME_All_TauSqFile
            WithinVarFile
            };
        run_imcalc(InputFiles, ME_All_WtName, directory, '1./(i1+i2)') % calculates total inverse variance (ME weighting, original analysis) for each study by 1 / (between + within study variances)
        
        InputFiles = {
            ME_Adj_TauSqFile
            WithinVarFile
            BinarisedFile
            };
        run_imcalc(InputFiles, ME_Adj_WtName, directory, '(1./(i1+i2)).*i3') % calculates total inverse variance (ME weighting, adjusted analysis) for each study by 1 / (between + within study variances)
        
        InputFiles = {
            EstimatesFile
            ME_All_WtFile
            };
        run_imcalc(InputFiles, ME_All_Wt_x_EstName, directory, 'i1.*i2') % multiples the effect size estimates by the ME weighting for the original analysis
        
        InputFiles = {
            EstimatesFile
            ME_Adj_WtFile
            };
        run_imcalc(InputFiles, ME_Adj_Wt_x_EstName, directory, 'i1.*i2') % multiples the effect size estimates by the ME weighting for the adjusted analysis
        
        ME_All_Wt_x_EstFile = [directory,'\',ME_All_Wt_x_EstName,'.nii,1']; % sets file paths
        ME_All_WtFile = [directory,'\',ME_All_WtName,'.nii,1'];
        ME_Adj_Wt_x_EstFile = [directory,'\',ME_Adj_Wt_x_EstName,'.nii,1'];
        ME_Adj_WtFile = [directory,'\',ME_Adj_WtName,'.nii,1'];
        
        if group == 0 % creates lists of the relevant files for group 0
            ME_All_Wt_x_EstListGroup0{group0count,1}=ME_All_Wt_x_EstFile;
            ME_All_WtListGroup0{group0count,1}=ME_All_WtFile;
            ME_Adj_Wt_x_EstListGroup0{group0count,1}=ME_Adj_Wt_x_EstFile;
            ME_Adj_WtListGroup0{group0count,1}=ME_Adj_WtFile;
            group0count = group0count + 1;
        elseif group == 1 % creates lists of the relevant files for group 1
            ME_All_Wt_x_EstListGroup1{group1count,1}=ME_All_Wt_x_EstFile;
            ME_All_WtListGroup1{group1count,1}=ME_All_WtFile;
            ME_Adj_Wt_x_EstListGroup1{group1count,1}=ME_Adj_Wt_x_EstFile;
            ME_Adj_WtListGroup1{group1count,1}=ME_Adj_WtFile;
            group1count = group1count + 1;
        end
        
    end
    
    %%
    
    InputFiles = ME_All_Wt_x_EstListGroup0;
    run_imcalc(InputFiles, 'Sum_ME_All_Estimates_Group0', directory, sumcalcGroup0) % sums the effect size estimates (original analysis) for group 0
    
    InputFiles = ME_All_Wt_x_EstListGroup1;
    run_imcalc(InputFiles, 'Sum_ME_All_Estimates_Group1', directory, sumcalcGroup1) % sums the effect size estimates (original analysis) for group 1
    
    InputFiles = ME_Adj_Wt_x_EstListGroup0;
    run_imcalc(InputFiles, 'Sum_ME_Adj_Estimates_Group0', directory, sumcalcGroup0) % sums the effect size estimates (adjusted analysis) for group 0
    
    InputFiles = ME_Adj_Wt_x_EstListGroup1;
    run_imcalc(InputFiles, 'Sum_ME_Adj_Estimates_Group1', directory, sumcalcGroup1) % sums the effect size estimates (adjusted analysis) for group 1
    
    InputFiles = ME_All_WtListGroup0;
    run_imcalc(InputFiles, 'Sum_ME_All_Wt_Group0', directory, sumcalcGroup0); % sums the weightings (original analysis) for group 0
    
    InputFiles = ME_All_WtListGroup1;
    run_imcalc(InputFiles, 'Sum_ME_All_Wt_Group1', directory, sumcalcGroup1); % sums the weightings (original analysis) for group 1
    
    InputFiles = ME_Adj_WtListGroup0;
    run_imcalc(InputFiles, 'Sum_ME_Adj_Wt_Group0', directory, sumcalcGroup0); % sums the weightings (adjusted analysis) for group 0
    
    InputFiles = ME_Adj_WtListGroup1;
    run_imcalc(InputFiles, 'Sum_ME_Adj_Wt_Group1', directory, sumcalcGroup1); % sums the weightings (adjusted analysis) for group 1
    
    Sum_ME_All_EstimatesFileGroup0 = [directory,'\','Sum_ME_All_Estimates_Group0.nii,1']; % sets file paths
    Sum_ME_All_WtFileGroup0 = [directory,'\','Sum_ME_All_Wt_Group0.nii,1'];
    Sum_ME_Adj_EstimatesFileGroup0 = [directory,'\','Sum_ME_Adj_Estimates_Group0.nii,1'];
    Sum_ME_Adj_WtFileGroup0 = [directory,'\','Sum_ME_Adj_Wt_Group0.nii,1'];
    
    Sum_ME_All_EstimatesFileGroup1 = [directory,'\','Sum_ME_All_Estimates_Group1.nii,1'];
    Sum_ME_All_WtFileGroup1 = [directory,'\','Sum_ME_All_Wt_Group1.nii,1'];
    Sum_ME_Adj_EstimatesFileGroup1 = [directory,'\','Sum_ME_Adj_Estimates_Group1.nii,1'];
    Sum_ME_Adj_WtFileGroup1 = [directory,'\','Sum_ME_Adj_Wt_Group1.nii,1'];
    
    InputFiles = {
        Sum_ME_All_EstimatesFileGroup0
        Sum_ME_All_WtFileGroup0
        };
    run_imcalc(InputFiles, 'ME_All_EffectSize_Group0', directory, 'i1./i2') % divides sum of the estimates by the sum of the weightings to give the mixed effects effect size image, original analysis, group 0
    
    InputFiles = {
        Sum_ME_All_EstimatesFileGroup1
        Sum_ME_All_WtFileGroup1
        };
    run_imcalc(InputFiles, 'ME_All_EffectSize_Group1', directory, 'i1./i2') % divides sum of the estimates by the sum of the weightings to give the mixed effects effect size image, original analysis, group 1
    
    InputFiles = {
        Sum_ME_Adj_EstimatesFileGroup0
        Sum_ME_Adj_WtFileGroup0
        };
    run_imcalc(InputFiles, 'ME_Adj_EffectSize_Group0', directory, 'i1./i2') % divides sum of the estimates by the sum of the weightings to give the mixed effects effect size image, adjusted analysis, group 0
    
    InputFiles = {
        Sum_ME_Adj_EstimatesFileGroup1
        Sum_ME_Adj_WtFileGroup1
        };
    run_imcalc(InputFiles, 'ME_Adj_EffectSize_Group1', directory, 'i1./i2') % divides sum of the estimates by the sum of the weightings to give the mixed effects effect size image, adjusted analysis, group 1
        
    ME_All_EffectSizeFileGroup0 = [directory,'\','ME_All_EffectSize_Group0.nii,1']; % sets file paths
    ME_All_EffectSizeFileGroup1 = [directory,'\','ME_All_EffectSize_Group1.nii,1'];
    ME_Adj_EffectSizeFileGroup0 = [directory,'\','ME_Adj_EffectSize_Group0.nii,1'];
    ME_Adj_EffectSizeFileGroup1 = [directory,'\','ME_Adj_EffectSize_Group1.nii,1'];
    ME_All_EffectSizeFileG1_G0 = [directory,'\','ME_All_EffectSize_G1-G0.nii,1'];
    ME_Adj_EffectSizeFileG1_G0 = [directory,'\','ME_Adj_EffectSize_G1-G0.nii,1'];
    
    InputFiles = {Sum_ME_All_WtFileGroup0};
    run_imcalc(InputFiles, 'ME_All_Variance_Group0', directory, '1./i1') % calculates the variance for the original analysis, group 0 (1 / sum of weightings)
    
    InputFiles = {Sum_ME_All_WtFileGroup1};
    run_imcalc(InputFiles, 'ME_All_Variance_Group1', directory, '1./i1') % calculates the variance for the original analysis, group 1 (1 / sum of weightings)
    
    InputFiles = {Sum_ME_Adj_WtFileGroup0};
    run_imcalc(InputFiles, 'ME_Adj_Variance_Group0', directory, '1./i1') % calculates the variance for the adjusted analysis, group 0 (1 / sum of weightings)
    
    InputFiles = {Sum_ME_Adj_WtFileGroup1};
    run_imcalc(InputFiles, 'ME_Adj_Variance_Group1', directory, '1./i1') % calculates the variance for the adjusted analysis, group 1 (1 / sum of weightings)
    
    ME_All_VarianceFileGroup0 = [directory,'\','ME_All_Variance_Group0.nii,1']; % sets file paths
    ME_Adj_VarianceFileGroup0 = [directory,'\','ME_Adj_Variance_Group0.nii,1'];
    ME_All_VarianceFileGroup1 = [directory,'\','ME_All_Variance_Group1.nii,1'];
    ME_Adj_VarianceFileGroup1 = [directory,'\','ME_Adj_Variance_Group1.nii,1'];
    ME_All_VarianceFileG1_G0 = [directory,'\','ME_All_Variance_G1-G0.nii,1'];
    ME_Adj_VarianceFileG1_G0 = [directory,'\','ME_Adj_Variance_G1-G0.nii,1'];
    
    InputFiles = {
        ME_All_EffectSizeFileGroup1
        ME_All_EffectSizeFileGroup0};
    run_imcalc(InputFiles, 'ME_All_EffectSize_G1-G0', directory, 'i1-i2') % calculates the effect size map for the difference in effect sizes for group 1 > group 0, original analysis
    
    InputFiles = {
        ME_All_VarianceFileGroup1
        ME_All_VarianceFileGroup0};
    run_imcalc(InputFiles, 'ME_All_Variance_G1-G0', directory, 'i1+i2') % calculates the combined variance map for group 1 and group 0, original analysis (used in the analysis of the difference between groups)
    
    InputFiles = {
        ME_Adj_EffectSizeFileGroup1
        ME_Adj_EffectSizeFileGroup0};
    run_imcalc(InputFiles, 'ME_Adj_EffectSize_G1-G0', directory, 'i1-i2') % calculates the effect size map for the difference in effect sizes for group 1 > group 0, adjusted analysis
    
    InputFiles = {
        ME_Adj_VarianceFileGroup1
        ME_Adj_VarianceFileGroup0};
    run_imcalc(InputFiles, 'ME_Adj_Variance_G1-G0', directory, 'i1+i2') % calculates the combined variance map for group 1 and group 0, adjusted analysis (used in the analysis of the difference between groups)
    
    ME_var_neg_vals = 1;
    for x = 1:6 % checks whether any voxels in the variance maps have negative values which may happen in the adjusted analysis
        tocheck = {'ME_All_Variance_Group1.nii', 'ME_All_Variance_Group0.nii', 'ME_Adj_Variance_Group1.nii', 'ME_Adj_Variance_Group0.nii', 'ME_All_Variance_G1-G0.nii', 'ME_Adj_Variance_G1-G0.nii'};
        varfile = tocheck{x};
        readin = spm_vol(varfile);
        getfile = spm_read_vols(readin);
        minimum = min(min(min(getfile)));
        if minimum < 0 % any negative values need to be removed before conversion to standard error (square root of variance)
            copyname = [varfile(1:end-4),'_NegValues.nii'];
            copyfile(varfile, copyname);
            NegValuesFile = [directory,'\',copyname,',1']; % creates a copy of the original variance map for the user to investigate which voxels the issues were in
            InputFiles = {NegValuesFile};
            run_imcalc(InputFiles, varfile(1:end-4), directory, 'i1.*(i1>0)')
            ME_var_neg_vals_log{ME_var_neg_vals,1} = varfile; % keeps a log of any files which have been edited to report at the end
            ME_var_neg_vals = ME_var_neg_vals + 1;
        else
        end
        
    end
    
    InputFiles = {ME_All_VarianceFileGroup0};
    run_imcalc(InputFiles, 'ME_All_StandardError_Group0', directory, 'sqrt(i1)') % converts variance to standard error for the original analysis, group 0
    
    InputFiles = {ME_All_VarianceFileGroup1};
    run_imcalc(InputFiles, 'ME_All_StandardError_Group1', directory, 'sqrt(i1)') % converts variance to standard error for the original analysis, group 1
    
    InputFiles = {ME_Adj_VarianceFileGroup0};
    run_imcalc(InputFiles, 'ME_Adj_StandardError_Group0', directory, 'sqrt(i1)') % converts variance to standard error for the adjusted analysis, group 0
    
    InputFiles = {ME_Adj_VarianceFileGroup1};
    run_imcalc(InputFiles, 'ME_Adj_StandardError_Group1', directory, 'sqrt(i1)') % converts variance to standard error for the adjusted analysis, group 1
    
    InputFiles = {ME_All_VarianceFileG1_G0};
    run_imcalc(InputFiles, 'ME_All_StandardError_G1-G0', directory, 'sqrt(i1)') % converts variance to standard error for the original analysis, group 1 > group 0
    
    InputFiles = {ME_Adj_VarianceFileG1_G0};
    run_imcalc(InputFiles, 'ME_Adj_StandardError_G1-G0', directory, 'sqrt(i1)') % converts variance to standard error for the adjusted analysis, group 1 > group 0
    
    ME_All_StandardErrorFileGroup0 = [directory,'\','ME_All_StandardError_Group0.nii,1']; % sets file paths
    ME_Adj_StandardErrorFileGroup0 = [directory,'\','ME_Adj_StandardError_Group0.nii,1'];
    ME_All_StandardErrorFileGroup1 = [directory,'\','ME_All_StandardError_Group1.nii,1'];
    ME_Adj_StandardErrorFileGroup1 = [directory,'\','ME_Adj_StandardError_Group1.nii,1'];
    ME_All_StandardErrorFileG1_G0 = [directory,'\','ME_All_StandardError_G1-G0.nii,1'];
    ME_Adj_StandardErrorFileG1_G0 = [directory,'\','ME_Adj_StandardError_G1-G0.nii,1'];
    
    InputFiles = {
        ME_All_EffectSizeFileGroup0
        ME_All_StandardErrorFileGroup0
        };
    run_imcalc(InputFiles, 'ME_All_Z_Group0', directory, 'i1./i2') % calculates the map of Z scores for the mixed effects for group 0, original analysis
    
    InputFiles = {
        ME_All_EffectSizeFileGroup1
        ME_All_StandardErrorFileGroup1
        };
    run_imcalc(InputFiles, 'ME_All_Z_Group1', directory, 'i1./i2') % calculates the map of Z scores for the mixed effects for group 1, original analysis
    
    InputFiles = {
        ME_Adj_EffectSizeFileGroup0
        ME_Adj_StandardErrorFileGroup0
        };
    run_imcalc(InputFiles, 'ME_Adj_Z_Group0', directory, 'i1./i2') % calculates the map of Z scores for the mixed effects for group 0, adjusted analysis
    
    InputFiles = {
        ME_Adj_EffectSizeFileGroup1
        ME_Adj_StandardErrorFileGroup1
        };
    run_imcalc(InputFiles, 'ME_Adj_Z_Group1', directory, 'i1./i2') % calculates the map of Z scores for the mixed effects for group 1, adjusted analysis
    
    InputFiles = {
        ME_All_EffectSizeFileG1_G0
        ME_All_StandardErrorFileG1_G0
        };
    run_imcalc(InputFiles, 'ME_All_Z_G1-G0', directory, 'i1./i2') % calculates the map of Z scores for the mixed effects for group 1 > group 0, original analysis
    
    InputFiles = {
        ME_Adj_EffectSizeFileG1_G0
        ME_Adj_StandardErrorFileG1_G0
        };
    run_imcalc(InputFiles, 'ME_Adj_Z_G1-G0', directory, 'i1./i2') % calculates the map of Z scores for the mixed effects for group 1 > group 0, adjusted analysis
        
    disp('Mixed effects analysis for means and comparisons between 2 groups complete.')
    
    if exist('ME_var_neg_vals_log') == 1 % if any variance files were changed due to negative values, reports them here
        disp('IMPORTANT - the following variance files contained negative values so were masked to exclude these before standard error could be calculated. See files *_Variance_*_NegValues.nii for the location of these voxels.')
        celldisp(ME_var_neg_vals_log)
    else
    end
end

