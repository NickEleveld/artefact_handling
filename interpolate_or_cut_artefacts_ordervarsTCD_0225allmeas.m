%% 1) Load data

%% For all measurements (TCD+ABP+NIRS & NIRS-only)
inputfolder = dir('MASKED');
inputfolder = inputfolder([inputfolder.isdir]==0);
outputfolder = ''MASKED';
outputfolder_combined = 'MASKED';

% 2a) Load artefactdata
inputfolder_artefacts = dir('MASKED');
inputfolder_artefacts = inputfolder_artefacts([inputfolder_artefacts.isdir]==0);

%% Settings
sampleFreq = 200; % in Hz
threshold_short_artefact = 3; % in sec
threshold_filelength = 101;%151; % in sec
% 2b) Make a conversion table between the annotated artefacts and the
% variables in the timetables
conversion_table = {
'idABP'    , 'ABP';
'idTCD1'   , 2;
'idTCD2'   , 3;
'idOxyHb_R', 'Rx1_Tx3O2Hb';
'idHHb_R'  , 'Rx1_Tx3HHb';
'idOxyHb_L', 'Rx8_Tx8O2Hb';
'idHHb_L'  , 'Rx8_Tx8HHb';
'idACC'    , 'ACC_X';
'idtau'    , 'tau final';
'idmax_xc' , 'max_xc_final'};
%%
for i=1:length(inputfolder)
%% Load data and corresponding artefact files
    data_segment_original     = load(fullfile(inputfolder(i).folder,inputfolder(i).name));
    fieldnames_data           = fieldnames(data_segment_original); %to automatically unstructure the dataset
    data_segment_original     = data_segment_original.(fieldnames_data{1});  
    index_artefacts  = contains({inputfolder_artefacts.name},inputfolder(i).name(1:end-4));
    artefact_segment = load(fullfile(inputfolder_artefacts(i).folder,inputfolder_artefacts(index_artefacts).name));
    fieldnames_artefact = fieldnames(artefact_segment); %to automatically unstructure the dataset
    artefact_segment    = artefact_segment.(fieldnames_artefact{1});
    
    % Linearly interpolate short artefacts
    for i_var_art = 1:width(artefact_segment)
        var_art = artefact_segment.Properties.VariableNames{i_var_art};
        artefact_selectedTEMP = artefact_segment(:,var_art);
        if artefact_selectedTEMP{1,:} ==1
            find_artefacts_onevarstart = 1;
        else
            find_artefacts_onevarstart = [];
        end
        if artefact_selectedTEMP{end,:} ==1
            find_artefacts_onevarend = height(artefact_selectedTEMP);
        else
            find_artefacts_onevarend = [];
        end
    
        find_artefacts_onevar = find([diff(artefact_selectedTEMP{:,:}) ~= 0]);  % TRUE if values change
        length(find_artefacts_onevar)
        find_artefacts_onevar_all = [find_artefacts_onevarstart; find_artefacts_onevar; find_artefacts_onevarend];
        find_artefacts_onevar = find_artefacts_onevar_all;

        difference_boolean_onevar_tworows = [find_artefacts_onevar(1:2:end),find_artefacts_onevar(2:2:end)];

        for i_artefacts_interp =1:size(difference_boolean_onevar_tworows,1)

        x = difference_boolean_onevar_tworows(i_artefacts_interp,:);
% % % % % % % %         VARIABLE = var_art -->conversion_table(find,1) --> conversion_table(find,2)
        variablename_data = conversion_table{find(strcmp(conversion_table(:,1)',var_art)),2};
        v = data_segment_original{x,variablename_data};
        vq= [x(1):1:x(2)];
        v_new = interp1(x,v,vq,'linear');
        data_segment_original{vq,variablename_data} = v_new';
        end
    end

    % 3) Combine artefacts_table into artefacts_any (use later)
    idBrush_combined = any(artefact_segment{:,:},2);
    
    % 4) Find artefact segments in individual columns < 3 sec
    difference_boolean = [true; diff(idBrush_combined) ~= 0; true];  % TRUE if values change
    n_length_segments = diff(find(difference_boolean));     % Number of repetitions
    Y = repelem(n_length_segments, n_length_segments); if isrow(Y); Y=Y';end
    Y2 = [Y,idBrush_combined];
    long_artefacts = (Y2(:,1)>= threshold_short_artefact*sampleFreq & Y2(:,2)==1);
    n = diff((long_artefacts));               % Number of repetitions
    found_changes_segment = find(n); 
    suffix = '_';
    if isempty(found_changes_segment);
        suffix = 'complete';
    elseif length(found_changes_segment)==1
        if long_artefacts(found_changes_segment-1) ==1
           found_changes_segment = [1; found_changes_segment];
        elseif long_artefacts(found_changes_segment-1) ==0
           found_changes_segment = [found_changes_segment; length(long_artefacts)-1];
        end
    end 
%     found_changes_segment_v2 = [1;found_changes_segment;height(data_segment)];
    iii=1;
%     for ii=1:2:length(found_changes_segment_v2)
    for ii=1:2:length(found_changes_segment)
        data_segment_original{found_changes_segment(ii):found_changes_segment(ii+1),:} = 9999; % NEW
        iii=iii+1;
    end
    find(strcmp(data_segment_original.Properties.VariableNames,'MCA_L_env') | strcmp(data_segment_original.Properties.VariableNames,'MCAL'))

%% Order TCD variables: TCD RIGHT: COLUMN 2, TCD LEFT: COLUMN 3
column_nr_left = find(strcmp(data_segment_original.Properties.VariableNames,'MCA_L_env') | strcmp(data_segment_original.Properties.VariableNames,'MCAL'));
column_nr_right = find(strcmp(data_segment_original.Properties.VariableNames,'MCA_R_env') | strcmp(data_segment_original.Properties.VariableNames,'MCAR'));

if isempty(column_nr_left) | isempty(column_nr_right)
    warning('No columns found for MCA_L_env or MCAL OR MCA_R or MCAR: Making placeholder CBFv variables')
    data_segment = data_segment_original;
    data_segment.CBFv_R = zeros(height(data_segment),1); data_segment = movevars(data_segment, "CBFv_R", "After", "ABP");
    data_segment.CBFv_L = zeros(height(data_segment),1); data_segment = movevars(data_segment, "CBFv_L", "After", "CBFv_R");
else
    data_segment = data_segment_original;
    data_segment(:,2) = data_segment_original(:,column_nr_right);
    data_segment.Properties.VariableNames(2) = "CBFv_R";data_segment_original.Properties.VariableNames(column_nr_right);
    data_segment(:,3) = data_segment_original(:,column_nr_left);
    data_segment.Properties.VariableNames(3) = "CBFv_L";data_segment_original.Properties.VariableNames(column_nr_left);
end

%% If no columns tau final and max_xc_final are present, add them manually
no_taucolumns = ~any(strcmp(data_segment_original.Properties.VariableNames,'tau final') | strcmp(data_segment_original.Properties.VariableNames,'max_xc_final'));
if no_taucolumns
    data_segment{:,"tau final"} = NaN(height(data_segment),1);
    data_segment{:,"max_xc_final"} = NaN(height(data_segment),1);  
    data_segment = movevars(data_segment,["tau final","max_xc_final"],"After","CBFv_L");
end
%%
% writetimetable(data_segment,fullfile(outputfolder,[inputfolder(i).name(1:end-4),'_segment',num2str(iii),suffix,'.txt']));% NOT SAVED TO SAVE SPACE
save(fullfile(outputfolder,[inputfolder(i).name(1:end-4),'_segment',num2str(iii),suffix,'.mat']),"data_segment"); % NEW

end
%% Combine files and remove excess columns to make the files smaller
excess_columns = ["ACC_Y","ACC_Z","GYR_X","GYR_Y","GYR_Z","Rx2_Tx5O2Hb","Rx2_Tx5HHb","Rx3_Tx6O2Hb","Rx3_Tx6HHb"]';

outputfolder_content = dir(outputfolder);
outputfolder_content = outputfolder_content(~(strcmp({outputfolder_content.name},'.') |strcmp({outputfolder_content.name},'..')));

dates_only = extractBefore({outputfolder_content.name},'intraclamp');
dates_only2 = extractBefore({outputfolder_content.name},'preclamp');
location_dates2 = ~cellfun('isempty',dates_only2);
dates_only(location_dates2) = dates_only2(location_dates2);
measurements = unique(dates_only');
timing = {'preclamp';'intraclamp'};
counter = 1;
files_with_artefactinlast150s = table();

for ii = 1:length(measurements)
    for iii=1:length(timing)
        loc_onemeasurement_onetime = contains({outputfolder_content.name},[measurements{ii},timing{iii}]);
        onemeasurement_onetime     = {outputfolder_content(loc_onemeasurement_onetime).name}';
        txt_onemeasurement_onetime = onemeasurement_onetime(contains(onemeasurement_onetime,'.txt'));
        mat_onemeasurement_onetime = onemeasurement_onetime(contains(onemeasurement_onetime,'.mat'));

        % Mat-files 
        for iiii= 1:length(mat_onemeasurement_onetime)
            loc_measurement_final = strcmp({outputfolder_content.name},mat_onemeasurement_onetime(iiii));
            data_segment_tobecombined = load(fullfile(outputfolder_content(loc_measurement_final).folder,outputfolder_content(loc_measurement_final).name));
            fieldnames_combined       = fieldnames(data_segment_tobecombined); %to automatically unstructure the dataset
            data_segment_tobecombined = data_segment_tobecombined.(fieldnames_combined{1});  

            % Remove excess columns
            excess_columns_all = [];
            for i_exce_colums = 1:length(excess_columns)
                excess_column_one = find(strcmp(data_segment_tobecombined.Properties.VariableNames,excess_columns(i_exce_colums)));
            excess_columns_all = [excess_columns_all,excess_column_one]; 
            end
            columns_tobe_maintained = true(width(data_segment_tobecombined),1);
            columns_tobe_maintained(excess_columns_all) = false;
            data_segment_tobecombined = data_segment_tobecombined(:,columns_tobe_maintained);            
             
            if iiii==1
                data_segment_combined = data_segment_tobecombined;
            else
                data_segment_tobecombined{1,:} = 8888;
                data_segment_combined = [data_segment_combined;data_segment_tobecombined];
            end
            % Remove excess columns
            excess_columns_all = [];
            for i_exce_colums = 1:length(excess_columns)
                excess_column_one = find(strcmp(data_segment_combined.Properties.VariableNames,excess_columns(i_exce_colums)));
            excess_columns_all = [excess_columns_all,excess_column_one]; 
            end
            columns_tobe_maintained = true(width(data_segment_combined),1);
            columns_tobe_maintained(excess_columns_all) = false;
            data_segment_combined = data_segment_combined(:,columns_tobe_maintained);

            save(fullfile(outputfolder_combined,[measurements{ii},timing{iii},'_combined.mat']),'data_segment_combined');
        end
    end
end