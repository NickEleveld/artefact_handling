%% Visualise data: ABP-TCD, NIRS, Accelerometer
%% For ABP + TCD + NIRS measurements
inputfolder_alldata = dir('MASKED');
outputfolder = 'MASKED';
%% For NIRS-only measurements (also contains measurements with bad ABP-signal)
inputfolder_alldata = dir('MASKED');
outputfolder = 'MASKED';

inputfolder_alldata = inputfolder_alldata([inputfolder_alldata.isdir]==0);
[~,~,ac] = fileparts({inputfolder_alldata.name});
inputfolder_alldata = inputfolder_alldata(strcmp(ac,'.mat'));
%%
for i=1:length(inputfolder_alldata)
syncedData = load(fullfile(inputfolder_alldata(i).folder,inputfolder_alldata(i).name));
disp(inputfolder_alldata(i).name)
fieldnames_syncedData    = fieldnames(syncedData); %to automatically unstructure the dataset
syncedData = syncedData.(fieldnames_syncedData{1});

close all
figure_name = tiledlayout(6,1,'TileSpacing','tight','Padding','tight');
ax(1) = nexttile;
plotABP =plot(syncedData.Timestamp,syncedData.ABP,'DisplayName','syncedData_combined.ABP');
ylim([0,250])
 legend
set(gca,'XTick',[])
ax(2) = nexttile;
plotTCD1 = plot(syncedData.Timestamp,syncedData.(2),'DisplayName',syncedData.Properties.VariableNames{2});hold on;
plotTCD2 = plot(syncedData.Timestamp,syncedData.(3),'DisplayName',syncedData.Properties.VariableNames{3});hold off;
ylim([0,250])
 legend
set(gca,'XTick',[])
ax(3) = nexttile;
plotHHb_R   = plot(syncedData.Timestamp,syncedData.Rx1_Tx3HHb-mean(syncedData.Rx1_Tx1HHb),'DisplayName','syncedData_combined.Rx1_Tx1HHb');hold on;
plotOxyHb_R = plot(syncedData.Timestamp,syncedData.Rx1_Tx3O2Hb-mean(syncedData.Rx1_Tx1O2Hb),'DisplayName','syncedData_combined.Rx1_Tx2O2Hb');
 legend
set(gca,'XTick',[])
ax(4) = nexttile;
plotHHb_R   = plot(syncedData.Timestamp,syncedData.Rx1_Tx3HHb-mean(syncedData.Rx1_Tx2HHb),'DisplayName','syncedData_combined.Rx1_Tx2HHb');hold on;
% % plotOxyHb_R = plot(syncedData.Timestamp,syncedData.Rx1_Tx3O2Hb-mean(syncedData.Rx1_Tx2O2Hb),'DisplayName','syncedData_combined.Rx1_Tx2O2Hb');

% plotHHb_L   = plot(syncedData.Timestamp,syncedData.Rx8_Tx8HHb-mean(syncedData.Rx8_Tx8HHb),'DisplayName','syncedData_combined.Rx8_Tx8HHb');hold on;
% plotOxyHb_L = plot(syncedData.Timestamp,syncedData.Rx8_Tx8O2Hb-mean(syncedData.Rx8_Tx8O2Hb),'DisplayName','syncedData_combined.Rx8_Tx8O2Hb');hold off;
 legend
set(gca,'XTick',[])
ax(5) = nexttile;
plotHHb_R   = plot(syncedData.Timestamp,syncedData.Rx1_Tx3HHb-mean(syncedData.Rx1_Tx3HHb),'DisplayName','syncedData_combined.Rx1_Tx3HHb');hold on;
plotOxyHb_R = plot(syncedData.Timestamp,syncedData.Rx1_Tx3O2Hb-mean(syncedData.Rx1_Tx3O2Hb),'DisplayName','syncedData_combined.Rx1_Tx3O2Hb');hold off;

syncedData_ACC = sqrt(syncedData.ACC_X.^2+syncedData.ACC_Y.^2+syncedData.ACC_Z.^2);
% plotACC = plot(syncedData.Timestamp,syncedData_ACC-mean(syncedData_ACC)    ,'DisplayName','syncedData_ACC');hold on;
% ylim([-10*std(syncedData_ACC-mean(syncedData_ACC)),10*std(syncedData_ACC-mean(syncedData_ACC))]);
 legend
set(gca,'XTick',[])
if any(contains(syncedData.Properties.VariableNames,'tau final'))
ax(6) = nexttile;
%% For NIRS-only
% plotHHb_R   = plot(syncedData.Timestamp,syncedData.Rx1_Tx3HHb-mean(syncedData.Rx1_Tx3HHb),'DisplayName','syncedData_combined.Rx1_Tx3HHb');hold on;
%% For TCD+ABP+NIRS
plottau = plot(syncedData.Timestamp,syncedData.("tau final")    ,'DisplayName','tau final (samples)');hold on;
yyaxis right
plotmax_xc = plot(syncedData.Timestamp,syncedData.max_xc_final, 'DisplayName','Cross-corr'); hold off
%%
legend
end
java.lang.System.gc()
linkaxes(ax,'x')
% legend
% figure_name.WindowState = 'maximized';
set(gcf, 'Position', [1,49,1920,955]);%get(gcf,'Position'));
brush on
fontsize(gcf,15,'points')
%% Marking per variable
% disp('select the artefacts and press ENTER')
pause
idABP     = boolean(get(plotABP, 'BrushData'))';     if isempty(idABP); idABP = false(height(syncedData),1); end 
idTCD1    = boolean(get(plotTCD1, 'BrushData'))';    if isempty(idTCD1); idTCD1 = false(height(syncedData),1); end 
idTCD2    = boolean(get(plotTCD2, 'BrushData'))';    if isempty(idTCD2); idTCD2 = false(height(syncedData),1); end 
idOxyHb_R = boolean(get(plotOxyHb_R, 'BrushData'))'; if isempty(idOxyHb_R); idOxyHb_R = false(height(syncedData),1); end 
idHHb_R   = boolean(get(plotHHb_R, 'BrushData'))';   if isempty(idHHb_R); idHHb_R = false(height(syncedData),1); end 
idOxyHb_L = boolean(get(plotOxyHb_L, 'BrushData'))'; if isempty(idOxyHb_L); idOxyHb_L = false(height(syncedData),1); end 
idHHb_L   = boolean(get(plotHHb_L, 'BrushData'))';   if isempty(idHHb_L); idHHb_L = false(height(syncedData),1); end 
idACC     = boolean(get(plotACC, 'BrushData'))';     if isempty(idACC); idACC= false(height(syncedData),1); end 
if any(contains(syncedData.Properties.VariableNames,'tau final'))
    idtau     = boolean(get(plottau, 'BrushData'))';     if isempty(idtau); idtau = false(height(syncedData),1); end 
    idmax_xc  = boolean(get(plotmax_xc, 'BrushData'))';  if isempty(idmax_xc); idmax_xc = false(height(syncedData),1); end 
else 
    idtau = false(height(syncedData),1);
    idmax_xc = false(height(syncedData),1);
end
idBrush_combined = any([idABP;idTCD1;idTCD2;idOxyHb_R;idHHb_R;idOxyHb_L;idHHb_L;idACC;idtau;idmax_xc],1);

variablenames_artefacts = {'idABP','idTCD1','idTCD2','idOxyHb_R','idHHb_R','idOxyHb_L','idHHb_L','idACC','idtau','idmax_xc'};
artefacts_table = timetable(syncedData.Properties.RowTimes,idABP,idTCD1,idTCD2,idOxyHb_R,idHHb_R,idOxyHb_L,idHHb_L,idACC,idtau,idmax_xc,'VariableNames',variablenames_artefacts);

[~,outputname,~] = fileparts(inputfolder_alldata(i).name);
% % save(fullfile(outputfolder,[outputname,'_artefacts.mat']),'artefacts_table');
% % saveas(gcf,fullfile(outputfolder,['\figs (helaas zonder brush waardes)\',outputname,'_artefacts.fig']))
end