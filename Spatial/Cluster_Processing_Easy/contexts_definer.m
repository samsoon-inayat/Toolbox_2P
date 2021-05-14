%% for context definitions (only need to be run in the start)
% for ii = 1:size(sel_T,1)
ii = 4;
    this_table_row = sT(ii,:);
    recording_folder = get_recording_folder_from_table_row(this_table_row);
    pd_folder = get_processed_data_folder_from_table_row(this_table_row);
    edit_define_contexts_file(ei{ii});
    plotMarkers(ei{ii}.b,ei{ii}.b.photo_sensor_f,[],100,0);
    plotMarkers(ei{ii}.b,ei{ii}.b.air_puff_r,[],101,0);
% end