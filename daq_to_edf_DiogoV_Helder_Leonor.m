function [key_matrix_daq, clinical_information_daq] = daq_to_edf_DiogoV_Helder_Leonor(directory_Origin, directory_Destin)
%%   Converts every .daq file in the directory directory_Origin and converts them
%%   to multiple type files and stores them in directory_Destin_Edf directory and creates a Identification Table
%
%   Inputs:
%   directory_Origin - Directory where the .daq files are
%   directory_Destin - Directory where the .edf files will be stored
%
%   Outputs:
%   key_matrix_daq - Identification Table (Name - Id)
%   clinical_information_daq - Clinical information about the patients
%
%   Author: Diogo Ventura, Maria Leonor, Helder Lopes

    id_Start = 17;
    
    genders = readtable('gender_Table.xlsx'); %	gender_Table.xlsx must be in the same folder as the file function daq_to_edf_DiogoV_Helder_Leonor.m
    function_Directory = pwd;
    
    %   Check directory existence
    if ~isfolder(directory_Origin)
        error("Directory 'directory_Origin' does not exist.");
    end
    
    %   Checks if the directory_Destin_Edf directory exists, if not creates it
    if exist(directory_Destin, 'dir') ~= 7
        mkdir (directory_Destin)
    end
    
    %   Get files names
    file_Names_daq = dir(fullfile(directory_Origin, '*.daq'));
    file_Names_daq = {file_Names_daq.name};
    file_Names_mat = dir(fullfile(directory_Origin, '*.mat'));
    file_Names_mat = {file_Names_mat.name};
    
    %   Creates the key_matrix and clinical_inf_table table using the number of files
    number_files = length(file_Names_daq);
    key_matrix_daq = table('Size', [number_files, 2], 'VariableTypes', {'string', 'string'}, ...
        'VariableNames', {'Name', 'case_ID'});
    
    variableTypes = {'string', 'double', 'double', 'double', 'string', 'double', ...
        'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string'};
    variableNames = {'Filename', 'Heigth', 'Weight', 'Age', 'Gender', 'SamplingFrequency', ...
        'Vi', 'Ve', 'Vc', 'Va', 'Vm', 'Vh', 'Vf', 'Vx', 'Vy', 'Vz'};
    clinical_information_daq = table('Size', [number_files, numel(variableNames)], ...
        'VariableTypes', variableTypes, 'VariableNames', variableNames);
    
    frank_labels = {'Vi', 'Ve', 'Vc', 'Va', 'Vm', 'Vh', 'Vf', 'Vx', 'Vy', 'Vz'};
    % frank_labels_filt = cellfun(@(x) [x, '_filtered'], frank_labels, 'UniformOutput', false);
    
    birthdayFormats = {'dd:MM:yyyy', 'dd.MM.yyyy', 'dd-MM-yyyy', 'dd.MM.yy', 'dd-MM-yy', 'dd:MM:yy'};
    
    
    %%   Convertion cycling all files
    for i = 1:number_files
        try
            fileName_daq = file_Names_daq{i};
        catch
            error("No .daq files in " + directory_Origin + "!");
        end
        try
            fileName_mat = file_Names_mat{i};
        catch
            error("No .mat files in " + directory_Origin + "!");
        end
        [data, time, ~, ~, daqinfo] = daqread(fullfile(directory_Origin, fileName_daq));
        load(fullfile(directory_Origin, fileName_mat), 'P_C_DAQ_S');
    
        samplingFrequency = daqinfo.ObjInfo.SampleRate;
    
        % Signal trimming first and last 20 seconds
        data = data((samplingFrequency * 20) : (end - samplingFrequency * 20), 1:7);
        time = time(1:(end - 2 * samplingFrequency * 20) + 1);
    
        data = data - mean(data);
    
    
        % Frank Derivations
        info.labels = frank_labels;
        Va = data(:,4); Vc = data(:,3); Ve = data(:,2); Vf = data(:,7);
        Vh = data(:,6); Vi = data(:,1); Vm = data(:,5);
        Vx = 0.61 * Va + 0.171 * Vc - 0.781 * Vi;
        Vy = 0.655 * Vf + 0.345 * Vm - Vh;
        Vz = 0.133 * Va + 0.736 * Vm - 0.264 * Vi - 0.374 * Ve - 0.231 * Vc;
        Vx = Vx - mean(Vx);
        Vy = Vy - mean(Vy);
        Vz = Vz - mean(Vz);
        data = [Vi, Ve, Vc, Va, Vm, Vh, Vf, Vx, Vy, Vz]*(1e6);
    
        % Process subject name
        parts = {P_C_DAQ_S.subjectfirstname, P_C_DAQ_S.subjectlastname};
        if ~(isstrprop(P_C_DAQ_S.subjectfirstname(1), 'upper') && isstrprop(P_C_DAQ_S.subjectlastname(1), 'upper'))
            parts = cellfun(@(x) [upper(x(1)), lower(x(2:end))], parts, 'UniformOutput', false);
        end
        % Concatenate name parts
        name = strjoin(parts, ' ');
    
        % Process Birthday
        info.age = 0;
        subject_birthday = P_C_DAQ_S.subjectbirthday;
        if ~isempty(subject_birthday)
            datetimeArray = subject_birthday;
    
            parsedBirthday = [];
            for p = 1:length(birthdayFormats)
                try
                    parsedBirthday = datetime(datetimeArray, 'InputFormat', birthdayFormats{p}, 'PivotYear', 2000);
                    break;  % Break out of the loop if successfully parsed
                catch
                    continue;  % Try the next format if parsing fails
                end
            end
            if parsedBirthday.Year < 100
                parsedBirthday.Year = parsedBirthday.Year + 1900;
            end
            birth_Day = datestr(parsedBirthday, 'dd-mmm-yyyy');
            info.birthDay = birth_Day;
        end
    
    
        % Process Age and Heigth
        info.heigth = [];
        info.weigth = [];
        heigth_weigth = P_C_DAQ_S.subjectcomment;
        if ~isempty(heigth_weigth)
            if ~isnan(str2double(heigth_weigth(1)))
                heigth_weigth = strsplit(heigth_weigth, ' ');
                if heigth_weigth{1}(end) == ';'
                    heigth_weigth{1}(end) = [];
                end
                if heigth_weigth{1}(end)=='m'
                    heigth = heigth_weigth{1};
                    weigth = heigth_weigth{2};
                else
                    heigth = heigth_weigth{2};
                    weigth = heigth_weigth{1};
                end
                info.weigth = weigth;
                info.heigth = heigth;
            end
        end
        
        new_fileName = sprintf('case_%02d', i + id_Start - 1);

        %   Creating header for EDF
        gender_case = genders.Gender{strcmp(genders.Name, name)};
        hdr = create_hdr(data, daqinfo, new_fileName + " " + gender_case, info);
        if ~isempty(hdr.StartDate) && ~isempty(info.birthDay)
            info.age = floor(years(datetime(hdr.StartDate, 'InputFormat', 'dd.MM.yy') ...
                - datetime(info.birthDay, 'InputFormat', 'dd-MMM-yyyy')));
        end

        % Baseline Removal
        new_data = baseline_remover(time, data, 'sin5');  
        data_table = array2table([data, new_data], 'VariableNames', [frank_labels, frank_labels_filt]);

        %   Save .edf file
        new_fileName_EDF = new_fileName + ".edf";
        cd (directory_Destin)
        edfwrite(new_fileName_EDF, hdr, data, 'InputSampleType', "physical");
        cd (function_Directory)
        
        %   Save .csv file
        new_fileName_CSV = new_fileName + ".csv";
        fullPath_csv = fullfile(directory_Destin, new_fileName_CSV);
        writetable(data_table, fullPath_csv);
        
        %   Save .mat file
        new_fileName_MAT = new_fileName + ".mat";
        fullPath_mat = fullfile(directory_Destin, new_fileName_MAT);
        save(fullPath_mat, "data_table", "time");
        
        %   Save .txt file
        new_fileName_TXT = new_fileName + ".txt";
        fullPath_txt = fullfile(directory_Destin, new_fileName_TXT);
        save(fullPath_txt, 'data', '-ascii','-double');
        
        %   Storing name and id in Key Matrix table
        key_matrix_daq(i, :) = {{string(name)}, new_fileName};
        if ~isempty(info.heigth) && strcmp(info.heigth(end), 'm')
            heigth_num = info.heigth(1:end-1);
            heigth_num = strrep(heigth_num, ',', '.');
        else
            heigth_num = 0;
        end
        if ~isempty(info.weigth) && strcmp(info.weigth(end-1:end), 'kg')
            weigth_num = info.weigth(1:end-2);
        else
            weigth_num = 0;
        end

        clinical_information_daq(i, :) = [{new_fileName}, {str2double(heigth_num) * 100}, ...
            {str2double(weigth_num)}, {info.age}, {gender_case}, {samplingFrequency}, hdr.PhysicalDimensions];

        disp("File " + i + "/" + number_files + " done.");
    end
    
    %   Save tables to Excel in each directory
    excel_FilePath = fullfile(directory_Origin, 'key_matrix_daq.xlsx');
    writetable(key_matrix_daq, excel_FilePath);
    excel_FilePath = fullfile(directory_Destin, 'clinical_information_daq.xlsx');
    writetable(clinical_information_daq, excel_FilePath);
    
    %   Return to original directory
    cd (function_Directory)
end


function hdr = create_hdr(data, daqinfo, new_fileName, info)
%   Creates an header using data, sampleRate, daqinfo and new_fileName inputs
%   for all the necessary information

    % Initialize header
    hdr = edfheader("EDF+");
    
    % Set patient information
    header_Patient = string(new_fileName);
    fieldsToConcatenate = {'birthDay', 'weigth', 'heigth'};
    for i = 1:length(fieldsToConcatenate)
        field = fieldsToConcatenate{i};
        try
            header_Patient = strcat(header_Patient, " ", info.(field));
        catch
        end
    end
    header_Patient = regexprep(header_Patient, '\s*$', '');
    hdr.Patient = header_Patient;
    

    % Set number of data records
    hdr.NumDataRecords = 1;
    
    % Compute duration of each data record
    recordDuration = seconds(size(data, 1) / daqinfo.ObjInfo.SampleRate);
    hdr.DataRecordDuration = recordDuration;
    
    % Set number of signals
    hdr.NumSignals = size(data, 2);
    
    % Set signal labels
    hdr.SignalLabels = info.labels;
    
    % Compute min and max values for all signals
    minData = min(data);
    maxData = max(data);
    hdr.PhysicalMin = floor(minData * 1e3) / 1e3; % Rounding 5th decimal down
    hdr.PhysicalMax = ceil(maxData * 1e3) / 1e3; % Rounding 5th decimal up
    % hdr.PhysicalMin = minData; % Rounding 5th decimal down
    % hdr.PhysicalMax = maxData; % Rounding 5th decimal up
    
    % Set physical dimensions
    % if daqinfo.ObjInfo.Channel(1).Units == "Volts"
    %     units = 'V';
    % end
    micro = 'uV';
    hdr.PhysicalDimensions = repelem({micro}, hdr.NumSignals);
    
    % Set digital min and max values
    hdr.DigitalMin = repelem(2^15 * (-1), hdr.NumSignals);
    hdr.DigitalMax = repelem(2^15 - 1, hdr.NumSignals);
    
    % Set start date and time
    datetimeArray = daqinfo.ObjInfo.EventLog(1).Data.AbsTime;
    datetimeObject = datetime(datetimeArray(1:6));
    startDate = datestr(datetimeObject, 'dd.mm.yy'); %#ok<*DATST>
    startDate2 = datestr(datetimeObject, 'dd-mmm-yyyy');
    startTime = datestr(datetimeObject, 'HH.MM.SS');
    hdr.StartDate = startDate;
    hdr.StartTime = startTime;
    
    % Set recording information
    hdr.Recording = "Stardate: " + startDate2 + " Startime: " + startTime;
end


function new_data = baseline_remover(time, data, fit_type)
%   Removes baseline wander from one or more signals using a fourier curve
%   fitting type.
%
%   Inputs:
%   time - time vector for the x axis.
%   data - signal data.
%
%   Output:
%   new_data = signal data without the baseline wande

    % Preallocating variable to store baseline-removed signals
    new_data = zeros(size(data));
    
    for i = 1:size(data, 2)
        [time, data(:, i)] = prepareCurveData(time, data(:, i));
        % Set up fittype and options.
        ft = fittype(fit_type);
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Normalize = 'on';
    
        % Fit baseline model to each signal
        baselineModel = fit(time, data(:, i), ft, opts);
    
        % Remove baseline from signal
        new_data(:, i) = data(:, i) - baselineModel(time);
    end
end