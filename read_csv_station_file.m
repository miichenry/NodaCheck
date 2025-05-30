function [stat_list, stat_lon, stat_lat] = read_csv_station_file(filename, dataLines)
%IMPORTFILE Import data from a text file
%  [STAT_LIST, STAT_LON, STAT_LAT] = IMPORTFILE(FILENAME) reads data
%  from text file FILENAME for the default selection.  Returns the data
%  as column vectors.
%
%  [STAT_LIST, STAT_LON, STAT_LAT] = IMPORTFILE(FILE, DATALINES) reads
%  data for the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  [stat_list, stat_lon, stat_lat] = importfile("/media/savardg/sandisk4TB/matlab-swant/scripts/synthetic_tests_to_design_network/NANT_stations_example_D20km_N18.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 25-May-2023 15:14:10

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["stat_list", "stat_lon", "stat_lat", "Var4", "Var5"];
opts.SelectedVariableNames = ["stat_list", "stat_lon", "stat_lat"];
opts.VariableTypes = ["double", "double", "double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var4", "Var5"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var4", "Var5"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
stat_list = tbl.stat_list;
stat_lon = tbl.stat_lon;
stat_lat = tbl.stat_lat;
end