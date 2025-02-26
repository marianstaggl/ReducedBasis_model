function [train_data, test_data] = load_solution_data(main_path, varargin)
%LOAD_DATA Summary of this function goes here
%   Detailed explanation goes here
ip = inputParser();
ip.addParameter('ext', '');
ip.parse(varargin{:});

c_ext = ip.Results.ext;
if ~isempty(c_ext), c_ext = ['_', c_ext]; end

train_data = load(fullfile(main_path, ['all_data_train', c_ext,'.mat']));
train_data = train_data.all_data;

test_data = load(fullfile(main_path, ['all_data_test', c_ext,'.mat']));
test_data = test_data.all_data;
end

