function [shape_funs, shape_bounds, test_funs, test_bounds, turb_field, ...
    diff_system, conv_system, rom_equations] = load_base_functions(varargin)
%LOAD_BASE_FUNCTIONS Summary of this function goes here
%   Detailed explanation goes here

if isempty(varargin), filename = 'rom_base.mat';
else, filename = varargin{1}; end

full_base = load(filename);
shape_funs = full_base.rom_base.shape_base;
shape_bounds = full_base.rom_base.shape_bounds;
test_funs = full_base.rom_base.test_base;
test_bounds = full_base.rom_base.test_bounds;
turb_field = full_base.rom_base.turb_base;
diff_system = full_base.rom_base.diff_system;
conv_system = full_base.rom_base.conv_system;
rom_equations = full_base.rom_base.equation_type;
end

