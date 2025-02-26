function store_basefunctions(shape_base, shape_bounds, ...
    test_base, test_bounds, turb_field, diff_sys, conv_sys, ...
    equation_type)
%STORE_BASEFUNCTIONS store the base functions
rom_base.test_bounds = test_bounds;
rom_base.test_base = test_base;
rom_base.shape_bounds = shape_bounds;
rom_base.shape_base = shape_base;
rom_base.turb_base = turb_field;
rom_base.diff_system = diff_sys;
rom_base.conv_system = conv_sys;
rom_base.equation_type = equation_type;
save('rom_base.mat', 'rom_base');
end

