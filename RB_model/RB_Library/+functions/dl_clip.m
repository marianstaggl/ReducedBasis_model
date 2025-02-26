function clipped_value = dl_clip(value, lower, upper)
%DL_CLIP clip a dl array between lower and upper limit with a relu function
clipped_value = relu(value - lower) - relu(value - upper);
end

