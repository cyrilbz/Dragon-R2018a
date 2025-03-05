function I_norm = imNormalise(I)
% IMNORMALISE  stretches the pixel values of an integer image type to fill
% the whole range
%
%   @input: I - Image with integer type
%
%   @output: I_norm - I shifted and scaled to fill the whole range of its data type

    type = class(I);
    if(~isinteger(type))
        error('Expected data type integer but type is: %s\n' + 'Integer types in MATLAB include: int8, int16, int32, int64, uint8, uint16, uint32, and uint64.', type);
    end
    
    I_min = I - min(I, [], 'all');
    I_norm = I_min * (intmax(type)/max(I_min, [], 'all'));
end