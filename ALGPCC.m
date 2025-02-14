function I_ALGPCC = ALGPCC(raw, mosaic, mask)

    % Initialize variables
    [row, col, ch] = size(mosaic);
    
    persistent filters_sp filters_po
    if isempty(filters_sp) || isempty(filters_po)
        filters_sp = cat(3, ...
            [-1, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0], ...
            [0, 0, -1, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0], ...
            [0, 0, 0, 0, -1; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0], ...
            [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; -1, 0, 1, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0], ...
            [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, -1; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0], ...
            [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 0, 0; -1, 0, 0, 0, 0], ...
            [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 0, 0; 0, 0, -1, 0, 0], ... 
            [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, -1] ...
        );
        filters_po = cat(3, ...
            [0, 0, 0; -1, 0, 1; 0, 0, 0], ...
            [0, -1, 0; 0, 0, 0; 0, 1, 0], ...
            [-1, 0, 0; 0, 0, 0; 0, 0, 1], ...
            [0, 0, -1; 0, 0, 0; 1, 0, 0] ...
        );
    end
    LG_sp_tmp = arrayfun(@(i) abs(imfilter(raw, filters_sp(:,:,i), 'replicate')), 1:size(filters_sp, 3), 'UniformOutput', false);
    LG_po_tmp = arrayfun(@(i) abs(imfilter(raw, filters_po(:,:,i), 'replicate')), 1:size(filters_po, 3), 'UniformOutput', false);
    
    LG_sp = 1/8 .* (LG_sp_tmp{1} + LG_sp_tmp{2} + LG_sp_tmp{3} + LG_sp_tmp{4} + LG_sp_tmp{5} + LG_sp_tmp{6} + LG_sp_tmp{7} + LG_sp_tmp{8});
    LG_po = 1/3 .* (LG_po_tmp{1} + LG_po_tmp{2} + 1/2 .* (LG_po_tmp{3} + LG_po_tmp{4}));
    LG = 1/2 .* (LG_sp + LG_po);
    LG_weight = 1./ (LG + eps);

    % Initialization
    Filter_bilinear = [1, 2, 1; 2, 4, 2; 1, 2, 1] / 4;
    I_initial = imfilter(LG_weight .* mosaic, Filter_bilinear, 'replicate') ./ ...
                imfilter(LG_weight .* mask, Filter_bilinear, 'replicate');

    % Calculate box filter size
    h = 7; v = 7;
    N = boxfilter(ones(row, col), h, v);

    % Calculate means and variances
    means = boxfilter(I_initial, h, v) ./ N;
    vars = abs(boxfilter(I_initial.^2, h, v) ./ N - means.^2) + eps;

    % Calculate NCC and APCC_weight
    channels = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    num_channels = size(channels, 1);
    NCC = zeros(row, col, num_channels);
    APCC_weight = zeros(row, col, num_channels);
    for i = 1 : num_channels
        ch1 = channels(i, 1);
        ch2 = channels(i, 2);
        NCC_tmp = (boxfilter(I_initial(:,:,ch1).*I_initial(:,:,ch2), h, v) ./ N - means(:,:,ch1).*means(:,:,ch2)) ./ sqrt(vars(:,:,ch1).*vars(:,:,ch2));
        NCC(:,:,i) = max(min(NCC_tmp, 1), -1);
        APCC_weight(:,:,i) = 1 ./ (sqrt((vars(:,:,ch1) + vars(:,:,ch2)) .* (1 - NCC(:,:,i).^2)) + eps);
    end
   
    % Calculate I_ALGPCC
    I_ALGPCC = zeros(row, col, ch);
    I_array = [1,2,3,4,1,2,3];
    APCC_array = [1,2,3; 4,5,1; 6,2,4; 3,5,6];
    for k = 1 : ch
        I_diff = mosaic(:,:,k) - I_initial(:,:,I_array(k+1:k+3)) .* mask(:,:,k);
        I_diff = imfilter(LG_weight .* I_diff, Filter_bilinear, 'replicate') ./ ...
                 imfilter(LG_weight .* mask(:,:,k), Filter_bilinear, 'replicate');
        I_combined = I_initial(:,:,I_array(k+1:k+3)) + I_diff;
        
        I_ALGPCC(:,:,k) = sum(APCC_weight(:,:,APCC_array(k,:)) .* I_combined, 3) ./ ...
                        sum(APCC_weight(:,:,APCC_array(k,:)), 3);
    end
    I_ALGPCC = max(min(I_ALGPCC, 65535), 0);
end