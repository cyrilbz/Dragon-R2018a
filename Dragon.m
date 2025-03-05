function net_props = Dragon(path)
% DRAGON  Extracts the properties of an actin network in 3D using images
%
%   @input: path - the file path in which all of the data resides. Make
%           sure to include the trailing backslash or it will not work.
%
%   @output net_props - a structure containing all of the data extracted
%                       from the image stack
%
%   Note that the input filenames are stored as variables to make things
%   simple to change
%
%   Default Filenames:
%       raw_stack_filename  = "actin_stack.tif";         - REQUIRED INPUT
%       projection_filename = "actin_original.tif";
%       mask_filename       = "mask.tif";                - REQUIRED INPUT
%       variables_filename  = "netProps.mat";
%
%   The ratio between pixel width and voxel depth is specified by z_scale
%   so will depend on your microscopy setup (in getParams.m)

    % All filenames used in this function
    raw_stack_filename  = "actin_stack.tif"; %existing file needed to run
    projection_filename = "actin_original.tif"; %gets written
    mask_filename       = "mask.tif"; %existing file needed to run
    variables_filename  = "netProps.mat"; %gets written
    
    filename = char(path + raw_stack_filename);
    im_stack = grabImageStack(filename);
    im_raw = maxProjection(im_stack);
    imwrite(im_raw, char(path + projection_filename));
    im_raw = cytoFilter(char(path));
    

    % Grab 2D Network
    im_mask = grabCellMask(char(path+mask_filename), im_raw);
    net_props = getNetworkProperties(im_raw, im_mask);


    % Extend to 3D
    im_stack = rotateStack(im_stack, rad2deg(-net_props.orientation));
    p = getParams();
    %z_scale = 5.5;       %ratio of step size in z to pixel width
    net_props = angleZ(net_props, im_stack, p.z_scale);

    
    % Plot, graph and save the results
    save(char(path + variables_filename), 'net_props');    % Store all of the raw data in .mat format
    %outDat(net_props, path);
    
end

function outDat(np, path)
    datafile = fopen(path + "netProps.dat", 'w');
    fprintf(datafile, "BranchRatio\tFilAngXY\tFilAngZ\tBranchAngXY\tBranchAngZ\tFilLen\tStructSize\tStructFreq\tskelDensity\tBPDensity\tavgDeviation\tavgCurvature\tavgCurvatureSigned\tavgFilWidth\r\n");
    fprintf(datafile, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", np.branchRatio, np.avgFilXY, np.avgFilZ, np.avgBranchXY, np.avgBranchZ, np.avgLen, np.Structures.mean, np.Structures.freq, np.skelDensity, np.cellBPDensity, mean(np.deviation(:), 'omitnan'), mean(np.curvature(:), 'omitnan'), mean(np.curvatureSigned(:), 'omitnan'), mean(np.filWidth(:), 'omitnan'));
    fclose(datafile);            
end
