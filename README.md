# DRAGoN Summary
The DRAGoN algorithm is a series of Matlab scripts which segments a filamentous network from a fluorescence microscopy image (or stack of images). The centreline of the network is estimated and replaced with labelled pixel-thick lines, which are then analysed to provide several mathematical, physical and biological properties.

The scripts are split up in such a way that they may be useful by themselves, or perhaps as a subset. This also enables the Matlab help function to quickly provide a summary of how each script works. If you don't care about that; however, you can just call the function in `Dragon.m` and everything will be called for you.

To learn about the development, testing, and outcomes, please read our paper here: [AUTOMATIC EXTRACTION OF ACTIN NETWORKS IN PLANTS](https://doi.org/10.1101/2023.01.18.524528 "preprint link")



# Table of Contents

- [Requirements](#requirements)
  - [File Requirements](#file-requirements)
  - [MATLAB Requirements](#matlab-requirements)
    - [Compatibility Details](#compatibility-details)
- [DRAGoN Function Information](#dragon-function-information)
  - [Algorithm Parameters and Behaviour](#algorithm-parameters-and-behaviour)
- [DRAGoN Tutorial](#dragon-tutorial)
  - [First Steps](#first-steps)
  - [Further Testing](#further-testing)

# Requirements
### File Requirements
You will need the image stack and the mask in the given path or the code will not run. ***Ensure their names match the filenames in Dragon.m***

**Default Input Filenames:**
```
raw_stack_filename = "actin_stack.tif";
mask_filename      = "mask.tif";
```

**Default Output Filenames:**
```
projection_filename = "actin_original.tif";
variables_filename  = "netProps.mat";
```


### MATLAB Requirements
**Toolboxes:**
  + Image Processing Toolbox
  + Statistics and Machine Learning Toolbox

Minimum Version: r2018b
***Recommended Version: r2022a***

#### Compatibility Details
For function compatibility *and* consistent behaviour, you will need *at least* MATLAB r2018b as well as the **Image Processing Toolbox** and the **Statistics and Machine Learning Toolbox**. Previous versions may run, but algorithm behavioural changes occurred and results may be sub-optimal. Tested and developed on r2019a, with some recent testing on r2022a also working.

# DRAGoN Function Information
```
@input: path - the file path in which all of the data resides
@output net_props - a structure containing all of the data extracted
		    from the image stack

Note that the input filenames are stored as variables to make things simple to change. 

Default Filenames:
    raw_stack_filename  = "actin_stack.tif";
    projection_filename = "actin_original.tif";
    mask_filename       = "mask.tif";
    variables_filename  = "netProps.mat";

The ratio between pixel width and voxel depth is specified by z_scale so will depend on your microscopy setup.
This parameter, along with the rest, is found in getParams.m
```

***Please note that the path being passed to `Dragon.m` should include the trailing slash(s) (type dependant on your OS), or the string concatenation process may not work.***

### Algorithm Parameters and Behaviour
All of the parameters are stored in `getParams.m` with a description of the effect of each one. `prctile_thresh` is the most sensitive and will likely need tweaking for your setup. `z_scale` is the ratio between the size of your pixel in 2D and the distance between your stacked images, make sure to adjust this where needed.


# DRAGoN Tutorial
### First Steps
1. Download the complete repository, including all .m files and the Examples_ActinHypocotyls folder
2. Using the navigation pane on the left side of MATLAB, find the folder containing all the .m files
 
	a. Right click on the folder, and 'add to path'
	
	b. This ensures that MATLAB can find all the script files when running, no matter which folder you are in
	
3. Navigate to the Examples_ActinHypocotyls folder, and enter it

	a. You should now see 3 folders (Arp2, Col0, Triple) and a .mat file
	
	b. Double click the .mat file to import the path lists into MATLAB (variable called `lst`)

4. Test that everything is working and in the correct locations with the following command: `Dragon(lst(1));`

	a. If everything is working, it should take a little under a minute to finish, and `ans` should contain the output.
	
	b. Typing in `imshow(labeloverlay(imadjust(ans.imRaw, stretchlim(ans.imRaw, 0)), ans.skelLabel));` will provide the original image with the network overlayed in various colours.



### Further Testing
If everything in [First Steps](#first-steps) worked, we can begin to run more tests. To run everything in one go and see those overlays, you can use:
```
for idx = 1:length(lst)
	Dragon(lst(idx));
	figure(idx);
	imshow(labeloverlay(imadjust(ans.imRaw, stretchlim(ans.imRaw, 0)), ans.skelLabel));
end
```

Inside `Dragon.m`, there is a function which outputs some of our key measures, and averages, into a file called `netProps.dat`, a tab delimited text file. While this is good for single runs, collating these for a large data series can be annoying. Combining all of the output structures into a structure array, can keep your data in one place and make it easier to loop over it and extract what is needed for downstream analysis. To do this, we edit the above to look something like this:
```
netProps = [];
for idx = 1:length(lst)
	ans = Dragon(lst(idx));
	netProps = [NetProps; ans];
	figure(idx);
	imshow(labeloverlay(imadjust(ans.imRaw, stretchlim(ans.imRaw, 0)), ans.skelLabel));
end
```
Now you have all of the output in an array of structures - each row a new image and each column a series of measurements or additional structures.
