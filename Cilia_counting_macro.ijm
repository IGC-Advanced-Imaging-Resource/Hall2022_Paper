/*
 * Macro counting cilia per segmented nucleus
 * Segments nuclei then enlarges ROI and searches within enlarged area for cilia and centrosomes
 * 
 * Expected Input:
 * The expected input for this script are non-timelapse stacks in tiff/tif/ND2/ims/oib format that are 3 channel, where channel 1 represents nuclear staining.
 * At the start of macro, user selects a parent directory where there images can be found (nested folders dealth with fine) and a directory where they wish the results to be saved
 * 
 * Description of script workflow:
 * The first channel is used to segment the nuclei (background subtraction, median filter, Huang threshold, fill holes, watershed, connectivity)
 * Each nucleus ROI is centered to view and a small window pops up for the user to indicate how many cilia they see for the central cell (could be used to count other cellular objects
 * 
 * Output:
 * An RGB image is saved with the nuclei outlines flattened
 * A .CSV files with each nucleus listed, the filename of the image it's from and how many cilia the user has said it has
 * 
 * Installation Notes: 
 * Download FIJI here - https://imagej.net/software/fiji/downloads                                                
 * How to use an ImageJ macro from Github - https://github.com/IGC-Advanced-Imaging-Resource/ImageJ_Macros
 * This script is hosted here with an example image - https://github.com/IGC-Advanced-Imaging-Resource/Hall2022_paper       
 * 
 * Written by Laura Murphy (laura.murphy@ed.ac.uk)                                                                                                                                                                                                                                                                                                                             
 * IGC Advanced Imaging Resource (https://www.ed.ac.uk/institute-genetics-cancer/facilities/advanced-imaging-resource)
 * First written: November 2021  Last updated April 2022                                                                                                       				
*/

//--------------------------------//-------------------------------------------------------------------------------------------------------
//-- Part 0: Preparation steps: get directories from user, set up arrays to store results that are being added to and some other setup
//--------------------------------//-------------------------------------------------------------------------------------------------------

// -- Getting input and output folders from user input 
inputFolder = getDirectory("Choose the folder containing your images");
outputFolder = getDirectory("Choose the folder where you want to save your results")

count = 0;

dirList = newArray();
dirList = getFileTree(inputFolder, dirList);

print("There are " + dirList.length + " files to be processed:");

Filename = newArray;
Nuclei = newArray;
Cilia = newArray;

//--------------------------------//-----------------------------------------------------------------------------------
//-- Part 1: Opening images, creating sum intensity projection, applying LUTs and renaming images 
//--------------------------------//-----------------------------------------------------------------------------------
	
for (i=0; i<dirList.length; i++){ 
	path = dirList[i];
	run("Bio-Formats Macro Extensions");
	run("Bio-Formats Importer", "open=[&path] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");           
	
	print("The active file being processed is: \n" + dirList[i]); 

	fileName = getTitle();
	imgName = File.nameWithoutExtension();
	getDimensions(width, height, channels, slices, frames);

	// Duplicates image to preserve original
	Stack.setChannel(1);
	run("Blue"); 
	Stack.setChannel(3);
	run("Red"); 
	run("Z Project...", "projection=[Max Intensity]");
	rename("Projection");
	selectWindow(fileName);
	run("Close");
	selectWindow("Projection");
	run("Duplicate...", "duplicate");	
	rename(imgName);
	selectWindow("Projection");
	run("Duplicate...", "duplicate");	
	rename("Active");

	// Splits channels and renames based on signal
	run("Split Channels");	
	selectWindow("C1-Active");
	rename("Nuclei");
	selectWindow("C2-Active");
	rename("Cilia");
	selectWindow("C3-Active");
	rename("Tubulin");

//--------------------------------//-----------------------------------------------------------------------------------
//-- Part 2: Isolating nuclei and adding to ROI
//--------------------------------//-----------------------------------------------------------------------------------

    // Select nuclei
	selectWindow("Nuclei");
	run("Duplicate...", "title=Mask_nuclei");
	run("Subtract Background...", "rolling=50 sliding");
	run("Enhance Contrast...", "saturated=0 normalize equalize");
	run("Median...", "radius=5");
	setAutoThreshold("Huang dark");
	run("Convert to Mask");
	run("Fill Holes");
	run("Watershed");
	run("Analyze Particles...", "size=100-Infinity exclude add");
	
	nuclei = roiManager("Count");

	for(j = 0; j<nuclei; j++){
		roiManager("Select", j);
		roiManager("Rename", "Nuclei_"+(j+1));
		
		currentNuclei = "Nuclei "+(j+1);
		currentImage = imgName;
	
		Nuclei = Array.concat(Nuclei, currentNuclei);
		Filename = Array.concat(Filename, currentImage);
		
		roiManager("Set Color", "cyan");
		selectWindow("Projection");
		roiManager("Select", j);
		run("Enlarge...", "enlarge=3");
		run("To Selection"); 
		run("Select None");	
		roiManager("deselect");

		Stack.setChannel(1);
		run("Enhance Contrast", "saturated=0.35");
		Stack.setChannel(2);
		run("Enhance Contrast", "saturated=0.1");
		Stack.setChannel(3);
		run("Enhance Contrast", "saturated=0.1");
	
		Dialog.create("Input for nuclei " + (j+1));
		Dialog.addNumber("Number of cilia for the central nucleus:", 0);
		Dialog.show();
		
		ciliaCurrent = Dialog.getNumber();
		
		Cilia = Array.concat(Cilia, ciliaCurrent);
		
		}

	selectWindow("Projection");
	run("Flatten");
	roiManager("Show all without labels");
	run("Flatten");
	saveAs("Tiff", outputFolder + imgName + ".tiff");
	
	run("Close All");
	roiManager("Reset");

	Array.show("Results", Filename, Nuclei, Cilia);
	saveAs("Results", outputFolder + File.separator + "Summary Results.csv");
		
	}

//--------------------------------//-----------------------------------------------------------------------------------
//-- Part 6: Results
//--------------------------------//-----------------------------------------------------------------------------------

Array.show("Results", Filename, Nuclei, Cilia);
saveAs("Results", outputFolder + File.separator + "Summary Results.csv");

//--------------------------------//-----------------------------------------------------------------------------------
//-- Part 7: Finish
//--------------------------------//-----------------------------------------------------------------------------------

if (isOpen("Results")) { 
	selectWindow("Results"); 
	run("Close"); 
} 

if (isOpen("Log")) { 
	selectWindow("Log"); 
    run("Close"); 
} 

setBatchMode(false);

run("Close All");
roiManager("reset");
run("Collect Garbage");

// Let user know it's done!
Dialog.create("Progress"); 
Dialog.addMessage("Macro Complete!");
Dialog.show; 

//--------------------------------//-----------------------------------------------------------------------------------
//-- Getting filenames from child folders, only searching for file types produced within facility
//--------------------------------//-----------------------------------------------------------------------------------


function getFileTree(dir , fileTree){
	list = getFileList(dir);

	for(f = 0; f < list.length; f++){
		if (matches(list[f], "(?i).*\\.(tif|tiff|nd2|lif|ndpi|mvd2|ims|oib)$"))
			fileTree = Array.concat(fileTree, dir + list[f]);
		if(File.isDirectory(dir + File.separator + list[f]))
			fileTree = getFileTree(dir + list[f],fileTree);
	}
	return fileTree;
}
