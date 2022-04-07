/*
 * Radial Intensity from Centrosomes
 * Macro for segmenting centrosomes and measuring intensity spreading out from them radially
 * 
 * Expected Input:
 * The expected input for this script are non-timelapse stacks in tiff/tif/ND2/ims/oib format that are 3 channel, where channel 1 represents nuclear staining, channel 3 represents centrosome staining and channel 2 is the channel of interest
 * At the start of macro, user selects a parent directory where there images can be found (nested folders dealth with fine) and a directory where they wish the results to be saved
 * 
 * Description of script workflow:
 * The third channel is used to segment the centrosome (RenyiEntropy threshold, watershed and size filtering)
 * For each centrosome, bands of a certain pixel size are created by first increase the selection size and using  the 'XOR' function with the previous selection to create a band
 * The mean intensity of the second channel is measured in the centrosome roi and for each band created
 * 
 * Output:
 * An RGB image is saved with the centrosome outlines flattened on and another with the bands flattened on
 * The full set of ROIs in ImageJ format are saved for each image
 * A .CSV is created for each image which contains all the centrosome and band mean intensities
 * One additional .CSV file contained the average of these measurements for each image processed.
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

// Get user to select folder with images for input and the output folder where they want results saved
inputFolder = getDirectory("Choose the folder containing your images");
outputFolder = getDirectory("Choose the folder where you want to save your results")

setBatchMode(true);

// Get list of images in nested folders with certain extensions  - function at bottom of macro
dirList = newArray();
dirList = getFileTree(inputFolder, dirList);

// -- Message to users about how many images will be processed
count = dirList.length;
print("There are " + count + " files to be processed");

// Set up empty arrays that will be filled during script and displayed and saved at the end
Filename = newArray;
Signal_Mean = newArray(count);
Centrosome_Mean = newArray(count);
Band_One_Mean = newArray(count);
Band_Two_Mean = newArray(count);
Band_Three_Mean = newArray(count);
Band_Four_Mean = newArray(count);
Band_Five_Mean = newArray(count);

projection = "Average Intensity"; //Choose from 'Average Intensity' and 'Sum Slices ' for measurements to be made on
radius = 0.5; //change this number to change band diameter that are created around centrosomes

//--------------------------------//-----------------------------------------------------------------------------------
//-- Part 1: Image loop begins,
//--------------------------------//-----------------------------------------------------------------------------------

// -- Using Bio-Formats to open images
for (i = 0; i < dirList.length; i++){ 
	path = dirList[i];
	run("Bio-Formats Macro Extensions");
	run("Bio-Formats Importer", "open=&path color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");           

	fileName = getTitle();
	Filename = Array.concat(Filename, fileName);
	imgName = File.nameWithoutExtension();
	getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);

	// Create max projection for segmentation channels, set channel LUTs
	Stack.setChannel(1);
	run("Blue"); 
	Stack.setChannel(3);
	run("Red"); 
	run("Z Project...", "projection=[Max Intensity]");
	rename("MProjection");

	// Create  projection for measured channel, set channel LUTs
	selectWindow(fileName);
	Stack.setChannel(1);
	run("Blue"); 
	Stack.setChannel(3);
	run("Red"); 
	run("Z Project...", "projection=[" + projection + "]");
	run("Duplicate...", "duplicate");	
	rename("Projection");
	
	selectWindow(fileName);
	run("Close");
	
	// Splits channels and renames based on signal
	selectWindow("MProjection");
	run("Duplicate...", "duplicate");	
	rename("Active");
	run("Split Channels");	
	selectWindow("C1-Active");
	rename("Nuclei");
	selectWindow("C2-Active");
	run("Close");
	selectWindow("C3-Active");
	rename("Centrosomes");

	selectWindow("Projection");
	run("Duplicate...", "duplicate");	
	rename("Active");
	run("Split Channels");	
	selectWindow("C1-Active");
	run("Close");
	selectWindow("C2-Active");
	rename("Signal");
	selectWindow("C3-Active");
	run("Close");
	
//--------------------------------//-----------------------------------------------------------------------------------
//-- Part 2: Segmenting centrosomes
//--------------------------------//-----------------------------------------------------------------------------------

	// -- Segmentation with centrosomes, may need changed (think about sum vs MIP projections, may want a different threshold and size limit)
	selectWindow("Centrosomes");
	setAutoThreshold("RenyiEntropy dark");
	run("Convert to Mask");
	run("Watershed");
	makeRectangle(15,15,getWidth-30, getHeight-30); //This ROI is created to make sure bands can be created instead of going off image edge
	run("Analyze Particles...", "size=0.1-Infinity exclude add");
	roiManager("Set Color", "magenta");
	run("Set Measurements...", "area centroid display redirect=None decimal=4");

	// -- Save image with centrosome ROIs
	selectWindow("MProjection");
	Stack.setChannel(1);
	run("Enhance Contrast...", "saturated=0.3 normalize");
	Stack.setChannel(2);
	run("Enhance Contrast...", "saturated=0.3 normalize");
	Stack.setChannel(3);
	run("Enhance Contrast...", "saturated=0.3 normalize");
	run("RGB Color");
	roiManager("Show All without labels");
	run("Flatten");
	saveAs("Tiff", outputFolder + File.separator + imgName + "_Centrosomes.tif");

	// -- Print message to tell user how many centrosomes have been detected
	centrosomes = roiManager("Count");
	print(imgName + " has a pixel size of " + pixelWidth + " and has " + centrosomes + " centrosomes"); 

//--------------------------------//-----------------------------------------------------------------------------------
//-- Part 3: Create bands from segmentation and measure them on the signal channel
//-- Currently five bands around centrosome are created
//--------------------------------//-----------------------------------------------------------------------------------

	selectWindow("Signal");
	rename(imgName);
	getStatistics(area, Signal_Mean[i], min, max, std, histogram);
	
	// -- Create arrays to store results
	Label = newArray(centrosomes);
	Centrosome = newArray(centrosomes);
	Band_One = newArray(centrosomes);
	Band_Two = newArray(centrosomes);
	Band_Three = newArray(centrosomes);
	Band_Four = newArray(centrosomes);
	Band_Five = newArray(centrosomes);

	// -- Band creation time
	for(j = 0; j < centrosomes; j++){
		roiManager("Select", j);	
		roiManager("Rename", "Centrosome_"+(j+1));
		Roi.setStrokeColor("Magenta");
		Label[j] = "Centrosome_"+(j+1);
		getStatistics(area, Centrosome[j], min, max, std, histogram);

		roiManager("Select", j);
		run("Enlarge...", "enlarge=" + radius);
		Roi.setName("Centrosome_" + (j+1) + ": Circle_1");
		Roi.setStrokeColor("Orange");
		roiManager("Add");
		roiManager("Select", newArray(j, roiManager("Count")-1));
		roiManager("XOR");
		Roi.setName("Centrosome_" + (j+1) + ": Band_1");
		Roi.setStrokeColor("Yellow");
		roiManager("Add");
		roiManager("Select", roiManager("Count")-1);
		getStatistics(area, Band_One[j], min, max, std, histogram);

		roiManager("Select", roiManager("Count")-2);
		run("Enlarge...", "enlarge=" + radius);
		Roi.setName("Centrosome_" + (j+1) + ": Circle_2");
		Roi.setStrokeColor("Orange");
		roiManager("Add");
		roiManager("Select", newArray(roiManager("Count")-3, roiManager("Count")-1));
		roiManager("XOR");
		Roi.setName("Centrosome_" + (j+1) + ": Band_2");
		Roi.setStrokeColor("Yellow");
		roiManager("Add");
		roiManager("Select", roiManager("Count")-1);
		getStatistics(area, Band_Two[j], min, max, std, histogram);

		roiManager("Select", roiManager("Count")-2);
		run("Enlarge...", "enlarge=" + radius);
		Roi.setName("Centrosome_" + (j+1) + ": Circle_3");
		Roi.setStrokeColor("Orange");
		roiManager("Add");
		roiManager("Select", newArray(roiManager("Count")-3, roiManager("Count")-1));
		roiManager("XOR");
		Roi.setName("Centrosome_" + (j+1) + ": Band_3");
		Roi.setStrokeColor("Yellow");
		roiManager("Add");
		roiManager("Select", roiManager("Count")-1);
		getStatistics(area, Band_Three[j], min, max, std, histogram);

		roiManager("Select", roiManager("Count")-2);
		run("Enlarge...", "enlarge=" + radius);
		Roi.setName("Centrosome_" + (j+1) + ": Circle_4");
		Roi.setStrokeColor("Orange");
		roiManager("Add");
		roiManager("Select", newArray(roiManager("Count")-3, roiManager("Count")-1));
		roiManager("XOR");
		Roi.setName("Centrosome_" + (j+1) + ": Band_4");
		Roi.setStrokeColor("Yellow");
		roiManager("Add");
		roiManager("Select", roiManager("Count")-1);
		getStatistics(area, Band_Four[j], min, max, std, histogram);

		roiManager("Select", roiManager("Count")-2);
		run("Enlarge...", "enlarge=" + radius);
		Roi.setName("Centrosome_" + (j+1) + ": Circle_5");
		Roi.setStrokeColor("Orange");
		roiManager("Add");
		roiManager("Select", newArray(roiManager("Count")-3, roiManager("Count")-1));
		roiManager("XOR");
		Roi.setName("Centrosome_" + (j+1) + ": Band_5");
		Roi.setStrokeColor("Yellow");
		roiManager("Add");
		roiManager("Select", roiManager("Count")-1);
		getStatistics(area, Band_Five[j], min, max, std, histogram);
	}

	// -- Save image with all ROIs drawn on
	selectWindow("MProjection");
	Stack.setChannel(1);
	run("Enhance Contrast...", "saturated=0.3 normalize");
	Stack.setChannel(2);
	run("Enhance Contrast...", "saturated=0.3 normalize");
	Stack.setChannel(3);
	run("Enhance Contrast...", "saturated=0.3 normalize");
	run("RGB Color");
	roiManager("Show All without labels");
	roiManager("Save", outputFolder + File.separator + imgName + "_ROIset.zip");
	run("Flatten");
	saveAs("Tiff", outputFolder + File.separator + imgName + "_Outlines.tif");

	run("Close All");
	roiManager("Reset");

	// -- Save the image specific results
	Array.show("Mean Intensity", Label, Centrosome, Band_One, Band_Two, Band_Three, Band_Four, Band_Five);
	saveAs("Results", outputFolder + File.separator + imgName + "_Results.csv");
	run("Close");

	// -- Take average to save overall mean results across images
	Array.getStatistics(Centrosome, min, max, Centrosome_Mean[i], stdDev);
	Array.getStatistics(Band_One, min, max, Band_One_Mean[i], stdDev);
	Array.getStatistics(Band_Two, min, max, Band_Two_Mean[i], stdDev);
	Array.getStatistics(Band_Three, min, max, Band_Three_Mean[i], stdDev);
	Array.getStatistics(Band_Four, min, max, Band_Four_Mean[i], stdDev);
	Array.getStatistics(Band_Five, min, max, Band_Five_Mean[i], stdDev);	

	run("Collect Garbage");	
}

//--------------------------------//-----------------------------------------------------------------------------------
//-- Part 4: Finishing off steps now images are processed
//--------------------------------//-----------------------------------------------------------------------------------

// -- Save the average results for each image to user directed output folder now all images are processed
Array.show("Mean Intensity", Filename , Signal_Mean, Centrosome_Mean, Band_One_Mean, Band_Two_Mean, Band_Three_Mean, Band_Four_Mean, Band_Five_Mean );
saveAs("Results", outputFolder + File.separator + "Mean intensity for centrosomes bands per image.csv");
run("Close");

// -- Save log window for user (tells them when it was run, what kind of projection they choose and what diameter they used for band creation
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
print("This macro was run on " + dayOfMonth + "/" + month+1 + "/" + year + ": " + hour + ":" + minute + "\nUsing a radius of " + radius + " " + unit + "\nOn " + projection + " projections");
selectWindow("Log");
saveAs("Text", outputFolder + "Details of macro run");
print("\\Clear");

// -- Let user know the macro is complete.
Dialog.create("Progress");
Dialog.addMessage("Macro Complete!");
Dialog.show;

//--------------------------------//-----------------------------------------------------------------------------------
//-- Epilogue: Functions
//--------------------------------//-----------------------------------------------------------------------------------

function getFileTree(dir , fileTree){
	list = getFileList(dir);

	for(f = 0; f < list.length; f++){
		if (matches(list[f], "(?i).*\\.(tif|tiff|nd2|ims|oib)$"))
			fileTree = Array.concat(fileTree, dir + list[f]);
		if(File.isDirectory(dir + File.separator + list[f]))
			fileTree = getFileTree(dir + list[f],fileTree);
	}
	return fileTree;
}
