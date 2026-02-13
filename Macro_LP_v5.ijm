//Overlap measurement macro for an open three channel stack = red, green, DAPI
//Ema Bobocioiu-Caracas/Jacques Brocard (PLATIM) for Louis Picq/Antoine Marcais (CIRI) @2025
//May also be used for whole image directories using "Process/Batch/Macros"

//Initialization
run("Set Measurements...", "area mean redirect=None decimal=3");
t=getTitle();
getPixelSize(unit, pixelWidth, pixelHeight);
run("Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1 unit=pixel");
dir=getDirectory("image");

//Segmenting RED vesicles
selectWindow(t);
setSlice(1);
run("Select All");
run("Copy"); 
run("Internal Clipboard");
run("Subtract Background...", "rolling=5");
setAutoThreshold("Otsu dark");
run("Convert to Mask");
//run("Options...", "iterations=1 count=1 black do=Open");
rename("red");

//Segmenting GREEN vesicles
selectWindow(t);
setSlice(2);
run("Select All");
run("Copy");
run("Internal Clipboard");
run("Subtract Background...", "rolling=5");
setAutoThreshold("Otsu dark");
run("Convert to Mask");
//run("Options...", "iterations=1 count=1 black do=Open");
rename("green");

//Producing YELLOW image
imageCalculator("AND create", "red","green");
rename("yellow");

//Segmenting CYAN nuclei
selectWindow(t);
setSlice(3);
run("Select All");
run("Copy"); 
run("Internal Clipboard");
run("Median...", "radius=5");
setAutoThreshold("Default dark");
setAutoThreshold("Otsu dark");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Options...", "iterations=2 count=1 black pad do=Open");
run("Fill Holes");
//run("Watershed");
rename("cyan");
run("Analyze Particles...", "size=100-Infinity pixel clear add exclude");
setForegroundColor(0,0,0);

//Erasing any red, green or yellow signal from the nuclei area
selectWindow("red");
roiManager("Fill");
selectWindow("green");
roiManager("Fill");
selectWindow("yellow");
roiManager("Fill");

//Create a composite view of the segmented areas
run("Merge Channels...", "c1=red c2=green c5=cyan c7=yellow create ignore");

//Delineate cell contours
selectWindow(t);
run("Z Project...", "projection=[Max Intensity]");
run("Median...", "radius=5");
setAutoThreshold("Default dark");
setAutoThreshold("Huang dark");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Options...", "iterations=2 count=1 black pad do=Close");
run("Fill Holes");
run("Watershed");
run("Analyze Particles...", "size=500-Infinity pixel clear add exclude");
close();

//From each cell contour, measure intensity of each segmented area
selectWindow("Composite");
nROIs=roiManager("Count");
area=newArray(nROIs);
meanR=newArray(nROIs);
meanG=newArray(nROIs);
meanC=newArray(nROIs);
meanY=newArray(nROIs);
for (r=0;r<nROIs;r++){
	setSlice(1);
	roiManager("Select",r);
	getStatistics(area[r],meanR[r]);
	setSlice(2);
	roiManager("Select",r);
	getStatistics(area[r],meanG[r]);
	setSlice(3);
	roiManager("Select",r);
	getStatistics(area[r],meanC[r]);
	setSlice(4);
	roiManager("Select",r);
	getStatistics(area[r],meanY[r]);
	temp="";
	if (r<10) temp="0";
	roiManager("Rename", "Cell"+temp+r);
}
print("\n"+t);
print("Cell \t%Red \t%Green \t%Yellow \trandom \t%coloc \tArea(um²) \tRed(um²) \tGreen(um²)");	
del=0;
for (r=0;r<nROIs;r++){
	roiManager("Select",r-del);
	if (meanC[r]==0){//Delete cells with no nucleus, i.e. Cyan = 0
		roiManager("Delete");
		del++;
	}else{
		pixR=floor(area[r]*meanR[r]/255+1); //#red segmented pixels
		pixG=floor(area[r]*meanG[r]/255+1); //#green segmented pixels
		pixY=floor(area[r]*meanY[r]/255+1); //#yellow segmented pixels
		pixC=area[r]-floor(area[r]*meanC[r]/255+1); //#pixels of the cytoplasm = whole area - #cyan segmented pixels
		percentR=pixR/pixC; //%occupancy of red / cytoplasm
		percentG=pixG/pixC; //%occupancy of green / cytoplasm
		percentY=pixY/pixC; //%occupancy of yellow / cytoplasm
		randomY=percentR*percentG; //random expected #yellow segmented pixels based on % occupancy red vs green
		//Note that last printed variable is called %coloc since it corresponds to measured %occupancy of yellow / cytoplasm minus random expected occupancy
		area[r]=pixelWidth*pixelHeight*area[r];
		print(Roi.getName() +" \t"+ percentR+" \t"+ percentG+" \t"+ percentY+" \t"+ randomY+" \t"+ percentY-randomY+" \t"+ area[r]+" \t"+ percentR*area[r]+" \t"+ percentG*area[r]);	
	}
}

close(t);
selectWindow("Composite");
run("Remove Overlay");
roiManager("Deselect");
run("From ROI Manager");
selectWindow("Log");
saveAs("Txt",dir+"Results.txt");
