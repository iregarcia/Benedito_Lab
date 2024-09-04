//This macro will identify all cells within a given area, and decide which are positive for given reporters/markers.
//Below we provide a description of channels used for reporters/markers/segmentation for each image:

//Fig.1d: IsoB4 / Cherry-iChr (Reporter1) / MbYFP (Reporter2) / none (Marker) / ERG (Segmentation)
//Fig.1i,k, Fig.2f, Fig.3e, Sup.Fig.2a: IsoB4 / MbTomato (Reporter1) / MbYFP (Reporter2) / none (Marker) / ERG (Segmentation)
//Fig.2f: IsoB4 / MbTomato (Reporter1) / MbYFP (Reporter2) / none (Marker) / ERG (Segmentation)
//Fig.5a, Sup.Fig.5c,d: any channel with mean threshold at 0 / none (Reporter1) / none (Reporter2) / Ki67 (Marker) / DAPI (Segmentation)
//Fig.6g: IsoB4 / MbTomato (Reporter1) / MbYFP (Reporter2) / H2B-V5 (Marker) / ERG (Segmentation)
//Sup.Fig.2b: IsoB4 / MbTomato (Reporter1) / MbYFP (Reporter2) / Ki67 (Marker) / No Segmentation
//Sup.Fig.2c: IsoB4 / MbTomato (Reporter1) / MbYFP (Reporter2) / p21 (Marker) / No Segmentation


// Dialog window for parameters choice ----------------
yesno = newArray("no", "yes");
ftypes = newArray("tif", "jpg");
chnumb = newArray("1","2","3","4","5", "6", "none");
Dialog.create("Channel definition");
Dialog.addMessage("This macro will identify all cells within a given area, and decide which are positive for given reporters/markers.");
Dialog.addMessage("Define below which signals are present in each channel.");
Dialog.addChoice("IsoB4 (to define area of analysis): channel", chnumb);
Dialog.addChoice("Reporter1: channel", chnumb);
Dialog.addChoice("Reporter2: channel", chnumb);
Dialog.addChoice("Marker: channel", chnumb);
Dialog.addChoice("Do you want to segment based on a given channel", yesno);
Dialog.addChoice("Segmentation: channel", chnumb);


Dialog.show();

//chn = Dialog.getChoice;
IsoB4 = Dialog.getChoice;
Rep1c = Dialog.getChoice;
Rep2c = Dialog.getChoice;
Markc = Dialog.getChoice;
has=Dialog.getChoice;
Erg=Dialog.getChoice;


// End of Dialog --------------------------------------	
	
//open all images within a given folder
//path = getDirectory("Choose folder with images");
masks = getDirectory("Select folder to save results");
//index = getFileList(path);
//number = lengthOf(index);

run("Clear Results");
run("Set Measurements...", "  redirect=None decimal=3");


requires("1.52o");
	
//define name for storing resulting files
name = File.nameWithoutExtension;
	
//Prepare roi Manager(reset) and original image (guarantee that it is 8bit)
roiManager("Reset");
if(bitDepth() != 8)run("Clear Results");
run("8-bit");

//run("Set Scale...", "distance=0 known=0 unit=pixel");
	
oriID = getImageID();

//Area of interest for analysis, vein and artery
waitForUser("Draw a polygon to define the area to be analysed. Then click OK");
roiManager("Add"); //index 0 for total ROI for A-V analysis
roiManager("Select",roiManager("Count")-1);
if (selectionType!=2){
	exit("Polygon required");
}

run("Select All");
roiManager("Add");
roiManager("Add");
roiManager("Add");
roiManager("Add");


						
//Identify IsoB4+ area
selectImage(oriID);
run("Select None");
run("Duplicate...", "duplicate channels=IsoB4");
run("Subtract Background...", "rolling=300");
run("Gaussian Blur...", "sigma=2");
getStatistics(bla, meanErg, bla, bla, stdErg);
run("Threshold..."); //open threshold box
waitForUser("With threshold set to 'Mean' mode, adjust threshold and then click OK. If you want to analyse the whole area, place threshold at 0.");
getThreshold(lowerE, upperE);
setThreshold(lowerE,255);
run("Convert to Mask");
run("Create Selection");
roiManager("Add");
close();

//find all nuclei within the IsoB4 area independently of being Rep1+, Rep2+ or dual
selectImage(oriID);
run("Select None");
run("Duplicate...", "duplicate channels=Rep1c");
rename("R1");
run("Subtract Background...", "rolling=300");
run("Gaussian Blur...", "sigma=1");
getStatistics(bla, meanR1, bla, bla, std);
run("Grays");
run("Threshold..."); //open threshold box
waitForUser("With threshold set to 'Mean' mode, adjust threshold and then click OK");
getThreshold(lowerR1, upperR1);
setThreshold(lowerR1,255);
if(has=="no"){
run("Find Maxima...");}
run("Convert to Mask");
run("Watershed");
run("Invert");
rename("Rep1");
roiManager("Select", roiManager("Count")-1);//select IsoB4 selection
setBackgroundColor(255, 255, 255);
run("Clear Outside"); //clear all Rep1+ signals outside of Isob4 area
//selectWindow("R1");
//close();

selectImage(oriID);
run("Select None");
run("Duplicate...", "duplicate channels=Rep2c");
rename("R2");
run("Subtract Background...", "rolling=300");
run("Gaussian Blur...", "sigma=1");
getStatistics(bla, meanR2, bla, bla, std);
if(Rep2c=="none"){
setThreshold(255, 255);	
}
else{
run("Grays");	
run("Threshold..."); //open threshold box
waitForUser("With threshold set to 'Mean' mode, adjust threshold and then click OK");
}
getThreshold(lowerR2, upperR2);
setThreshold(lowerR2,255);
if(has=="no"){
run("Find Maxima...");}
run("Convert to Mask");
run("Watershed");
run("Invert");
rename("Rep2");
roiManager("Select", roiManager("Count")-1);//select IsoB4 selection
setBackgroundColor(255, 255, 255);
run("Clear Outside"); //clear all Rep2+ signals outside of Isob4 area
//selectWindow("R2");
//close();
roiManager("Select", roiManager("Count")-1);
roiManager("Delete");

	if (has=="no"){
	imageCalculator("AND create", "Rep1","Rep2"); //add the signals of Rep1+ nuclei and Rep2+ nuclei together
	rename("Nuc");
	run("Invert");
	run("Watershed");
	}else{
	selectImage(oriID);
	run("Select None");
	run("Duplicate...", "duplicate channels=Erg");
	run("Subtract Background...", "rolling=300");
	run("Gaussian Blur...", "sigma=1");
	getStatistics(bla, meanErg, bla, bla, stdErg);
	run("Threshold..."); //open threshold box
	waitForUser("With threshold set to 'Mean' mode, adjust threshold and then click OK");
	getThreshold(lowerERG, upperERG);
	run("Find Maxima...");
	//run("Convert to Mask");
	run("Watershed");	
	rename("Nuc");
	}
	selectWindow("Rep1");
	close();
	selectWindow("Rep2");
	close();

selectImage(oriID);
run("Select None");
run("Duplicate...", "duplicate channels=Markc");
run("Subtract Background...", "rolling=300");
run("Gaussian Blur...", "sigma=1");
getStatistics(bla, meanM, bla, bla, std);
if(Markc=="none"){
setThreshold(255, 255);	
}
else{
run("Threshold..."); //open threshold box
waitForUser("With threshold set to 'Mean' mode, adjust threshold and then click OK");
}
getThreshold(lowerM, upperM);
MarkFactor = lowerM/meanM;
close();

//eliminate all nuclei outside ROI of interest
selectWindow("Nuc");
run("Select None");
roiManager("Select", 0); // to identify ROI between artery and vein
setBackgroundColor(255, 255, 255);
run("Analyze Particles...", "size=30-2000 add");
countCell = roiManager("Count")-6;
run("Clear Outside");



//Get centroid for each nuclei and store info if Rep1+, Rep2+, or both, and if marker positive or not and draw selections
centerX = newArray(roiManager("Count")-6); //x coordinates for all ECs
centerY = newArray(roiManager("Count")-6); //y coordinates for all ECs

Rep1 = newArray(roiManager("Count")-6); //reports if Reporter1 positive or not
Rep2 = newArray(roiManager("Count")-6); //reports if Reporter2 positive or not
Mark = newArray(roiManager("Count")-6); //reports if marker positive or not

MeanRep1 = newArray(roiManager("Count")-6); //reports mean intensity of Reporter1
MeanRep2 = newArray(roiManager("Count")-6); //reports mean intensity of Reporter2
MeanMark = newArray(roiManager("Count")-6); //reports mean intensity of Marker

Region = newArray(roiManager("Count")-6); //reports to which region in retina EC belongs
ECarea = newArray(roiManager("Count")-6); //reports EC area
ECround = newArray(roiManager("Count")-6); //reports EC roundness


setBatchMode(true);

selectImage(oriID);
run("Select None");
run("Duplicate...", "duplicate channels=Rep1c");
run("Subtract Background...", "rolling=300");
run("Gaussian Blur...", "sigma=1");
rename("Rep1");

selectImage(oriID);
run("Select None");
run("Duplicate...", "duplicate channels=Rep2c");
run("Subtract Background...", "rolling=300");
run("Gaussian Blur...", "sigma=1");
rename("Rep2");

selectImage(oriID);
run("Select None");
run("Duplicate...", "duplicate channels=Markc");
run("Subtract Background...", "rolling=300");
run("Gaussian Blur...", "sigma=1");
rename("Mark");


selectImage(oriID);
run("Select None");
run("Duplicate...", "duplicate");
rename("Temp");
nBins = 256;

// Dialog window for cutoff limits ----------------
//yesno = newArray("no", "yes");
//ftypes = newArray("tif", "jpg");
chnumb_T = newArray("0","10","20","30","40","50", "60", "70", "80","90","100");
chnumb_M = newArray("0","10","20","30","40","50", "60", "70", "80","90","100");
Dialog.create("Cutoff limits");
Dialog.addMessage("To determine colocalization of reporter/marker with segmented cells, define the maximal % of pixels allowed");
Dialog.addMessage("to have an intensity smaller than the defined threshold.");
Dialog.addChoice("Maximal % of negative pixels in Reporter 1", chnumb_T);
Dialog.addChoice("Maximal % of negative pixels in Reporter 2", chnumb_T);
Dialog.addChoice("Maximal% of negative pixels in Marker", chnumb_M);
Dialog.show();

//chn = Dialog.getChoice;
Cutoff_R1 = Dialog.getChoice;
Cutoff_R2 = Dialog.getChoice;
Cutoff_M = Dialog.getChoice;


// End of Dialog --------------------------------------	


//get area and intensities of cells in um2
for(j=6;j<roiManager("Count");j++){
	selectWindow("Rep1");
	run("Select None");
	roiManager("Select", j)
	List.setMeasurements();
	ECarea[j-6] = List.getValue("Area");
	getStatistics(bla, meanR1, bla, bla, std);
	MeanRep1[j-6]=meanR1;
	selectWindow("Rep2");
	run("Select None");
	roiManager("Select", j)
	getStatistics(bla, meanR2, bla, bla, std);
	MeanRep2[j-6]=meanR2;
	selectWindow("Mark");
	run("Select None");
	roiManager("Select", j)
	getStatistics(bla, meanM, bla, bla, std);
	MeanMark[j-6]=meanM;
}
run("Set Scale...", "distance=0 known=0 unit=pixel");


//get roundness and coordinates in px of cells
for(j=6;j<roiManager("Count");j++){
	//selectWindow("Rep1");
	run("Select None");
	roiManager("Select", j)
	List.setMeasurements();
	//toScaled(Results);
	ECround[j-6] = List.getValue("Round");
	X=List.getValue("X");//To get scaled values independently of resolution
	toScaled(x);
	centerX[j-6] = X;
	Y=List.getValue("Y");//To get scaled values independently of resolution
	toScaled(x);
	centerY[j-6] = Y;
	centerX[j-6] = Math.round(centerX[j-6]);
	centerY[j-6] = Math.round(centerY[j-6]);
}



pix_below_cutoff_R1=newArray(roiManager("Count")-6); //reports no. of pixels below Cutoff_R1
pix_below_cutoff_R2=newArray(roiManager("Count")-6); //reports no. of pixels below Cutoff_R2
pix_below_cutoff_M=newArray(roiManager("Count")-6); //reports no. of pixels below Cutoff_M

Per_below_cutoff_R1=newArray(roiManager("Count")-6); //reports %. of pixels below Cutoff_R1
Per_below_cutoff_R2=newArray(roiManager("Count")-6); //reports %. of pixels below Cutoff_R2
Per_below_cutoff_M=newArray(roiManager("Count")-6); //reports %. of pixels below Cutoff_M

//for Rep1
for(j=6;j<roiManager("Count");j++){
	selectWindow("Rep1");
	run("Select None");
	roiManager("Select", j);
	getHistogram(values, counts, nBins);
	pixel_num=0;
	for(i=0;i<lowerR1;i++){
		pix_below_cutoff_R1[j-6]=pix_below_cutoff_R1[j-6]+counts[i];
	}
	for(i=0;i<counts.length;i++){
		pixel_num=pixel_num+counts[i];
	}
	Per_below_cutoff_R1[j-6] = pix_below_cutoff_R1[j-6]/pixel_num*100;
	if (Per_below_cutoff_R1[j-6]<Cutoff_R1){ //Reporter1+
		Rep1[j-6]="Yes";
	}else{
		Rep1[j-6]="No";
	}
}

//for Rep2
for(j=6;j<roiManager("Count");j++){
	selectWindow("Rep2");
	run("Select None");
	roiManager("Select", j);
	getHistogram(values, counts, nBins);
	pixel_num=0;
	for(i=0;i<lowerR2;i++){
		pix_below_cutoff_R2[j-6]=pix_below_cutoff_R2[j-6]+counts[i];
	}
	for(i=0;i<counts.length;i++){
		pixel_num=pixel_num+counts[i];
	}
	Per_below_cutoff_R2[j-6] = pix_below_cutoff_R2[j-6]/pixel_num*100;
	if (Per_below_cutoff_R2[j-6]<Cutoff_R2){ //Reporter2+
		Rep2[j-6]="Yes";
	}else{
		Rep2[j-6]="No";
	}
}

//for Mark
for(j=6;j<roiManager("Count");j++){
	selectWindow("Mark");
	run("Select None");
	roiManager("Select", j);
	getHistogram(values, counts, nBins);
	pixel_num=0;
	for(i=0;i<lowerM;i++){
		pix_below_cutoff_M[j-6]=pix_below_cutoff_M[j-6]+counts[i];
	}
	for(i=0;i<counts.length;i++){
		pixel_num=pixel_num+counts[i];
	}
	Per_below_cutoff_M[j-6] = pix_below_cutoff_M[j-6]/pixel_num*100;
	if (Per_below_cutoff_M[j-6]<Cutoff_M){ //Reporter2+
		Mark[j-6]="Yes";
	}else{
		Mark[j-6]="No";
	}
}

//Decide color code of cell
rep1=0;
rep2=0;
rep12=0;
repNo=0;
rep1M=0;
rep2M=0;
rep12M=0;
repNoM=0;

for(j=6;j<roiManager("Count");j++){
	if((Rep1[j-6]=="No")&(Rep2[j-6]=="No")&(Mark[j-6]=="No")){ //Cell is rep1-,Rep2-,Marker-
				selectWindow("Temp");
				repNo+=1;
				run("Select None");
				roiManager("Select", j);
				roiManager("Set Color", "#848484");//Grey selection
				roiManager("Set Line Width", 2);
				run("Add Selection...");
			}else{
				if((Rep1[j-6]=="No")&(Rep2[j-6]=="No")&(Mark[j-6]=="Yes")){ //Cell is rep1-,Rep2-,Marker+
				selectWindow("Temp");
				repNoM+=1;
				run("Select None");
				roiManager("Select", j);
				roiManager("Set Color", "#0040FF"); //Blue selection
				roiManager("Set Line Width", 2);
				run("Add Selection...");
			}else{
				if((Rep1[j-6]=="Yes")&(Rep2[j-6]=="Yes")&(Mark[j-6]=="Yes")){ //Cell is rep1+,Rep2+,Marker+
				selectWindow("Temp");
				rep12M+=1;
				run("Select None");
				roiManager("Select", j);
				roiManager("Set Color", "white");
				roiManager("Set Line Width", 2);
				run("Add Selection...");
			}else{
				if((Rep1[j-6]=="Yes")&(Rep2[j-6]=="Yes")&(Mark[j-6]=="No")){ //Cell is rep1+,Rep2+,Marker-
					selectWindow("Temp");
					rep12+=1;
					run("Select None");
					roiManager("Select", j);
					roiManager("Set Color", "yellow");
					roiManager("Set Line Width", 2);
					run("Add Selection...");
				}else{
					if((Rep1[j-6]=="Yes")&(Rep2[j-6]=="No")&(Mark[j-6]=="No")){ //Cell is rep1+,Rep2-,Marker-
						selectWindow("Temp");
						rep1+=1;
						run("Select None");
						roiManager("Select", j);
						roiManager("Set Color", "red");
						roiManager("Set Line Width", 2);
						run("Add Selection...");
					}else{
						if((Rep1[j-6]=="Yes")&(Rep2[j-6]=="No")&(Mark[j-6]=="Yes")){ //Cell is rep1+,Rep2-,Marker+
							selectWindow("Temp");
							rep1M+=1;
							run("Select None");
							roiManager("Select", j);
							roiManager("Set Color", "magenta");
							roiManager("Set Line Width", 2);
							run("Add Selection...");
						}else{
							if((Rep1[j-6]=="No")&(Rep2[j-6]=="Yes")&(Mark[j-6]=="Yes")){ //Cell is rep1-,Rep2+,Marker+
								selectWindow("Temp");
								rep2M+=1;
								run("Select None");
								roiManager("Select", j);
								roiManager("Set Color", "cyan");
								roiManager("Set Line Width", 2);
								run("Add Selection...");
							}else{
								if((Rep1[j-6]=="No")&(Rep2[j-6]=="Yes")&(Mark[j-6]=="No")){ //Cell is rep1-,Rep2+,Marker+
									selectWindow("Temp");
									rep2+=1;
									run("Select None");
									roiManager("Select", j);
									roiManager("Set Color", "green");
									roiManager("Set Line Width", 2);
									run("Add Selection...");
								}
							}	
						}
					}
				}
			}
		}
	}		
}	

selectWindow("Temp");
saveAs(".tif", masks + name + "_Image Plus segmented objects");
close();





Table.create("All_Objects_Parameters");
Table.setColumn("Area in um2", ECarea);
Table.setColumn("Roundness", ECround);
Table.setColumn("Reporter1 mean", MeanRep1);
Table.setColumn("Reporter2 mean", MeanRep2);
Table.setColumn("Marker mean", MeanMark);
Table.setColumn("Reporter1 positivity", Rep1);
Table.setColumn("Reporter2 positivity", Rep2);
Table.setColumn("Marker positivity", Mark);
	
Table.save(masks + name + "_All_Objects_Parameters");


				
//SHOW results
totalrep1=0;
totalrep1=rep1+rep1M;
totalrep2=0;
totalrep2=rep2+rep2M;
totaldou=0;
totaldou=rep12+rep12M;	
totalneg=0;
totalneg=repNo+repNoM;	
setResult("Parameter", 0, "Mean value used in threshold for IsoB4");
setResult(name, 0, lowerE);
setResult("Parameter", 1, "Mean value used in threshold for Reporter1");
setResult(name, 1, lowerR1);
setResult("Parameter", 2, "Mean value used in threshold for Reporter2");
setResult(name, 2, lowerR2);
setResult("Parameter", 3, "Mean value used in threshold for Marker");
setResult(name, 3, lowerM);
setResult("Parameter", 4, "Histogram cut-off used for Reporter1");
setResult(name, 4, Cutoff_R1);
setResult("Parameter", 5, "Histogram cut-off used for Reporter2");
setResult(name, 5, Cutoff_R2);
setResult("Parameter", 6, "Histogram cut-off used for Marker");
setResult(name, 6, Cutoff_M);
setResult("Parameter", 7, "Reporter 1 + Reporter 2 -");
setResult(name, 7, totalrep1);
setResult("Parameter", 8, "Reporter 1 - Reporter 2 +");
setResult(name, 8, totalrep2);
setResult("Parameter", 9, "Reporter 1 + Reporter 2 +");
setResult(name, 9, totaldou);
setResult("Parameter", 10, "Negative");
setResult(name, 10, totalneg);
setResult("Parameter", 11, "Rep1+ Rep2+ Mark+");
setResult(name, 11, rep12M);
setResult("Parameter", 12, "Rep1+ Rep2+ Mark-");
setResult(name, 12, rep12);
setResult("Parameter", 13, "Rep1+ Rep2- Mark-");
setResult(name, 13, rep1);
setResult("Parameter", 14, "Rep1+ Rep2- Mark+");
setResult(name, 14, rep1M);
setResult("Parameter", 15, "Rep1- Rep2+ Mark-");
setResult(name, 15, rep2);
setResult("Parameter", 16, "Rep1- Rep2+ Mark+");
setResult(name, 16, rep2M);
setResult("Parameter", 17, "Rep1- Rep2- Mark-");
setResult(name, 17, repNo);		
setResult("Parameter", 18, "Rep1- Rep2- Mark+");
setResult(name, 18, repNoM);
setResult("Parameter", 19, "Total cell count");
setResult(name, 19, countCell);	


run("Close All");

	
updateResults();
saveAs("Results", masks + name + "_User_Def_Settings");

