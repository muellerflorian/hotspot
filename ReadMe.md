# HotSpot
Allows the quantification of translation/transcription sites.

We developed this tool for cases, where an automated analysis is difficult, e.g. due to 
poor signal-to-noise ratio. With this tool, the user can manually annotate

1. translation/transcription sites (we refer to them as hotspots)
2. corresponding individual proteins/RNAs (referred to as spots).

The average signal of the individual spots is then calculated. By comparing integrated intensity of the hotpsots and individual spots, we approximate how many individual
spots are in these local aggregations. 

Tool requires to open an image pair (e.g. RNA and protein image) side-by-side, but only one of the images will be quantififed. This can be useful when analyzing protein data (e.g. for translation where the presence of RNA can be helpful to identify hotspots). However, you can also simply open the same image twice, e.g. if you want to analyze transcriptinal activity from smFISH or MS2

# System Requirements

## Hardware requirements
HotSpot requires only a standard computer.

## Software requirements

### OS Requirements
This package has been tested on *Win 10 Pro*.

### Matlab dependencies
RNA detection is performed with our prevously published Matlab package
<a href="https://bitbucket.org/muellerflorian/fish_quant" target="_blank">**FISH-quant.**</a>

FISH-quant requires the following **toolboxes**:
* Optimization toolbox
* Statistics toolbox
* Image processing toolbox
* (Optional) Parallel processing toolbox

FISH-quant has been tested on **Matlab 2019b**.

# Installation guide

## FISH-quant
You can obtain the latest version of FISH-quant
<a href="https://bitbucket.org/muellerflorian/fish_quant" target="_blank">**here.**</a>

Installation instruction are provided. Installation time is rapid and requires only
to download the most recent version.

## HotSpot
Download or clone this repository, and add the local folder and all subfolders to the Matlab path definition.

# Instructions for use: general remarks

## Starting HotSpot

Before starting HotSpot, make sure that the folder containing the source-code, and all subfolders are on the Matlab path definition.

Go to Matlab, and type `HotSpot` and press enter.

## Loading data

__IMPORTANT__ Protein and RNA images have to be in the same folder.

To load the protein image, press the button Load protein data. The RNA data will be loaded automatically.

To setup automated data loading, the first time you press this button, the program asks you for the **identifier for the protein and RNA data**. This is the unique string that distinguishes the two, in general you will specify the channel identifier that is used by the microscope to distinguish the images.

For instance, if the protein image is called `test_w1GPF_01.tif`, and the RNA image is called `test_w2CY3_01.tif`, then the unique identifiers are `w1GFP` and `w2CY3`. Once these identifiers are specified, you can select the protein image that you want to open, the RNA image will be loaded automatically.

## Navigating images

The two images are displayed side-by-side (protein left, RNA right). You can ZOOM and MOVE in the image after pressing the dedicated buttons. The two images are linked, e.g. the same zoom and move will be applied ot each image. You can also go through the z-stack with the dedicated controls, again the two images are linked. The contrast, can be changed for each image individually.

* Zoom-out is by shift-click.
* You can also enlarge the GUI by clicking on the corresponding controls.
* Sometimes changing the contrast undoes the zoom. So it is recommended to change the contrast first.

# Instructions for use: analysis modes

In the pulldown menu you can choose between different analysis modes

* Hotspots
* CoLocalization analysis

## HotSpots

### 2D vs 3D analysis

Spots (individual molecules, and the aggregated hotspots) are cropped in 3D with a user-defined crop-size. If the crop size is set to 0 for z, you can also analyse 2D images (e.g. MIPs). You can then also analyze movies, since in the results files the z-position is recorded (for a movie this will correspond to time).

### Quantification parameters

The cropping area is by default set to +/-3 pixel in xy, and +/-1 pixel in z. To change these values, press on the menu item SETTINGS.

### Definining transcriptional hotspots

You can define a hotspot in the protein image by clicking once on the corresponding control. When you move over the image, the cursor will change to a cross. Press on the spot that you want to highlight. It will then be displayed as a red circle, and added to the listbox. Already existing hotspots will be shown in yellow. Don&#39;t define the spots too close to the boundary (not within +/- 3 pixels in XY and +/-1 pixels in Z), otherwise the program will crash.

The position of the hotspot will be **automatically corrected slightly**. The program will look for the pixel with the maximum intensity within a region that has been specified for the analysis (usually +/- a few pixels in each direction (see below). These modified coordinates in XYZ will be saved.

You can also **delete** a hotspot, by selecting it from the listbox and pressing the corresponding button. You can also delete all hotspots with the corresponding button.

### Definining individual proteins

**IMPORTANT** : every hotspot has to have a list of proteins assigned

**IMPORTANT** : you have to terminate the selection of each list of protein as described above. Otherwise the protein list will not be saved, and the program will crash!

Once you have defined a hotspot, you can define the associated individual proteins by clicking on the corresponding button. Again, the cursor will transform in a cross when moving over the protein image. This cursor will stay active to select multiple proteins, each time you click another protein will be selected. To **terminate** the selection, click either enter when the cursor is a cross, or click outside the protein image. You can select more proteins for the same locus, by pressing on the Define proteins button again.

- It might be important to select proteins from different locations around the nucleus in order to get a somewhat homogeneous background when averaging them together.

- Don't define the spots to close to the boundary (not within +/- 3 pixels in XY and +/-1 pixels in Z), otherwise the program will crash.

- The precise location of a protein is modified as described for the hotspots.

- You can reuse the proteins from other hotspots. Two options are available
  - From the hotspot that you just defined before.
  - From another hotspots that you defined any moment. For this select in the lest the hotspot from which you want to use the proteins, then press the button Reuse proteins, and select in the listbox the hotspot to which you want to associate the proteins.

- You can delete either last defined protein or all proteins associated to the currently selcted hotspot with the corresponding buttons.

### Performing quantification

You simply have to press the button Analyze. For each translation spot, the averaged image of the specified proteins is calculated. Then the image of the protein and the hotspot are projected along Z (using the above mentioned cropping area for Z). The resulting image is fit with a 2D Gaussian. The integrated intensity below this Gaussian function and above background is calculated. The number of proteins at the hotspots is inferred by calculating the ratio between the integrated intensities.

The program also saves automatically a PNG image for each analyzed hotspot with a comparison of the projected image of the hotspot and the averaged protein image. These images can be found in subfolder `_HotSpot_`.

### Saving & Loading

#### Saving results

To save the results, choose Save results from the menu. The program proposes a default file-name, which is the name of the image followed by \_HOTspots\_RESULTS.txt. This file contains a detailed summary for each analyzed hotspot. Specifically, it contains its name, the nuymber of inferred proteins, and the fitting results for the hotspots and the averaged protein image.

#### Saving locations

To save the localization results, choose Save spots from the menu. The programm proposes a default file-name, which is the name of the image followed by \_HOTspots\_SPOTS.txt. This file contains a detailed summary for each analyzed hotspot. Specifically, it contains its name, the nuymber of inferred proteins, and the fitting results for the hotspots and the averaged protein image.

#### Loading locations

To load localization results, first laod the image, then choose Load spots from the menu. You can select a file containing the saved spot locations. The programm will then load the specified hotspots and associated proteins.

## CoLocalization

The analysis is performed on the MIP.

### Define cells

You have to define cells, before you can define proteins or RNAs. Cells are defined with a freehand form, simply draw the oultine with a pressed mouse button, and release when done.

You can delete a cell in case you are not satisfied.

### Define Proteins / RNAs

They are define in their respective images. Points marked as proteins are considere to be co-localized with RNA, while points marked as RNA are considered to NOT be co-localized and in the summary file.

You can delete either last marked protein/RNA or all proteins/RNAs associated to the currently selected cell.

### Saving results - summary

This will save a tab-delimited text file, where each cell is defined by one row listing the name of the image, the name of the cell, the number of defined proteins (localized spots), and the number of RNAs (not-localized spots).

### Saving results – spots

Here a text-delimited text file is saved specifying the XY coordinates of all elements, e.g. the outline of a cell with the associated proteins and RNAs.

# License
MIT License.