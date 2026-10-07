/*
 * UM-MAP, QuPath step: cell detection, cell classification and per-cell export for one
 * H&E whole-slide image.
 *
 * SPDX-License-Identifier: MIT
 *
 * What it does, for the image in a single-image QuPath project (see mkproj.groovy):
 *   1. Image type H&E, with QuPath's "H&E default" stain vectors (hematoxylin 0.65111 0.70119
 *      0.29049; eosin 0.2159 0.8012 0.5581; background 255 255 255).
 *   2. Region to analyse:
 *        full                 one annotation covering the whole image (default), or
 *        <pixel classifier>   tissue annotations from a pixel classifier JSON (e.g. TissueROI.json),
 *                             minimum area 1.0E6 um^2, minimum hole area 10 um^2, split into
 *                             separate objects.
 *   3. Watershed cell detection on hematoxylin OD: pixel size 0.5 um, background radius 4 um
 *      (by reconstruction), median radius 0, sigma 1.5 um, nucleus area 10-400 um^2, maximum
 *      background 4, watershed post-processing, cell expansion 5 um. The detection threshold is
 *      an argument (default 0.05).
 *   4. Object classification of every cell with an object classifier JSON.
 *   5. Export of one tab-separated file per image, <output_dir>/<image name>.tsv, where every
 *      character of the image name other than A-Z, a-z and 0-9 is replaced by "_". Columns:
 *      Class, Centroid X µm, Centroid Y µm (one row per cell). This is the input of the
 *      G-cross step.
 *
 * Usage (QuPath --args are positional, in this order; trailing ones may be omitted, and an
 * empty value also means "use the default"):
 *
 *   QuPath script -p <project.qpproj> classify_export.groovy \
 *       --args <output_dir> --args <classifier> --args <tissue_mode> --args <threshold>
 *
 *   output_dir   folder for the TSV, created if needed                         (required)
 *   classifier   object classifier JSON, file name in the classifier folder    (Combo-12-08-24.json)
 *   tissue_mode  "full", or a pixel classifier JSON file name in the
 *                classifier folder, e.g. TissueROI.json                        (full)
 *   threshold    watershed detection threshold                                 (0.05)
 *
 *   The defaults are the UMICH2 settings. TCGA settings: classifier
 *   Tumor-Fibroblast-Lymphoid-Myeloid-2.json, tissue_mode TissueROI.json, threshold 0.1.
 *
 * Classifier folder: $UMMAP_CLASSIFIER_DIR if that environment variable is set, otherwise
 * the classifiers/ folder of this repository (found from this script's location).
 * qupath/run_slide.sh builds the project, runs this script and removes the project again.
 * Errors (for example a missing classifier file) stop the script with a message in the QuPath
 * log; QuPath 0.4.3 still exits with status 0 when a script run with -p fails, so check that
 * the table was written (run_slide.sh does this).
 *
 * Inputs: the image of the QuPath project given with -p, and classifier JSON files from the
 * classifier folder. The classifiers/ folder of this repository holds:
 *   Combo-12-08-24.json                      OpenCV ANN_MLP; classes Fibroblast, Garbage,
 *                                             Lymphoid, Myeloid, Tumor
 *   Tumor-Fibroblast-Lymphoid-Myeloid-2.json  OpenCV random trees; classes Fibroblast,
 *                                             Lymphoid, Myeloid, Tumor
 *   TissueROI.json                            pixel classifier for tissue detection
 *
 * Output: <output_dir>/<image name>.tsv, one row per cell (step 5).
 *
 * Dependencies: QuPath 0.4.3 (tested version; run from the command line).
 */

import static qupath.lib.gui.scripting.QPEx.*
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.opencv.ml.pixel.PixelClassifiers
import java.nio.file.FileSystems
import groovy.io.FileType
import java.awt.image.BufferedImage
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.gui.commands.ProjectCommands
import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathDetectionObject
import qupath.lib.objects.PathCellObject

// ---------------------------------------------------------------------------
// Arguments and settings
// ---------------------------------------------------------------------------
def DEFAULT_CLASSIFIER = "Combo-12-08-24.json"
def DEFAULT_TISSUE_MODE = "full"
def DEFAULT_THRESHOLD = "0.05"
def USAGE = "Usage: QuPath script -p <project.qpproj> classify_export.groovy --args <output_dir> " +
        "[--args <classifier>] [--args <tissue_mode: full | pixel classifier JSON>] [--args <threshold>]"

def argList = (binding.hasVariable('args') && args != null) ? (args as List).collect { it == null ? "" : it.trim() } : []
if (argList.size() < 1 || argList.size() > 4 || argList[0].isEmpty()) {
    throw new IllegalArgumentException(USAGE)
}
def argOrDefault = { int i, String fallback -> (argList.size() > i && !argList[i].isEmpty()) ? argList[i] : fallback }

def analysisFolder = argList[0]
def classifierFile = argOrDefault(1, DEFAULT_CLASSIFIER)
def tissueMode = argOrDefault(2, DEFAULT_TISSUE_MODE)
def thresholdArg = argOrDefault(3, DEFAULT_THRESHOLD)

// Watershed detection threshold. Its text is placed in the detection parameters below.
double threshold
try {
    threshold = Double.parseDouble(thresholdArg)
} catch (NumberFormatException e) {
    throw new IllegalArgumentException("Threshold must be a number, got '" + thresholdArg + "'. " + USAGE)
}
if (Double.isNaN(threshold) || Double.isInfinite(threshold) || threshold < 0) {
    throw new IllegalArgumentException("Threshold must be a finite number >= 0, got '" + thresholdArg + "'")
}
def thresholdText = String.valueOf(threshold)

// Folder holding the classifier JSON files: $UMMAP_CLASSIFIER_DIR, else <repository>/classifiers.
// QuPath 0.4.3 passes the path of the running script as the binding variable 'qupath.script.file'.
def configDir
def envClassifierDir = System.getenv('UMMAP_CLASSIFIER_DIR')
if (envClassifierDir != null && !envClassifierDir.trim().isEmpty()) {
    configDir = new File(envClassifierDir.trim()).getAbsolutePath()
} else if (binding.hasVariable('qupath.script.file') && binding.getVariable('qupath.script.file') != null) {
    def scriptFile = new File(binding.getVariable('qupath.script.file').toString()).getCanonicalFile()
    configDir = new File(scriptFile.getParentFile().getParentFile(), "classifiers").getAbsolutePath()
} else {
    throw new IllegalStateException("Cannot locate the classifier folder: set the environment variable UMMAP_CLASSIFIER_DIR")
}

// Classifiers are given as file names inside the classifier folder.
def checkJson = { String name, String what ->
    if (name.contains("/") || name.contains("\\")) {
        throw new IllegalArgumentException(what + " must be a file name inside the classifier folder (" + configDir +
                "), got '" + name + "'. To use another folder, set UMMAP_CLASSIFIER_DIR.")
    }
    if (!new File(configDir, name).isFile()) {
        def available = new File(configDir).list()?.findAll { it.toLowerCase().endsWith(".json") }?.sort()
        throw new IllegalArgumentException(what + " not found: " + new File(configDir, name) + ". JSON files in the classifier folder: " + available)
    }
}
checkJson(classifierFile, "Object classifier")
def usePixelClassifier = !tissueMode.equalsIgnoreCase("full")
def roiFile = tissueMode
if (usePixelClassifier) {
    checkJson(roiFile, "Tissue pixel classifier")
}

def qupathVersion = "unknown"
try {
    qupathVersion = String.valueOf(qupath.lib.common.GeneralTools.getVersion())
} catch (Throwable ignored) {
}
print "QuPath version: " + qupathVersion + " (UM-MAP was validated with 0.4.3)\n"
print "Analysis folder: " + analysisFolder + "\n"
print "Classifier folder: " + configDir + "\n"
print "Classifier file: " + classifierFile + "\n"
print "Tissue region: " + (usePixelClassifier ? "pixel classifier " + roiFile : "full image") + "\n"
print "Watershed threshold: " + thresholdText + "\n"

// ---------------------------------------------------------------------------
// Processing
// ---------------------------------------------------------------------------
// Load the project which is specified as an argument to the script.
def project = getProject()
def imageList = project.getImageList()

print "Image: " + imageList[0] + "\n"
def imageData = getCurrentImageData()

setImageType('BRIGHTFIELD_H_E');
setColorDeconvolutionStains('{"Name" : "H&E default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049", "Stain 2" : "Eosin", "Values 2" : "0.2159 0.8012 0.5581", "Background" : " 255 255 255"}');

if (!usePixelClassifier) {
    // Tissue region = the whole image
    // Clear any existing selections
    resetSelection()

    // Create a full image annotation instead of using pixel classifier
    print "Creating full image annotation...\n"
    createFullImageAnnotation(true)

    // Select the annotation we just created
    selectAnnotations();
} else {
    // Tissue region = annotations from a pixel classifier
    // Load the classifier from the file
    def classifierPath = FileSystems.getDefault().getPath(configDir + "/" + roiFile)
    print "Classifier Path: " + classifierPath + "\n"
    def classifier = PixelClassifiers.readClassifier(classifierPath)

    // Check if classifier is loaded
    if (classifier == null) {
        print 'Error loading classifier!'
        return
    }

    print "Classifier: $classifier\n"

    createAnnotationsFromPixelClassifier(classifier, 1.0E6, 10.0, "SPLIT", "DELETE_EXISTING", "SELECT_NEW")
}

// Get the annotations that were just created
def annotations = getAnnotationObjects()
print "Number of annotations created: " + annotations.size() + "\n"

// Perform cell detection using selected ROI
print "Running cell detection...\n"
def cellDetectionParams = '{"detectionImageBrightfield":"Hematoxylin OD","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":4.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.5,"minAreaMicrons":10.0,"maxAreaMicrons":400.0,"threshold":' + thresholdText + ',"maxBackground":4.0,"watershedPostProcess":true,"cellExpansionMicrons":5.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}'
print "Cell detection parameters: " + cellDetectionParams + "\n"
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', cellDetectionParams)

// Run object classification
print "Running object classification...\n"
runObjectClassifier(configDir + "/" + classifierFile)

def entry_proj = getProjectEntry()
print "Proj Entry: " + entry_proj + "\n"
entry_proj.saveImageData(imageData)

// Ensure the directory exists (create if it doesn't)
new File(analysisFolder).mkdirs()

// Separate each measurement value in the output file with a tab ("\t")
def separator = "\t"

// Choose the columns that will be included in the export
def columnsToInclude = new String[]{"Class", "Centroid X µm", "Centroid Y µm"}

// Choose the type of objects that the export will process
def exportType = PathCellObject.class

def exporter = new MeasurementExporter()
.separator(separator) // Character that separates values
.includeOnlyColumns(columnsToInclude) // Columns are case-sensitive
.exportType(exportType) // Type of objects to export

// Iterate through each image and export measurements
def imagesToExport = project.getImageList()
imagesToExport.each { image ->
    def imageName = image.getImageName().replaceAll("[^a-zA-Z0-9]", "_")
    def outputPath = "${analysisFolder}/${imageName}.tsv"
    def outputFile = new File(outputPath)

    print "Exporting to: " + outputPath + "\n"

    exporter.imageList([image])             // Set current image
            .exportMeasurements(outputFile) // Start the export process for this image

    print "Export completed for: " + imageName + "\n"
}
