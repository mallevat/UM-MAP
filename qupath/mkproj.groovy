/*
 * UM-MAP, QuPath step: create a single-image QuPath project for one slide.
 *
 * SPDX-License-Identifier: MIT
 *
 * Usage:   QuPath script mkproj.groovy --args <slide>
 *
 * Input:   <slide>, a whole-slide image that QuPath can open (e.g. .svs).
 * Output:  a QuPath project in "<slide>-proj/" next to the slide.
 *
 * The project folder is placed next to the slide's canonical path (symbolic links resolved),
 * so the folder that holds the real file must be writable. The image name stored in the
 * project is the slide's file name; classify_export.groovy names its output table after it.
 * The image type set here (BRIGHTFIELD_H_DAB) is replaced by BRIGHTFIELD_H_E in
 * classify_export.groovy before any processing, so it does not affect the results.
 * qupath/run_slide.sh calls this script, then classify_export.groovy, then deletes the project.
 *
 * Dependencies: QuPath 0.4.3 (tested version; run from the command line).
 */

import groovy.io.FileType
import java.awt.image.BufferedImage
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.gui.commands.ProjectCommands
import java.io.File

if (args.size() == 0) {
    print "Please provide image name\n"
    return
}

File imagePath = new File(args[0])
print "Image path: $imagePath\n"
String canonicalImagePath = imagePath.getCanonicalPath()
print "Canonical path: $canonicalImagePath\n"
String fileName = imagePath.getName()
print("Filename: $fileName\n")

//Check if we already have a QuPath Project directory in there...
def projectName = canonicalImagePath + "-proj"
File directory = new File(projectName)

if (!directory.exists())
{
    print("No project directory, creating one!\n")
    directory.mkdirs()
}

// Create project
def project = Projects.createProject(directory , BufferedImage.class)

// Get serverBuilder
def support = ImageServerProvider.getPreferredUriImageSupport(BufferedImage.class, canonicalImagePath)
print "$support\n"
def builder = support.builders.get(0)

// Make sure we don't have null 
if (builder == null) {
    print "Image not supported: " + canonicalImagePath + "\n"
    return
}

// Add the image as entry to the project
print "Adding: " + canonicalImagePath + "\n"
entry = project.addImage(builder)

// Set a particular image type
def imageData = entry.readImageData()
imageData.setImageType(ImageData.ImageType.BRIGHTFIELD_H_DAB)
entry.saveImageData(imageData)

// Write a thumbnail if we can
var img = ProjectCommands.getThumbnailRGB(imageData.getServer());
entry.setThumbnail(img)

// Add an entry name (the filename)
entry.setImageName(fileName)

// Changes should now be reflected in the project directory
project.syncChanges()
