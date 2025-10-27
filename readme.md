READ ME


```
## Advanced Material & Object Aligner (AMOA)

## Overview

The **Advanced Material & Object Aligner (AMOA)** is a Blender add-on designed to streamline your 3D workflow by automating two key processes: material transfer and object alignment. This tool makes it easy to match materials and align objects across different collections, saving you significant time and effort, especially when working with complex scenes or imported models.


## Features

### Material Transfer

* **Automatic Matching**: Finds the best source object for each target object based on multiple criteria like volume, surface area, and mesh topology.
* **Per-Face & Layout-Based Transfer**: Accurately transfers multi-material assignments, even between meshes with different face counts, by using either a layout-based approach or an advanced spatial grid for complex meshes.
* **Intelligent Fallbacks**: Employs multiple strategies (including z-banding and object-level assignment) to ensure a successful material transfer, even when a perfect match isn't found.

### Object Alignment

* **Collection Alignment**: Aligns an entire target collection to a source collection, maintaining the relative positions of objects within the target collection.
* **Individual Mesh Alignment**: A dedicated option allows you to align each matched target mesh to its corresponding source mesh individually.
* **Advanced Alignment Algorithms**: Uses sophisticated methods like **Principal Component Analysis (PCA)** and **Iterative Closest Point (ICP)** to achieve precise and robust alignment.
* **RANSAC Filtering**: Employs a **Random Sample Consensus (RANSAC)** algorithm to filter out mismatched pairs and find stable alignments.

## How to Use
 

1.  **Installation**: Install the add-on in Blender via `Edit > Preferences > Add-ons > Install...`.
2.  **Enable**: Search for "AMOA" in the add-ons list and enable it.
3.  **Access**: The AMOA panel can be found in the 3D Viewport sidebar (`N` key) under its dedicated tab.
4.  **Select Collections**:
    * Choose your **Source Collection** (the collection with the materials and desired positions).
    * Choose your **Target Collection** (the objects you want to align and transfer materials to).
5.  **Adjust Settings**: (Found in the "Settings & Extras" dropdown)
    * **Align Individual Meshes**: Toggle this to align each matched mesh separately, rather than aligning the entire collection as one unit.
    * **Performance Mode (Destructive)**: Speeds up operations by modifying mesh data directly. **Warning**: This will affect all linked duplicates (Alt+D objects) of your target meshes.
    * **Animate Alignment**: Automatically creates keyframes to animate the alignment from the original to the final position.
6.  **Run Operations**:
    * **Transfer Materials**: Click this button to apply materials from the source collection to the target collection.
    * **Align Objects**: Click this button to align the target collection to the source collection.




## Changelog

v1.11.0

* Added a toggle for **individual mesh alignment**.
* Improved object matching with a multi-pass approach, including RANSAC filtering.
* Introduced an advanced spatial grid method for multi-material transfer, which is now applied **automatically** for complex meshes rather than via a manual toggle.
* Added fallback alignment options (e.g., multi-point SVD, single-point PCA, and global centroid alignment) to handle challenging cases.

## License:

To keep the addon's filesize low i have moved the c++-, build-, library-, and source-files to:   
https://github.com/HesliGIT/AMOA_Addon_build/


## Disclaimer:
 
High amount of (complex) objects within collections will take a lot of time to compute. Please note that blender may freeze for a few minutes, while using this addon on complex objects and high volume of objects, depending on your system.

## About Me: 

I'm Hesli Reiling, a 3D developer & designer focused on making our work in architectural visualization better. With a background in interactive UX, I design addons to solve the real-world challenges we all face, helping you optimize your workflow and free up more time for creativity.

Thank you for trying out my addon! I hope it will save you a lot of time during your projects. Feel free to send me any questions or report problems you are experiencing with the Advanced Material & Object Aligner addon.

```
