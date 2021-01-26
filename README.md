# breast_density
Source code to segment and measure breast parenchymal tissue in MRI based on Dixon method 
containing water-only and fat-only weighted images. 

### Source code
The project was initially written in Matlab 2011a. Some parts are in the process to be migrated to Python 3. 
To reproduce results from papers mentioned below, use the original Matlab source code. 

`/matlab` subfolder contains all Matlab code
`/python` subfolder contains the parts of the pipeline that have been migrated to Python 3 already. 

### External programs
The pipeline relies on external programs that need to be installed and be part of the environment path.
    * dcm2niix: [https://github.com/rordenlab/dcm2niix]. Tested version: 1.0.20190410
    * ANTs registration: [https://github.com/ANTsX/ANTs]. Tested version: 2.3.4
    * NiftyReg: [https://sourceforge.net/p/niftyreg/git/?source=navbar]. Tested versions: 20180328, v2.3.9
    * FSL: [https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/]. Tested versions: 5.0.11, 6.0.4

### Pipeline
The processing pipeline consists of three steps:
  1. Pre-processing: Converting dicom to nifti, reorient images, do bias-field correction and combine images
  2. breast segmentation: segment breast from background and other body parts using a template based registration approach. First the best matching template is determined by affine registration. Then the best matching template is non-rigidly registered to match the shape of subjects body and breasts. 
  3. segment parenchymal and fatty tissue in breast and compute summary measures. 

### Usage (Matlab)
Add all folders of the project to Matlab search path. Ensure that external programs can be called by Matlab `system` command by adding them to the system PATH. In case you get a `libstdc++` error during calling of external programs, it might help to start matlab with system libstdc++ instead of too old libstdc++ packaged with Matlab. Asuming that system `libstdc++` is in `/usr/lib64`:
    LD_PRELOAD=/usr/lib64/libstdc++.so.6 matlab -desktop

Water-only and fat-only dicom image files needs to be in a separate folders. 

The pipeline is started then as:
    pipeline (path_to_water_dicom, path_to_fat_dicom, path_to_working_folder, path_to_output_folder, path_to_templates)
    
In the working folder intermediate results are stored, such as volumes converted to Nifti format, bias-field correction, registration parameters. The folder is not deleted at the end of pipeline. The output folder contains the segmentation into fat and parenchymal tissue (`assignment.nii.gz`) and the density measures (`density.mat`). 

### Templates
The templates used in publication are not publicly available. 
To create a set of template breasts you need fat only and water only images and corresponding segmentation label of the breast (0 background, 1 foreground). In the following `f.nii.gz`, `w.nii.gz` and `s.nii.gz` are the fat only, water only and segmentation volumes in Nifti format. 
disp ('Calculate wf decision (combine w and f)');
    # First normalize images using ANTs ImageMath
    ImageMath 3 fn.nii.gz Normalize f.nii.gz
    ImageMath 3 wn.nii.gz Normalize w.nii.gz
    
    # Combine using ImageMath Decision function from ANTs
    ImageMath 3 wfdec.nii.gz Decision wn.nii.gz fn.nii.gz
    ImageMath 3 wfdec.nii.gz Normalize wfdec.nii.gz

Create for each template breast a subfolder containing `wfdec.nii.gz` and corresponding `s.nii.gz`. During segmentation of subject breast the program iterates over all sub-folders of given template folder and determines the best-matching template breast segmentation automatically. 

### Citation
Please use one of the following citations when using the program for research

    * Wengert, Georg Johannes; Helbich, Thomas H; Vogl, Wolf-Dieter; Baltzer, Pascal; Langs, Georg; Weber, Michael; Bogner, Wolfgang; Gruber, Stephan; Trattnig, Siegfried; Pinker, Katja; ",Introduction of an automated user–independent quantitative volumetric magnetic resonance imaging breast density measurement system using the Dixon sequence: Comparison with mammographic breast density assessment,Investigative radiology,50,2,73-80,2015
    * Wengert, Georg J; Pinker‐Domenig, Katja; Helbich, Thomas H; Vogl, Wolf‐Dieter; Clauser, Paola; Bickel, Hubert; Marino, Maria‐Adele; Magometschnigg, Heinrich F; Baltzer, Pascal A; ",Influence of fat–water separation and spatial resolution on automated volumetric MRI measurements of fibroglandular breast tissue,NMR in Biomedicine,29,6,702-708,2016
    * Wengert, GJ; Helbich, TH; Woitek, Ramona; Kapetas, P; Clauser, P; Baltzer, PA; Vogl, WD; Weber, M; Meyer-Baese, A; Pinker, Katja; ",Inter-and intra-observer agreement of BI-RADS-based subjective visual estimation of amount of fibroglandular breast tissue with magnetic resonance imaging: comparison to automated quantitative assessment,European radiology,26,11,3917-3922,2016

###### Bibtex
@article{wengert2015introduction,
  title={Introduction of an automated user--independent quantitative volumetric magnetic resonance imaging breast density measurement system using the Dixon sequence: Comparison with mammographic breast density assessment},
  author={Wengert, Georg Johannes and Helbich, Thomas H and Vogl, Wolf-Dieter and Baltzer, Pascal and Langs, Georg and Weber, Michael and Bogner, Wolfgang and Gruber, Stephan and Trattnig, Siegfried and Pinker, Katja},
  journal={Investigative radiology},
  volume={50},
  number={2},
  pages={73--80},
  year={2015},
  publisher={LWW}
}

@article{wengert2016influence,
  title={Influence of fat--water separation and spatial resolution on automated volumetric MRI measurements of fibroglandular breast tissue},
  author={Wengert, Georg J and Pinker-Domenig, Katja and Helbich, Thomas H and Vogl, Wolf-Dieter and Clauser, Paola and Bickel, Hubert and Marino, Maria-Adele and Magometschnigg, Heinrich F and Baltzer, Pascal A},
  journal={NMR in Biomedicine},
  volume={29},
  number={6},
  pages={702--708},
  year={2016}
}

@article{wengert2016inter,
  title={Inter-and intra-observer agreement of BI-RADS-based subjective visual estimation of amount of fibroglandular breast tissue with magnetic resonance imaging: comparison to automated quantitative assessment},
  author={Wengert, GJ and Helbich, TH and Woitek, Ramona and Kapetas, P and Clauser, P and Baltzer, PA and Vogl, WD and Weber, M and Meyer-Baese, A and Pinker, Katja},
  journal={European radiology},
  volume={26},
  number={11},
  pages={3917--3922},
  year={2016},
  publisher={Springer Berlin Heidelberg}
}
 
