spire-src-cal
==================

Source Calibration code and example scripts for SPIRE

***This respository is for development purposes only***

For _official_ releases of code and products, please see official SPIRE websites and documents, such as:
 + [The SPIRE calibration Web](http://herschel.esac.esa.int/twiki/bin/view/Public/SpireCalibrationWeb)
 + [The SPIRE Handbook](http://herschel.esac.esa.int/Docs/SPIRE/spire_handbook.pdf)
 + [SPIRE Photometer Beam Analysis pages](http://herschel.esac.esa.int/twiki/bin/view/Public/SpirePhotometerBeamProfileAnalysis)
 + [The HIPE online help pages](http://herschel.esac.esa.int/hcss-doc-12.0/)

Description
-----------
Currently all scripts are in Jython for use in HIPE, but it is expected that
in the future scripts in additional languages will be added.

Summary of directory structure and scripts:
-------------------------------------------
 + README.md				[*This file*]
 + HIPE					[*HIPE scripts for users*]
    + SpireHandbookBundle.py [*HIPE scripts for point/extended calibration, released with HIPE 12*]
 + HIPE/dev				[*Scripts in development, not intended for public use*]
    + SpireHandbookBundle_dev.py [*Updated/cleaned version of point/extended calibration calculations*]
    + SemiExtendedBundle_dev.py [*Functions for partially-extended source calibrations*]
    + sources_dev.py [*Classes and functions for Source Profiles*]
    + beamfunctions.py [*Functions for dealing with beam profiles, including convolution with sources*]
 + CalScripts				[*Directory of stand-alone HIPE scripts for making calibration products (not intended for public use)*]
    + makeSCalRadialCorrBeam.py [*Creation of radial beam profiles (core, constant and normArea)*]
    + makeSCalPhotColorK.py [*Calculation of colour corrections (point and extended)*]
    + makeSCalPhotColorBeam.py [*Calculation of beam colour corrections*]
    + makeSCalPhotColorHfi.py [*Calculation of SPIRE-HFI cross-calibration product*]
    + makeSCalPhotColorAperture.py [*Calculate of aperture corrections*]
    + cal2OMtabs.py [*read calibration tables and output TeX tables (for OM)*]
    
