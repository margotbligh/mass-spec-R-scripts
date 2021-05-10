# Scripts for analysis of DI-MS, GC-MS and LC-MS data recorded for studies of marine sugars.
Working space of Margi and Hagi. Scripts are categorised as either:
*  **clean**: 
    *  script has no errors from start to finish
    *  no "extra" code for testing
    *  everything commented nicely
    *  **code should be able to work with any dataset**. either comments make clear what needs to be changed, or all variables are controlled by data.
*  **working**: 
    *  messy code
    *  lots of test code
    *  code can be specific for the dataset the script is written for
    *  **should not be downloaded by others**

**Table of contents**
*  [LC-MS](#LC-MS-(liquid-chromatography-mass-spectrometry)-scripts)
    * [Clean scripts](#LC-MS-clean-scripts)
    * [Working scripts](#LC-MS-working-scripts) 
*  [DI-MS](#DI-MS-(direct-infusion-mass-spectrometry)-scripts)
    * [Clean scripts](#DI-MS-clean-scripts)
    * [Working scripts](#DI-MS-working-scripts)
*  [GC-MS](#GC-MS-(gas-chromatography-mass-spectrometry)-scripts)
    * [Clean scripts](#GC-MS-clean-scripts)
    * [Working scripts](#GC-MS-working-scripts)    

## LC-MS (liquid chromatography mass spectrometry) scripts
### LC-MS clean scripts

None yet :upside_down_face:

### LC-MS working scripts

*  [_**workflow-test_analysis.R**_](./LC-MS/working-scripts/workflow-test_analysis.R)

Data: standards and blanks for workflow including procainamide derivatisation. Final output: extracted ion chromatograms for identified features and fluorescence chromatogram. Margot only


*  [_**hagi-extractions-proca_20210503.R**_](./LC-MS/working-scripts/hagi-extractions-proca_20210503.R)

Data: standards, blanks and samples from POS531, both with and without procainamide derivatisation. Output: volcano plots, heatmaps, PCA plots, extracted ion chromatograms. Very messy right now. MS/MS extracted, but collision was at 200... Margot only


*  [_**glycan diversity for TMC.R**_](./LC-MS/working-scripts/glycan%20diversity%20for%20TMC.R)

Hagi only. Need description

*  [_**fitdog_analysis_v3.R**_](./LC-MS/working-scripts/fitdog_analysis_v3.R)

Data: glycans digested by FITDOG and enzymatic digestion. Final output: extracted ion chromatograms for identified features and fluorescence chromatogram, dot plot summarising peak detection, MS/MS plot comparison. Margot only

*  [_**a-mannan_v5.R**_](./LC-MS/working-scripts/a-mannan_v5.R)

Data: a-mannan digested with GH99. Final output: extracted ion chromatograms for identified features and fluorescence chromatogram, MS/MS plot. Margot only

*  [_**Analysis and visualization of LC-IM-MS data Nantes Oct 2019 for pitch.R**_](./LC-MS/working-scripts/Analysis%20and%20visualization%20of%20LC-IM-MS%20data%20Nantes%20Oct%202019%20for%20pitch.R)

Hagi only. Need description


## DI-MS (direct infusion mass spectrometry) scripts
### DI-MS clean scripts

*  [**_DI-MS_average_spectra.R_**](./DI-MS/clean-scripts/DI-MS_average_spectra.R)

Average spectra by file within a specific time range. Assumes only MS1 data. For the sample set the script was written for, an air bubble was in the capillary while one sample was infused. This is confirmed by plotting the TIC for that sample. Margot only

### DI-MS working scripts

*  [_**GCC-SPE_yield_analysis_final.R**_](./DI-MS/clean-scripts/GCC-SPE_yield_analysis_final.R)

Peaks were picked for averaged spectra output by DI-MS_average_spectra.R in mMass. This script works with the peak lists output by mMass (.txt format). It filters for ions of interest, then performs analyses based on linear regression models. Margot only

*  [_**FTICRMS_LnSt_GCCSPE_allsamples.R**_](./DI-MS/clean-scripts/FTICRMS_LnSt_GCCSPE_allsamples.R)

For analysis of FT-ICR-MS data. Samples were mixes of small polar metabolite standards at different concentrations with different extraction methods. Peaks picked with mass spec wavelet algorithm, grouping with clustering algorithm, and peaks filled (using XCMS/MSnBase). Isotopes and adducts detected with CAMERA - treat all spectra as one pseudospectra (no retention time dimension to separate on). Collapse isotopes to only retain M isotope. Annotate with predicted oligosaccharides (predictions made with [my tool](https://github.com/margotbligh/sugarMassesPredict)). Get per sample counts of standards and plot representation of counts. Margot only

*  [_**FITDOG_neutral-standards_MS-direct-injection.R**_](./DI-MS/clean-scripts/FITDOG_neutral-standards_MS-direct-injection.R)

No data processing. Just extracts and plots (base R) chromatograms for hexose saccharides with DP2-13 in one loop. Margot only

## GC-MS (gas chromatography mass spectrometry) scripts
### GC-MS clean scripts

None yet :upside_down_face:

### GC-MS working scripts

*  [_**GCMS-data-analysis-script_ASSEMBLE.r**_](./DI-MS/clean-scripts/GCMS-data-analysis-script_ASSEMBLE.r)

GC-MS data for Finland samples. Peak picking, grouping, RT alignment etc done with extensive parameter testing and checking by Hagi / Mona. Now trying to do deconvolution and database matching. Hagi and Margot
