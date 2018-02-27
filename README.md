
# NIMS\_to\_BIDS Conversion

## Summary

Takes neuroimaging downloaded from the Neurobiological Image Managament System at Stanford University (https://cni.stanford.edu/nims/) and converts it to BIDS (http://bids.neuroimaging.io/)

## Usage

python NIMS\_to\_BIDS.py NIMS\_dir BIDS\_dir

This version of NIMS\_to\_BIDS is intended to work through the command line on Sherlock. The script takes two arguments, the path to the NIMS_data and the path where the BIDS data should be saved.

This command requires that the sequences in NIMS_dir are named in accordance with BIDS standards to begin with.

For instance, in the protocol tasks should be named "task\_[task\_name]\_run\_[run\_num]". BIDS format 
requires that the name end up being "task-[task\_name]\_run-[run\_num]" (notice the underscores replaced by dashes). NIMS doesn't support such naming conventions. NIMS\_to\_BIDS.py corrects this.

NIMS normally uses a calibration scan at the beginning of any EPI scan. If a separate calibration scan is desired, it can be signaled by adding "sbref" to the end of the sequence name. Though these details won't affect any analyses going forward (NIMS has already made use of the calibration scan to reconstruct the EPI), NIMS_to_BIDS
will transfer such sequences appropriately.

BIDS is organized by subjects and sessions. To automate this process you must name your subjects
"[subject_ID]_[session_ID]", e.g. "s001\_1". This will be used to appropriately organize the data

```
NIMS Data Format (e.g)

|-- NIMS_data
    |-- 20150511_2026_16549
        |-- 16549_1_1_3Plane_Loc_SSFSE
        |-- 16549_2_1_HO_Shim
        |-- 16549_3_1_task_rest_run_1_sbref
        |-- 16549_4_1_task_rest_run_1_ssg
        |-- 16549_10_1_fmap_fieldmap
        |-- 16549_11_1_anat_T1w

BIDS Data Format (e.g) http://bids.neuroimaging.io/


|-- BIDS_data
    |-- task-rest_bold.json
    |-- dataset_description.json
    |-- sub-s001
        |-- ses-1
            |-- anat
                |-- sub-s001_T1w.nii.gz
            |-- func
                |-- sub-s001_task-rest_run-1_bold.nii.gz
```

## Dependencies

Pydeface must be installed. Follow instructions here: https://github.com/poldracklab/pydeface


