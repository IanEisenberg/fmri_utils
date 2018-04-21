import argparse
import glob
import json
import nibabel as nib
import nipype.interfaces.fsl as fsl
import numpy as np
import os
import re
import shutil
import subprocess
import sys
import tarfile

# *****************************
# *** Constant Variables
# *****************************

cardiac_bids = {
    "StartTime": -30.0,
    "SamplingFrequency": .01,
    "Columns": ["cardiac"]
}

respiratory_bids = {
    "StartTime": -30.0,
    "SamplingFrequency":  .04,
    "Columns": ["respiratory"]
}


# *****************************
# *** HELPER FUNCTIONS
# *****************************
def clean_file(f):
    f = f.replace('task_', 'task-').replace('run_','run-').replace('_ssg','')
    return f

def cleanup(path):
    for f in glob.glob(os.path.join(path, '*')):
        new_name = clean_file(f)
        os.rename(f,new_name)

def bids_anat(sub_id, anat_dir, anat_path):
    """
    Moves and converts a anatomical folder associated with a T1 or T2
    to BIDS format. Assumes files are named appropriately
    (e.g. "[examinfo]_anat_[T1w,T2w]")
    """
    #T1 or T2?
    anat_type = anat_dir[(anat_dir.index('anat_')+5):]
    anat_file = glob.glob(os.path.join(anat_dir,'*nii.gz'))
    assert len(anat_file) == 1, "More than one anat file found in directory %s" % anat_dir
    new_file = os.path.join(anat_path, sub_id + '_' + anat_type + '.nii.gz')
    if not os.path.exists(new_file):
        # deface
        print('\tDefacing...')
        subprocess.call("pydeface %s --outfile %s" % (anat_file[0], new_file), shell=True)
    else:
        print('Did not save anat because %s already exists!' % new_file)

def bids_fmap(sub_id, fmap_dir, fmap_path, func_path):
    """
    Moves and converts an epi folder associated with a fieldmap
    to BIDS format. Assumes files are named appropriately
    (e.g. "[examinfo]_fmap_fieldmap")
    """
    fmap_files = glob.glob(os.path.join(fmap_dir,'*.nii.gz'))
    assert len(fmap_files) == 2, "Didn't find the correct number of files in %s" % fmap_dir
    fmap_index = [0,1]['fieldmap.nii.gz' in fmap_files[1]]
    fmap_file = fmap_files[fmap_index]
    mag_file = fmap_files[1-fmap_index]
    fmap_name = sub_id + '_' + 'fieldmap.nii.gz'
    mag_name = sub_id + '_' + 'magnitude.nii.gz'
    if not os.path.exists(os.path.join(fmap_path, fmap_name)):
        shutil.copyfile(fmap_file, os.path.join(fmap_path, fmap_name))
        shutil.copyfile(mag_file, os.path.join(fmap_path, mag_name))
        func_runs = [os.sep.join(os.path.normpath(f).split(os.sep)[-3:]) for f in glob.glob(os.path.join(func_path,'*task*bold.nii.gz'))]
        fieldmap_meta = {'Units': 'Hz', 'IntendedFor': func_runs}
        json.dump(fieldmap_meta,open(os.path.join(fmap_path, sub_id + '_fieldmap.json'),'w'))
    else:
        print('Did not save fmap_epi because %s already exists!' % os.path.join(fmap_path, fmap_name))

def bids_sbref(sub_id, sbref_dir, func_path, bids_dir):
    """
    Moves and converts an epi folder associated with a sbref
    calibration scan to BIDS format. Assumes tasks are named appropriately
    (e.g. "[examinfo]_task_[task]_run_[run_number]_sbref")
    """
    task_index = sbref_dir.index('task')
    filename = sbref_dir[task_index:]
    sbref_files = glob.glob(os.path.join(sbref_dir,'*.nii.gz'))
    # remove files that are sometimes added, but are of no interest
    sbref_files = [i for i in sbref_files if 'phase' not in i]
    assert len(sbref_files) <= 1, "More than one func file found in directory %s" % sbref_dir
    if len(sbref_files) == 0:
        print('Skipping %s, no nii.gz file found' % sbref_dir)
        return

    # bring to subject directory and divide into sbref and bold
    sbref_file = clean_file(os.path.join(func_path, sub_id + '_' + filename + '.nii.gz'))
    # check if file exists. If it does, check if the saved file has more time points
    if os.path.exists(sbref_file):
        print('%s already exists!' % sbref_file)
        saved_shape = nib.load(sbref_file).shape
        current_shape = nib.load(sbref_files[0]).shape
        print('Dimensions of saved image: %s' % list(saved_shape))
        print('Dimensions of current image: %s' % list(current_shape))
        if (current_shape[-1] <= saved_shape[-1]):
            print('Current image has fewer or equivalent time points than saved image. Exiting...')
            return
        else:
            print('Current image has more time points than saved image. Overwriting...')
    # save sbref image to bids directory
    shutil.copyfile(sbref_files[0], sbref_file)
    # get metadata
    sbref_meta_path = clean_file(os.path.join(bids_dir, re.sub('_run[-_][0-9]','',filename) + '.json'))
    if not os.path.exists(sbref_meta_path):
        try:
            meta_file = [x for x in glob.glob(os.path.join(sbref_dir,'*.json')) 
                            if 'qa' not in x][0]
            func_meta = get_functional_meta(meta_file, filename)
            json.dump(func_meta,open(sbref_meta_path,'w'))
        except IndexError:
            print("Metadata couldn't be created for %s" % sbref_file)


def bids_task(sub_id, task_dir, func_path, bids_dir):
    """
    Moves and converts an epi folder associated with a task
    to BIDS format. Assumes tasks are named appropriately
    (e.g. "[examinfo]_task_[task]_run_[run_number]_ssg")
    """
    task_index = task_dir.index('task')
    taskname = task_dir[task_index:]
    task_file = [f for f in glob.glob(os.path.join(task_dir,'*.nii.gz')) if "fieldmap" not in f]
    assert len(task_file) <= 1, "More than one func file found in directory %s" % task_dir
    if len(task_file) == 0:
        print('Skipping %s, no nii.gz file found' % task_dir)
        return

    # bring to subject directory and divide into sbref and bold
    bold_file = os.path.join(func_path, sub_id + '_' + taskname + '.nii.gz')
    bold_file = bold_file.replace('_ssg', '_bold')
    # check if file exists. If it does, check if the saved file has more time points
    if os.path.exists(bold_file):
        print('%s already exists!' % bold_file)
        saved_shape = nib.load(bold_file).shape
        current_shape = nib.load(task_file[0]).shape
        print('Dimensions of saved image: %s' % list(saved_shape))
        print('Dimensions of current image: %s' % list(current_shape))
        if (current_shape[-1] <= saved_shape[-1]):
            print('Current image has fewer or equal time points than saved image. Exiting...')
            return
        else:
            print('Current image has more time points than saved image. Overwriting...')
    # save bold image to bids directory
    shutil.copyfile(task_file[0], bold_file)
    # get epi metadata
    bold_meta_path = os.path.join(bids_dir, re.sub('_run[-_][0-9]','',taskname) + '_bold.json')
    bold_meta_path = clean_file(bold_meta_path)
    if not os.path.exists(bold_meta_path):
        meta_file = [x for x in glob.glob(os.path.join(task_dir,'*.json')) if 'qa' not in x][0]
        func_meta = get_functional_meta(meta_file, taskname)
        json.dump(func_meta,open(bold_meta_path,'w'))
    # get physio if it exists
    physio_file = glob.glob(os.path.join(task_dir, '*physio.tgz'))
    if len(physio_file)>0:
        assert len(physio_file)==1, ("More than one physio file found in directory %s" % task_dir)
        tar = tarfile.open(physio_file[0])
        tar.extractall(func_path)
        # extract the actual filename of the physio data
        physio_file = os.path.basename(physio_file[0])[:-4]
        for pfile in glob.iglob(os.path.join(func_path, physio_file, '*Data*')):
            pname = 'respiratory' if 'RESP' in pfile else 'cardiac'
            new_physio_file = bold_file.replace('_bold.nii.gz', 
                                '_recording-' + pname + '_physio.tsv.gz')
            f = np.loadtxt(pfile)
            np.savetxt(new_physio_file, f, delimiter = '\t')
        shutil.rmtree(os.path.join(func_path, physio_file))

def get_functional_meta(json_file, taskname):
    """
    Returns BIDS meta data for bold 
    """
    meta_file = json.load(open(json_file,'r'))
    meta_data = {}
    mux = meta_file['num_bands']
    nslices = meta_file['num_slices'] * mux
    tr = meta_file['tr']
    n_echoes = meta_file['acquisition_matrix_y'] 

    # fill in metadata
    meta_data['TaskName'] = taskname.split('_')[1]
    meta_data['EffectiveEchoSpacing'] = meta_file['effective_echo_spacing']
    meta_data['EchoTime'] = meta_file['te']
    meta_data['FlipAngle'] = meta_file['flip_angle']
    meta_data['RepetitionTime'] = tr
    # slice timing
    meta_data['SliceTiming'] = get_slice_timing(nslices, tr, mux = mux)
    total_time = (n_echoes-1)*meta_data['EffectiveEchoSpacing']
    meta_data['TotalReadoutTime'] = total_time
    meta_data['PhaseEncodingDirection'] = ['i','j','k'][meta_file['phase_encode']] + '-'        
    return meta_data

def get_fmap_epi_meta(json_file, intended_list):
    """
    Returns BIDS meta data for epi fieldmaps
    """
    meta_file = json.load(open(json_file,'r'))
    meta_data = {}
    mux = meta_file['num_bands']
    nslices = meta_file['num_slices'] * mux
    n_echoes = meta_file['acquisition_matrix_y'] 
    # fill in metadata
    meta_data['EffectiveEchoSpacing'] = meta_file['effective_echo_spacing']
    total_time = (n_echoes-1)*meta_data['EffectiveEchoSpacing']
    meta_data['TotalReadoutTime'] = total_time
    meta_data['PhaseEncodingDirection'] = ['i','j','k'][meta_file['phase_encode']]    
    meta_data['IntendedFor'] = intended_list        
    return meta_data

def get_slice_timing(nslices, tr, mux = None, order = 'ascending'):
    """
    nslices: int, total number of slices
    tr: float, repetition total_time
    mux: int, optional mux factor
    """
    if mux:
        assert nslices%mux == 0
        nslices = nslices//mux
        mux_slice_acq_order = list(range(0,nslices,2)) + list(range(1,nslices,2))
        mux_slice_acq_time = [float(s)/nslices*tr for s in range(nslices)]
        unmux_slice_acq_order = [nslices*m+s for m in range(mux) for s in mux_slice_acq_order]
        unmux_slice_acq_time = mux_slice_acq_time * mux
        slicetimes = list(zip(unmux_slice_acq_time,unmux_slice_acq_order))
    else:
        slice_acq_order = list(range(0,nslices,2)) + list(range(1,nslices,2))
        slice_acq_time = [float(s)/nslices*tr for s in range(nslices)]
        slicetimes = list(zip(slice_acq_time,slice_acq_order))
    #reorder slicetimes by slice number
    sort_index = sorted(enumerate([x[1] for x in slicetimes]), key= lambda x: x[1])
    sorted_slicetimes = [slicetimes[i[0]][0] for i in sort_index]
    return sorted_slicetimes

def get_subj_path(nims_path, bids_dir, id_correction_dict=None):
    """
    Takes a nims path and returns a subject id
    If a dictionary specifying id corrections is provided
    (in the form of a json file with exam numbers as keys and
    ids as values), the function will return the corrected id number
    """
    exam_number = nims_path.split('_')[-1]
    try:
        meta_json = glob.glob(os.path.join(nims_path,'*','*1.json'))[0]
        meta_file = json.load(open(meta_json,'r'))
        json_exam_number = str(meta_file['exam_number'])
        sub_session = str(meta_file['patient_id'].split('@')[0])
        assert json_exam_number == exam_number
    except IndexError:
        print('No meta json found for %s' % nims_path)
        sub_session = None
    
    # correct session if provided
    if id_correction_dict:
        sub_session = id_correction_dict.get(exam_number,sub_session)
    if sub_session is not None:
        sub_id = sub_session.split('_')[0]
        session = '1'
        if '_' in sub_session:
            session = sub_session.split('_')[1]
        subj_path = os.path.join(bids_dir, 'sub-'+sub_id, 'ses-'+session)
        return subj_path
    else:
        return None

def mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
        print('Directory %s already existed' % path)
    return path

# *****************************
# *** Main BIDS function
# *****************************
def bids_subj(subj_path, bids_dir, nims_path):
    """
    Takes a subject path (the BIDS path to the subject directory),
    a data path (the path to the BIDS directory), and a 
    nims_path (the path to the subject's data in the original format) 
    and moves/converts that subject's data to BIDS
    """
    if os.path.exists(subj_path):
        print("Path %s already exists. Skipping." % subj_path)
    else:
        print("********************************************")
        print("BIDSifying %s" % subj_path)
        print("Using nims path: %s" % nims_path)
        print("********************************************")
        # extract subject ID
        split_path = os.path.normpath(subj_path).split(os.sep)
        sub_id = [x for x in split_path if 'sub' in x][0]
        ses_id = [x for x in split_path if 'ses' in x][0]
        base_file_id = sub_id + '_' + ses_id
        # split subject path into a super subject path and a session path
        session_path = subj_path
        subj_path = os.path.split(subj_path)[0]
        mkdir(subj_path)
        mkdir(session_path)
        anat_path = mkdir(os.path.join(session_path,'anat'))
        func_path = mkdir(os.path.join(session_path,'func'))
        fmap_path = mkdir(os.path.join(session_path,'fmap'))
        # strip paths for rsync transfer
        stripped_anat_path = os.sep.join(anat_path.split(os.sep)[-3:-1])
        stripped_func_path = os.sep.join(func_path.split(os.sep)[-3:-1])
        stripped_fmap_path = os.sep.join(fmap_path.split(os.sep)[-3:-1])

        print(anat_path)
        print('BIDSifying anatomy...')
        # anat files
        anat_dirs = sorted(glob.glob(os.path.join(nims_path,'*anat*')))[::-1]
        for anat_dir in anat_dirs:
            print('\t' + anat_dir)
            bids_anat(base_file_id, anat_dir, anat_path)

        print('BIDSifying sbref...')
        # task files
        sbref_dirs = sorted(glob.glob(os.path.join(nims_path,'*sbref*')))[::-1]
        for sbref_dir in sbref_dirs:
            print('\t' + sbref_dir)
            bids_sbref(base_file_id, sbref_dir, func_path, bids_dir)

        print('BIDSifying task...')
        # task files
        task_dirs = sorted(glob.glob(os.path.join(nims_path,'*task*')))[::-1]
        for task_dir in [x for x in task_dirs if 'sbref' not in x]:
            print('\t' + task_dir)
            bids_task(base_file_id, task_dir, func_path, bids_dir)

        # cleanup
        cleanup(func_path)

        print('BIDSifying fmap...')
        # fmap files
        fmap_dirs = sorted(glob.glob(os.path.join(nims_path,'*fieldmap*')))[::-1]
        for fmap_dir in fmap_dirs:
            print('\t' + fmap_dir)
            bids_fmap(base_file_id, fmap_dir, fmap_path, func_path)


        
        
# *****************************
# *** ORGANIZE IN BIDS
# *****************************

# parse arguments
parser = argparse.ArgumentParser(description='fMRI Analysis Entrypoint Script.')

parser.add_argument('nims_dir', help='Directory of the non-BIDS fmri data')
parser.add_argument('bids_dir', help='Directory of the BIDS fmri data')
parser.add_argument('--study_id', default=None, help='Study ID. If not supplied, the directory above the nims_dir will be used')
parser.add_argument('--id_correction', help='JSON file that lists subject id corrections for fmri scan IDs')
parser.add_argument('--nims_paths', nargs='+', default=None)
args, unknown = parser.parse_known_args()

# directory with bids data
nims_dir = args.nims_dir
# bids directory
bids_dir = args.bids_dir
#study name
if args.study_id:
    study_id == args.study_id
else:
    study_id = nims_dir.strip(os.sep).split(os.sep)[-2]
# set id_correction_dict if provided
id_correction_dict = None
if args.id_correction:
    id_correction_dict = json.load(open(args.id_correction,'r'))
print('Using ID correction json file: %s' % args.id_correction)
#header file
header = {'Name': study_id, 'BIDSVersion': '1.51-rc1'}
json.dump(header,open(os.path.join(bids_dir, 'dataset_description.json'),'w'))
# error file
error_file = os.path.join(bids_dir, 'error_record.txt')

# bidsify all subjects in path
if args.nims_paths is None:
    nims_paths = glob.glob(os.path.join(nims_dir, '*'))
else:
    nims_paths = [glob.glob(os.path.join(nims_dir, '*'+path))[0] for path in args.nims_paths]
for i, nims_path in enumerate(sorted(nims_paths)):
    print("BIDSifying path %s out of %s" % (str(i+1), str(len(nims_paths))))
    subj_path  = get_subj_path(nims_path, bids_dir, id_correction_dict)
    if subj_path == None:
        print("Couldn't find subj_path for %s" % nims_path)
        with open(error_file, 'a') as filey:
            filey.write("Couldn't find subj_path for %s\n" % nims_path)
        continue
    bids_subj(subj_path, bids_dir, nims_path)

# add physio metadata
if not os.path.exists(os.path.join(bids_dir, 'recording-cardiac_physio.json')):
    if len(glob.glob(os.path.join(bids_dir, 'sub-*', 'ses-*', 'func', '*cardiac*'))) > 0:
        json.dump(cardiac_bids, open(os.path.join(bids_dir, 'recording-cardiac_physio.json'),'w'))
if not os.path.exists(os.path.join(bids_dir, 'recording-respiratory_physio.json')):
    if len(glob.glob(os.path.join(bids_dir, 'sub-*', 'ses-*', 'func', '*respiratory*'))) > 0:
        json.dump(respiratory_bids, open(os.path.join(bids_dir, 'recording-respiratory_physio.json'),'w'))

# *****************************
# *** Cleanup
# *****************************
cleanup(bids_dir)
