import dicom
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
def bids_anat(subj_id, anat_dir, anat_path):
	#T1 or T2?
	anat_type = anat_dir[(anat_dir.index('anat_')+5):]
	anat_file = glob.glob(os.path.join(anat_dir,'*nii.gz'))
	assert len(anat_file) == 1, "More than one anat file found in directory %s" % anat_dir
	new_file = os.path.join(anat_path, subj_id + '_' + anat_type + '.nii.gz')
	if not os.path.exists(new_file):
		shutil.copyfile(anat_file[0], new_file)
	else:
		print('Did not save anat because %s already exists!' % new_file)

def bids_fmap(subj_id, fmap_dir, fmap_path, func_path):
	fmap_files = glob.glob(os.path.join(fmap_dir,'*.nii.gz'))
	assert len(fmap_files) == 2, "Didn't find the correct number of files in %s" % fmap_dir
	fmap_index = [0,1]['fieldmap.nii.gz' in fmap_files[1]]
	fmap_file = fmap_files[fmap_index]
	mag_file = fmap_files[1-fmap_index]
	fmap_name = subj_id + '_' + 'fieldmap.nii.gz'
	mag_name = subj_id + '_' + 'magnitude.nii.gz'
	if not os.path.exists(os.path.join(fmap_path, fmap_name)):
		shutil.copyfile(fmap_file, os.path.join(fmap_path, fmap_name))
		shutil.copyfile(mag_file, os.path.join(fmap_path, mag_name))
		func_runs = [os.path.join('func',os.path.basename(f)) for f in glob.glob(os.path.join(func_path,'*task*bold.nii.gz'))]
		fieldmap_meta = {'Units': 'Hz', 'IntendedFor': func_runs}
		json.dump(fieldmap_meta,open(os.path.join(fmap_path, subj_id + '_fieldmap.json'),'w'))
	else:
		print('Did not save fmap_epi because %s already exists!' % os.path.join(fmap_path, fmap_name))

def bids_fmap_epi(subj_id, fmap_epi_dir, fmap_path, func_path):
	fmap_epi_file = glob.glob(os.path.join(fmap_epi_dir,'*.nii.gz'))
	assert len(fmap_epi_file) <= 1, "More than one func file found in directory %s" % fmap_epi_dir
	if len(fmap_epi_file) == 0:
		print('Skipping %s, no nii.gz file found' % fmap_epi_dir)
		return
	# bring to subject directory and divide into sbref and bold
	new_fmap_file = os.path.join(fmap_path, subj_id + '_epi.nii.gz')
	shutil.copyfile(fmap_epi_file[0], new_fmap_file)
	func_runs = [os.path.join('func',os.path.basename(f)) for f in glob.glob(os.path.join(func_path,'*task*bold.nii.gz'))]
	meta_path = os.path.join(fmap_path,subj_id + 'epi.json')
	if not os.path.exists(meta_path):
		meta_file = [x for x in glob.glob(os.path.join(fmap_epi_dir,'*.json')) if 'qa' not in x][0]
		fmap_epi_meta = get_fmap_epi_meta(meta_file,func_runs)
		json.dump(fmap_epi_meta,open(meta_path,'w'))
	else:
		print('Did not save fmap because %s already exists!' % meta_path)

def bids_sbref(subj_id, sbref_dir, func_path, data_path):
	task_index = sbref_dir.index('task')
	filename = sbref_dir[task_index:]
	sbref_files = glob.glob(os.path.join(sbref_dir,'*.nii.gz'))
	assert len(sbref_files) <= 1, "More than one func file found in directory %s" % sbref_dir
	if len(sbref_files) == 0:
		print('Skipping %s, no nii.gz file found' % sbref_dir)
		return

	# bring to subject directory and divide into sbref and bold
	sbref_file = os.path.join(func_path, subj_id + '_' + filename + '.nii.gz')
	# check if file exists. If it does, check if the saved file has more time points
	if os.path.exists(sbref_file):
		print('%s already exists!' % sbref_file)
		saved_shape = nib.load(sbref_file).shape
		current_shape = nib.load(sbref_files[0]).shape
		print('Dimensions of saved image: %s' % list(saved_shape))
		print('Dimensions of current image: %s' % list(current_shape))
		if (current_shape[-1] < saved_shape[-1]):
			print('Current image has fewer time points than saved image. Exiting...')
			return
		else:
			print('Current image has more time points than saved image. Overwriting...')
	# save sbref image to bids directory
	shutil.copyfile(sbref_files[0], sbref_file)
	# get metadata
	sbref_meta_path = os.path.join(data_path, re.sub('_run[-_][0-9]','',filename) + '.json')
	if not os.path.exists(sbref_meta_path):
		try:
			meta_file = [x for x in glob.glob(os.path.join(sbref_dir,'*.json')) 
							if 'qa' not in x][0]
			func_meta = get_functional_meta(meta_file, filename)
			json.dump(func_meta,open(sbref_meta_path,'w'))
		except IndexError:
			print("Metadata couldn't be created for %s" % sbref_file)


def bids_task(subj_id, task_dir, func_path, data_path):
	task_index = task_dir.index('task')
	taskname = task_dir[task_index:]
	task_file = glob.glob(os.path.join(task_dir,'*.nii.gz'))
	assert len(task_file) <= 1, "More than one func file found in directory %s" % task_dir
	if len(task_file) == 0:
		print('Skipping %s, no nii.gz file found' % task_dir)
		return

	# bring to subject directory and divide into sbref and bold
	bold_file = os.path.join(func_path, subj_id + '_' + taskname + '.nii.gz')
	bold_file = bold_file.replace('_ssg', '_bold')
	# check if file exists. If it does, check if the saved file has more time points
	if os.path.exists(bold_file):
		print('%s already exists!' % bold_file)
		saved_shape = nib.load(bold_file).shape
		current_shape = nib.load(task_file[0]).shape
		print('Dimensions of saved image: %s' % list(saved_shape))
		print('Dimensions of current image: %s' % list(current_shape))
		if (current_shape[-1] < saved_shape[-1]):
			print('Current image has fewer time points than saved image. Exiting...')
			return
		else:
			print('Current image has more time points than saved image. Overwriting...')
	# save bold image to bids directory
	shutil.copyfile(task_file[0], bold_file)
	# get epi metadata
	bold_meta_path = os.path.join(data_path, re.sub('_run[-_][0-9]','',taskname) + '_bold.json')
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
        	new_physio_file = taskname + '_recording-' + pname
        	f = np.loadtxt(pfile)
        	np.savetxt(os.path.join(func_path, new_physio_file + '_physio.tsv.gz'), f, delimiter = '\t')
        shutil.rmtree(os.path.join(func_path, physio_file))
		

def cleanup(path):
	for f in glob.glob(os.path.join(path, '*')):
		new_name = f.replace('task_', 'task-').replace('run_','run-')
		os.rename(f,new_name)

def get_functional_meta(json_file, taskname):
	meta_file = json.load(open(json_file,'r'))
	meta_data = {}
	mux = meta_file['num_slices']
	nslices = meta_file['num_slices'] * mux
	tr = meta_file['tr']
	n_echoes = meta_file['acquisition_matrix_y'] 

	# fill in metadata
	meta_data['TaskName'] = taskname
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
	meta_file = json.load(open(json_file,'r'))
	meta_data = {}
	mux = meta_file['num_slices']
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
        nslices = nslices/8
        mux_slice_acq_order = range(0,nslices,2) + range(1,nslices,2)
        mux_slice_acq_time = [float(s)/nslices*tr for s in xrange(nslices)]
        unmux_slice_acq_order = [nslices*m+s for m in xrange(mux) for s in mux_slice_acq_order]
        unmux_slice_acq_time = mux_slice_acq_time * mux
        slicetimes = zip(unmux_slice_acq_time,unmux_slice_acq_order)
    else:
        slice_acq_order = range(0,nslices,2) + range(1,nslices,2)
        slice_acq_time = [float(s)/nslices*tr for s in xrange(nslices)]
        slicetimes = zip(slice_acq_time,slice_acq_order)
    #reorder slicetimes by slice number
    sort_index = sorted(enumerate([x[1] for x in slicetimes]), key= lambda x: x[1])
    sorted_slicetimes = [slicetimes[i[0]][0] for i in sort_index]
    return sorted_slicetimes

def get_subj_path(nims_file, id_correction_dict=None):
	meta_json = glob.glob(os.path.join(nims_file,'*','*1.json'))[0]
	meta_file = json.load(open(meta_json,'r'))
	exam_number = meta_file['exam_number']
	sub_session = str(meta_file['patient_id'].split('@')[0])
	subj_id = sub_session.split('_')[0]
	# correct sub_id if provided
	if id_correction_dict:
		subj_id = id_correction_dict.get(exam_number,subj_id)
	session = '1'
	if '_' in sub_session:
		session = sub_session.split('_')[1]
	subj_path = os.path.join(data_path, 'sub-'+subj_id, 'session_'+session)
	return subj_path

def mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
    	print('Directory %s already existed' % path)
    return path

# *****************************
# *** Main BIDS function
# *****************************
def bids_subj(subj_path, data_path, nims_path):
	if os.path.exists(subj_path):
		print("Path %s already exists. Skipping." % subj_path)
	else:
		print("********************************************")
		print("BIDSifying %s" % subj_path)
		print("********************************************")
		# extract subject ID
		split_path = os.path.normpath(subj_path).split(os.sep)
		subj_id = [x for x in split_path if 'sub' in x][0]
		# split subject path into a super subject path and a session path
		session_path = subj_path
		subj_path = os.path.split(subj_path)[0]
		mkdir(subj_path)
		mkdir(session_path)
		anat_path = mkdir(os.path.join(session_path,'anat'))
		func_path = mkdir(os.path.join(session_path,'func'))
		fmap_path = mkdir(os.path.join(session_path,'fmap'))

		#header file
		header = {'Name': study_id, 'BIDSVersion': '1.51-rc1'}
		json.dump(header,open(os.path.join(data_path, 'dataset_description.json'),'w'))


		print('BIDSifying anatomy...')
		# anat files
		anat_dirs = sorted(glob.glob(os.path.join(nims_path,'*anat*')))[::-1]
		for anat_dir in anat_dirs:
			print('\t' + anat_dir)
			bids_anat(subj_id, anat_dir, anat_path)

		print('BIDSifying sbref...')
		# task files
		sbref_dirs = sorted(glob.glob(os.path.join(nims_path,'*sbref*')))[::-1]
		for sbref_dir in sbref_dirs:
			print('\t' + sbref_dir)
			bids_sbref(subj_id, sbref_dir, func_path, data_path)

		print('BIDSifying task...')
		# task files
		task_dirs = sorted(glob.glob(os.path.join(nims_path,'*task*')))[::-1]
		for task_dir in [x for x in task_dirs if 'sbref' not in x]:
			print('\t' + task_dir)
			bids_task(subj_id, task_dir, func_path, data_path)
		# cleanup
		cleanup(func_path)

		print('BIDSifying fmap...')
		#fmap_epi files
		fmap_epi_dirs = sorted(glob.glob(os.path.join(nims_path,'*fmap_epi*')))[::-1]
		for fmap_epi_dir in fmap_epi_dirs:
			print('\t' + fmap_epi_dir)
			bids_fmap_epi(subj_id, fmap_epi_dir, fmap_path, func_path)

		# fmap files
		fmap_dirs = sorted(glob.glob(os.path.join(nims_path,'*fieldmap*')))[::-1]
		for fmap_dir in fmap_dirs:
			print('\t' + fmap_dir)
			bids_fmap(subj_id, fmap_dir, fmap_path, func_path)

		
		
# *****************************
# *** ORGANIZE IN BIDS
# *****************************

nims_paths = sys.argv[1:]
#study name
study_id = nims_paths[0].split('/')[3]
# get data directory to save bids in
data_path = os.path.join('/data',study_id)
mkdir(data_path)
# set id_correction_dict if provided
id_correction_dict = None
if 'json' in sys.argv[-1]:
	id_correction_dict = json.load(open(sys.argv[-1],'r'))
	nims_paths = nims_paths[:-1]
print('Using ID correction json file: %s' % sys.argv[-1])
# bidsify all subjects in path
for nims_file in sorted(nims_paths):
	print(nims_file)
	#try:
	subj_path  = get_subj_path(nims_file, id_correction_dict)
	if subj_path == None:
		print("Couldn't find subj_path for %s" % nims_file)
		continue
	bids_subj(subj_path, data_path, nims_file)
	#except UnboundLocalError:
	#	print("BIDSifying failed on %s!" % nims_file)

# add physio metadata
if not os.path.exists(os.path.join(data_path, 'recording-cardiac_physio.json')):
	if len(glob.glob(os.path.join(data_path, 'sub-*', 'func', '*cardiac*'))) > 0:
		json.dump(cardiac_bids,open(os.path.join(data_path, 'recording-cardiac_physio.json'),'w'))
if not os.path.exists(os.path.join(data_path, 'recording-respiratory_physio.json')):
	if len(glob.glob(os.path.join(data_path, 'sub-*', 'func', '*respiratory*'))) > 0:
		json.dump(respiratory_bids,open(os.path.join(data_path, 'recording-respiratory_physio.json'),'w'))
# *****************************
# *** Cleanup
# *****************************
cleanup(data_path)
