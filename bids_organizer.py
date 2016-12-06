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
	new_name = subj_id + '_' + anat_type + '.nii.gz'
	shutil.copyfile(anat_file[0], os.path.join(anat_path, new_name))

def bids_fmap(subj_id, fmap_dir, fmap_path, func_path):
	fmap_files = glob.glob(os.path.join(fmap_dir,'*.nii.gz'))
	assert len(fmap_files) == 2, "Didn't find the correct number of files in %s" % fmap_dir
	fmap_index = [0,1]['fieldmap.nii.gz' in fmap_files[1]]
	fmap_file = fmap_files[fmap_index]
	mag_file = fmap_files[1-fmap_index]
	fmap_name = subj_id + '_' + 'fieldmap.nii.gz'
	mag_name = subj_id + '_' + 'magnitude.nii.gz'
	shutil.copyfile(fmap_file, os.path.join(fmap_path, fmap_name))
	shutil.copyfile(mag_file, os.path.join(fmap_path, mag_name))
	func_runs = [os.path.join('func',os.path.basename(f)) for f in glob.glob(os.path.join(func_path,'*task*bold.nii.gz'))]
	fieldmap_meta = {'Units': 'Hz', 'IntendedFor': func_runs}
	json.dump(fieldmap_meta,open(os.path.join(fmap_path, subj_id + '_fieldmap.json'),'w'))

def bids_fmap_epi(subj_id, fmap_epi_dir, fmap_path, func_path):
	fmap_epi_file = glob.glob(os.path.join(fmap_epi_dir,'*.nii.gz'))
	assert len(fmap_epi_file) <= 1, "More than one func file found in directory %s" % fmap_epi_dir
	if len(fmap_epi_file) == 0:
		print('Skipping %s, no nii.gz file found' % fmap_epi_dir)
		return
	# bring to subject directory and divide into sbref and bold
	new_fmap_file = os.path.join(func_path, subj_id + '_epi.nii.gz')
	shutil.copyfile(fmap_epi_file[0], new_fmap_file)
	epi_img = nib.load(new_fmap_file)
	func_runs = [os.path.join('func',os.path.basename(f)) for f in glob.glob(os.path.join(func_path,'*task*bold.nii.gz'))]
	meta_path = os.path.join(fmap_path,'epi.json')
	if not os.path.exists(meta_path):
		fmap_epi_meta = get_fmap_epi_meta(epi_img.header,func_runs)
		json.dump(fmap_epi_meta,open(meta_path,'w'))
	else:
		print('Did not save fmap_epi because %s already exists!' % meta_path)

def bids_task(subj_id, task_dir, func_path, data_path):
	task_index = task_dir.index('task')
	taskname = task_dir[task_index:]
	task_file = glob.glob(os.path.join(task_dir,'*.nii.gz'))
	assert len(task_file) <= 1, "More than one func file found in directory %s" % task_dir
	if len(task_file) == 0:
		print('Skipping %s, no nii.gz file found' % task_dir)
		return
	# bring to subject directory and divide into sbref and bold
	new_task_file = os.path.join(func_path, subj_id + '_' + taskname + '.nii.gz')
	shutil.copyfile(task_file[0], new_task_file)
	epi_img = nib.load(new_task_file)
	bold_name = subj_id + '_' + taskname + '_bold.nii.gz'
	sbref_name = subj_id + '_' + taskname + '_sbref.nii.gz'
	# extracting sbref 
	fsl.utils.ExtractROI(in_file = new_task_file,
						roi_file = os.path.join(func_path, sbref_name),
						t_min = 2, t_size = 1).run()
	# extracting bold activity
	fsl.utils.ExtractROI(in_file = new_task_file,
						roi_file = os.path.join(func_path, bold_name),
						t_min = 3, t_size = epi_img.shape[3]-2).run()
	os.remove(new_task_file)
	# get metadata
	meta_path = os.path.join(data_path, re.sub('_run[-_][0-9]','',taskname) + '_bold.json')
	if not os.path.exists(meta_path):
		func_meta = get_functional_meta(epi_img.header, taskname)
		json.dump(func_meta,open(meta_path,'w'))
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



def get_functional_meta(header, taskname):
	meta_data = {}
	descrip = str(header['descrip'])[0:].split(';')
	if descrip[0:2] == "b'":
		descrip = descrip[2:]
	descrip = {k:v for k,v in [i.split('=') for i in descrip]}
	acq = json.loads(descrip['acq'])
	phase_dim = header.get_dim_info()[1]
	tr = header['pixdim'][4]
	mux = 8
	nslices = header['dim'][3]
	# fill in metadata
	meta_data['TaskName'] = taskname
	meta_data['EffectiveEchoSpacing'] = float(descrip['ec'])/1000
	meta_data['EchoTime'] = float(descrip['te'])/1000 # check this with Chris/units
	meta_data['FlipAngle'] = descrip['fa']
	meta_data['RepetitionTime'] = round(tr,4)
	# slice timing
	meta_data['SliceTiming'] = get_slice_timing(nslices, tr, mux = mux)
	total_time = (acq[phase_dim]-1)*meta_data['EffectiveEchoSpacing']
	meta_data['TotalReadoutTime'] = total_time
	meta_data['PhaseEncodingDirection'] = ['i','j','k'][phase_dim] + '-'	
	return meta_data

def get_fmap_epi_meta(header, intended_list):
	meta_data = {}
	descrip = str(header['descrip'])[0:].split(';')
	if descrip[0:2] == "b'":
		descrip = descrip[2:]
	descrip = {k:v for k,v in [i.split('=') for i in descrip]}
	acq = json.loads(descrip['acq'])
	phase_dim = header.get_dim_info()[1]
	# fill in metadata
	meta_data['EffectiveEchoSpacing'] = float(descrip['ec'])/1000
	total_time = (acq[phase_dim]-1)*meta_data['EffectiveEchoSpacing']
	meta_data['TotalReadoutTime'] = total_time
	meta_data['PhaseEncodingDirection'] = ['i','j','k'][phase_dim] + '+'
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

def get_subj_path(raw_path, temp_dir, data_path):
	tar = tarfile.open(glob.glob(os.path.join(raw_path,'*3Plane_Loc*/*dicoms.tgz*'))[0])
	tar.extractall(temp_dir)
	dcms = glob.iglob(os.path.join(temp_dir,'*dicoms/*.dcm'))
	dicominfo = dicom.read_file(next(dcms))
	shutil.rmtree(glob.glob(os.path.join(temp_dir,'*dicoms'))[0])
	if '@' in dicominfo.PatientID:
		subj_id = dicominfo.PatientID.split('@')[0]
	elif 'russpold' not in dicominfo.PatientID:
		subj_id = dicominfo.PatientID
	assert len(subj_id) > 0, "No subject file found"
	subj_path = os.path.join(data_path,'sub-' + subj_id)
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
		subj_id = os.path.basename(os.path.normpath(subj_path))
		mkdir(subj_path)
		anat_path = mkdir(os.path.join(subj_path,'anat'))
		func_path = mkdir(os.path.join(subj_path,'func'))
		fmap_path = mkdir(os.path.join(subj_path,'fmap'))

		#header file
		header = {'Name': study_id, 'BIDSVersion': '1.51-rc1'}
		json.dump(header,open(os.path.join(data_path, 'dataset_description.json'),'w'))

		# anat files
		anat_dirs = glob.iglob(os.path.join(nims_path,'*anat*'))
		for anat_dir in anat_dirs:
			bids_anat(subj_id, anat_dir, anat_path)

		# task files
		task_dirs = glob.iglob(os.path.join(nims_path,'*task*'))
		for task_dir in task_dirs:
			bids_task(subj_id, task_dir, func_path, data_path)
		# cleanup
		cleanup(func_path)

		#fmap_epi files
		fmap_epi_dirs = glob.iglob(os.path.join(nims_path,'*fmap_epi*'))
		for fmap_epi_dir in fmap_epi_dirs:
			bids_fmap_epi(subj_id, fmap_epi_dir, fmap_path, func_path)

		# fmap files
		fmap_dirs = glob.iglob(os.path.join(nims_path,'*fieldmap*'))
		for fmap_dir in fmap_dirs:
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

#make temp directory 
temp_dir = os.path.join(data_path,'temp')
mkdir(temp_dir)

# bidsify all subjects in path
for nims_file in nims_paths:
	raw_path = os.path.join('/nimsfs', 'raw', nims_file.split('nimsfs/')[1])
	try:
		subj_path  = get_subj_path(raw_path, temp_dir, data_path)
		bids_subj(subj_path, data_path, nims_file)
	except UnboundLocalError:
		print("BIDSifying failed on %s!" % nims_file)

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
os.rmdir(temp_dir)
cleanup(data_path)
