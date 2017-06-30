# TACC
# destination=weisian@login1.corral.tacc.utexas.edu:/corral-repl/utexas/poldracklab/data/uh2
# record=/data/uh2/TACC_completed_files.txt
# sherlock
destination=sherlock:/oak/stanford/groups/russpold/data/uh2
record=/data/uh2/sherlock_completed_files.txt

# loop through all uh2 files
f=/nimsfs/russpold/uh2/*$1*
if grep -Fxq "$f" $record; then
	echo $f already transferred
else		
	# run bids_organizer, which moves files to /data/uh2
	python bids_organizer.py  $f --rsync_output test --id_correction uh2_id_correction.json --record $record

	rm -r /data/uh2/sub*
fi

sort $record -o $record

# move non subject specific folders (e.g.task meta data)
# rsync -avz --progress /data/uh2/task* $destination
# rsync -avz --progress /data/uh2/data* $destination

