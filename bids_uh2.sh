
destination=sherlock:/oak/stanford/groups/russpold/data/uh2
record=/data/uh2/sherlock_completed_files.txt

# loop through all uh2 files
for f in /nimsfs/russpold/uh2/*
do
	if grep -Fxq "$f" $record; then
		echo $f already transferred
	else		
		# run bids_organizer, which moves files to /data/uh2
		python bids_organizer.py  $f --rsync_output $record --id_correction uh2_id_correction.json --record $record

		rm -r /data/uh2/sub*
	fi
done

sort $record -o $record


# move non subject specific folders (e.g.task meta data)
# rsync -avz --progress /data/uh2/task* $destination
# rsync -avz --progress /data/uh2/data* $destination

