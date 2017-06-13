destination=sherlock:/scratch/PI/russpold/data/uh2
record=/data/uh2/sherlock_completed_files.txt

# loop through all uh2 files
for f in /nimsfs/russpold/uh2/*
do
	if grep -Fxq "$f" $record; then
		echo $f already transferred
	else		
		# run bids_organizer, which moves files to /data/uh2
		python bids_organizer.py  $f uh2_id_correction.json
		# move data
		if rsync -avz --progress /data/uh2/sub* $destination ; then
			echo $f >> $record
			# rm bids folder so we don't waste space
		else
			echo "Transfer failed"
		fi
		rm -r /data/uh2/sub*
	fi
done


# move non subject specific folders (e.g.task meta data)
rsync -avz --progress /data/uh2/task* $destination
rsync -avz --progress /data/uh2/data* $destination

