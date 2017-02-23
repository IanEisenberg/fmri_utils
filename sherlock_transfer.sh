data_folder=$1

rsync -aHvz --progress /data/$data_folder sherlock:/scratch/users/ieisenbe/Data
