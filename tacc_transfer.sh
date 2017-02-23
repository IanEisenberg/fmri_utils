data_folder=$1

rsync -aHvz --progress /data/$data_folder ls5:/corral-repl/utexas/poldracklab/data/ian
