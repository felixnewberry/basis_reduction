#!/bin/bash

function setPhtField {
  echo "Writing to $3: $1 --> $2"
  sed -i "s/$1\=\".*\"/$1\=\"$2\"/g" "$3"
}

function writePht {
  TAGS=("number_of_steps" "start_index")
  FILE="$1"
  VALUES=("${@:2}") # All parameters but the first
  local j
  for (( j=0; j<${#TAGS[@]}; j++ )); do
    setPhtField "${TAGS[$j]}" "${VALUES[$j]}" "$FILE"
  done
}

### BEGIN SCRIPT ###

# Sanitize input
NPROCS="16"
RUN_ID="felixTest"
PHT_FILE="flowPosix_extract_1.phts"
DIR_INCREMENT=1
NUM_START=183
NUM_RUN_DIRS=18
n_time="1"
time_spacing=50

if [ -z "$1" ] || (("$1" < 1)) || (("$1" > "$DIR_INCREMENT")); then
  echo "Need an input between 1 and $DIR_INCREMENT; this is the starting offset for which directory to run."
  exit
else
  DIR_OFFSET="$1"
fi

#
# Loop over run directories
#
for (( i=$(($NUM_START+$DIR_OFFSET-1)); i<$(($NUM_START+$NUM_RUN_DIRS)); i+="$DIR_INCREMENT" )); do

  # (***FORMAT SPECIFIES LEADING ZEROS IN Run-0001 etc. DIRECTORIES***)
  d=$(printf "Run-%03i" $i)
  echo "Extracting run data from $d (offset is $DIR_OFFSET)"
  cd "$d"

    cd "$NPROCS-procs_case"
    

      # Determine most recently-written restart file (temporally, through 'ls -rt')
      lastTimeStep_orig=$(ls -rt restart-dat.* | tail -n1 | sed 's/\./\n/g' | head -n2 | tail -n1)
    cd ..

    for  (( i_time=0;i_time<=$n_time-1;i_time++));do 
    lastTimeStep=$((lastTimeStep_orig-$i_time*$time_spacing))

    vals=("1" "$lastTimeStep")
    writePht "$PHT_FILE" "${vals[@]}"
    pvpython ../extract_cylinder_temperature.py $(printf "%s-ts%i-%03i" $RUN_ID $lastTimeStep $i) "$PHT_FILE" "$lastTimeStep"

    writePht "$PHT_FILE" "${vals[@]}"
    pvpython ../extract_end_temperature.py $(printf "%s-ts%i-%03i" $RUN_ID $lastTimeStep $i) "$PHT_FILE" "$lastTimeStep"
           
    writePht "$PHT_FILE" "${vals[@]}"
    pvpython ../extract_mid_temperature.py $(printf "%s-ts%i-%03i" $RUN_ID $lastTimeStep $i) "$PHT_FILE" "$lastTimeStep"
    done

  cd ..

done
