#!/bin/bash
#


home=$(pwd)

{

rm -r passed
rm -r pass

mkdir pass
mkdir fail
mkdir running
mkdir xyz
} &>/dev/null

# check Usage
USAGE="Usage: ML_NMR_process.sh <Dataset directory>"
if [ ! $# -eq 1 ];
then
	echo $USAGE
	exit
fi
dataset=$1


if [ -f nmr_process.log ]
then
  rm nmr_process.log
fi
if [ -f running/errors.log ]
then
  rm running/errors.log
fi

echo NMR_process.sh > nmr_process.log

for file_path in $dataset/*NMR*; do
  file=$(basename $file_path)
#for file in HC064_CUQMOS_OPT_NMR_58.log; do
  cd $home
  cp $dataset/$file running/
  echo running $file >> nmr_process.log
  echo running $file
  cd running/

  Error_count=$(grep -c "Error" $file)
  size_count=$(du -k $file | awk '{print $1}')
  if [ $Error_count -gt 0 ]
  then
      rm $file
      cd $home
      echo $file failed, NMR calculation failed previously
      echo $file failed, NMR calculation failed previously >> nmr_process.log
      continue
  fi
  if [ $size_count -lt 100 ]
  then
      rm $file
      cd $home
      echo $file failed, NMR calculation failed previously
      echo $file failed, NMR calculation failed previously >> nmr_process.log
      continue
  fi

  cut1=$(grep -n "Center     Atomic      Atomic             Coordinates (Angstroms)" $file | tail -1 | awk '{print $1}' | sed 's/://')
  cut2=$(tail -n +$cut1 $file |  grep -n "\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-" | tail -1 | awk '{print $1}' | sed 's/://')
  cut3=$(($cut1))
  cut4=$(($cut2+$cut1-2))

  xyz_name=$( echo $file | sed 's/.log/.xyz/')

  sed -n $cut3,${cut4}p $file > $xyz_name

  python ../ML_NMR_process.py $xyz_name
  code=$?

  if [ $code -eq 0 ]
  then
    {
    mv *.txt ../pass/
    mv *.log ../pass/
    rm $file
    mv $xyz_name ../xyz/
    } &>/dev/null
    cd $home
    echo $file successful
    echo successful >> nmr_process.log
  fi
  if [ $code -eq 3 ]
  then
    {
    rm $file
    mv $xyz_name ../xyz/
    } &>/dev/null
    cd $home
    echo $file failed, correct bond order could not be established
    echo failed, correct bond order could not be established >> nmr_process.log

  fi
  if [ $code -eq 1 ]
  then
    {
    rm $file
    mv $xyz_name ../xyz/
    } &>/dev/null
    cd $home
    echo failed, python script error
  fi
  #rm $xyz_name
  #rm $file

done

total_ran=$(grep -c "running" nmr_process.log)
total_good=$(grep -c "successful" nmr_process.log)
total_nan=$(grep -c "previously" nmr_process.log)
total_bond=$(grep -c "bond order" nmr_process.log)

possible=$(($total_ran-$total_nan))


echo RESULT =================================================
echo TOTAL RAN: $total_ran
echo TOTAL POSSIBLE: $possible
echo BAD BOND ORDER: $total_bond
echo SUCCESSFUL: $total_good
