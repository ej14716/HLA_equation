# ! /bin/bash
# Will Gerrard  wg12385@bristol.ac.uk
#
#
#
# bash script to prepare data_train and data_test datasets for running in a linear regression algorithm
#
# Needs ML_NMR_process to already have been run in the same directory, so there is a folder named pass containing the dangle_mat files
# Needs J coupling matrix files inside J_matrix folder

{
  rm -r ML_raw_files
  mkdir ML_raw_files
  cp pass/*dangle.txt ML_raw_files/
  cp pass/*atno.txt ML_raw_files/
  cp J_matrix/*Raw.txt ML_raw_files/
  rm D_J.Data
} &>/dev/null


for dfile in $( ls -v ML_raw_files/*dangle.txt); do
  #get filename
  molname=$(echo $dfile | tr "_" "\n" )
  molname=$(echo $molname | awk '{print $4}')

  echo processing $molname

  nfile=$( echo $dfile | sed 's,dangle,atno,')
  jfile="ML_raw_files/${molname}_J_Raw.txt"

  if [ ! -f $jfile ]
  then
    echo "No matching J coupling file for $molname (${molname}_J_Raw.txt)"
    continue
  fi
  if [ ! -f $nfile ]
  then
    echo "No matching J coupling file for $molname (${molname}_J_Raw.txt)"
    continue
  fi
  python reg_dataset.py $molname $dfile $jfile $nfile
  code=$?
  if [ $code -eq 1 ]
  then
    exit 1
  fi


done
