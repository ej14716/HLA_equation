###################################################################################
## J-Coupling
echo "Starting J coupling calculation"



# Awk script (Claire Dickson)
read -d "" scriptvariable << EOF
BEGIN{}

{
  if ((\$1=="Copyright"))
  {
    ++m
	  for (k=1; k<=100; k++)
    {
        getline
	      if (\$2=="NMR") {n=\$1}
    }
  }

  if ((\$4=="coupling") && (\$5=="J"))
  {
    printf molnam " " n >> molnam"_J_Summary.txt"
    print "" >> molnam"_J_Summary.txt"
    for (i=1; i<500; i++)
    {
			getline
      if (\$1=="End") {imax=i-1; i=510}
      else {a[i]=\$1 ; b[i]=\$2 ; c[i]=\$3 ; d[i]=\$4 ; e[i]=\$5 ; f[i]=\$6}
		}

		for (i=1; i<=imax; i++)
		{
		    if (f[i]==null && e[i]>1) {printf "\\\t" >> molnam"_J_Summary.txt"}
		    printf a[i] "\\\t" >> molnam"_J_Summary.txt"
		    printf b[i] "\\\t" >> molnam"_J_Summary.txt"
		    printf c[i] "\\\t" >> molnam"_J_Summary.txt"
		    printf d[i] "\\\t" >> molnam"_J_Summary.txt"
		    printf e[i] "\\\t" >> molnam"_J_Summary.txt"
		    printf f[i] "\\\t" >> molnam"_J_Summary.txt"
        print "" >> molnam"_J_Summary.txt"
		}
	}
}
END{}
EOF

mkdir J_matrix

# for each NMR logfile in the directory, run the awk script to get the J coupling summary
for logfile in $( ls -v pass/*NMR*.log); do
  #get filename
  molname=$(echo $logfile | tr "_" "\n" )
  molname=$(echo $molname | awk '{print $2}')

  {
  rm ${molname}_J_Summary.txt
  rm ${molname}_J_Raw.txt
  } &>/dev/null


  gawk -v molnam=$molname "$scriptvariable" $logfile

  # Convert the D exponential values to E, python doesnt recognise the former
  sed -i 's/D/E/g' ${molname}_J_Summary.txt
  # This line is stupid, i have no idea why this file gets made, but for now this just deletes it
  #rm ${molname}_J_Summary.txt-e

  # Run the coupling.py python script to convert the basic output to matrix format, and calculate boltzmann weighted values
  python J_coupling.py $molname
  code=$?
  if [ $code -eq 1 ]
  then
    exit 1
  fi

  mv ${molname}_J_Summary.txt J_matrix/
  mv ${molname}_J_Raw.txt J_matrix/

done


echo "J coupling processing done"
