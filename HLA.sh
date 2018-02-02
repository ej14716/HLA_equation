for file in $( ls -v test/ ); do
	# get molname from filename using awk or something
	python HLA.py $name
