install:
	cpanm BioPerl --force
	cpanm Data::Dumper
	cpanm Getopt::Long
	cpanm File::Temp
clean:
	find ./ -name "*.DS_Store" -delete
	find ./ -name "*.log" -delete
	find ./ -name "*.gz" -delete
