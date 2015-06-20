#!/usr/bin/perl -w

#Cohn lab script for adding unified locus [chromosome]:[location], line number, and better ESP data fields to
#seattleseq outputs.  Also autogenerates the mySQL script required to read the data into a mySQL table.
#Copyright 2014 Michael Weinstein, Daniel H. Cohn laboratory, UCLA.

use strict;
use warnings;
use Getopt::Std;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use List::MoreUtils;
$|=1;


sub main(){

    my %opts;  #setup the hash for commandline options
    getopts('f:d:t:', \%opts);  #write the commandline options to the array
    unless(checkopts(\%opts)){  #check the options for validity
        print "Invalid or non-existent filename specified. \n";
        usage();  #if they are not valid (function returns false), print the usage instructions
    }
    my $seqin = $opts{f};  #sets up the variable for the file to be annotated using the commandline option
    my $seqout = $seqin."\.SQLready\.txt";  #sets up the variable for the output file for annotated data
    my $sqlout = $seqin.".sql";  #sets up the variable for the sql script file
    open(SSI, $seqin) or die "Couldn't open file $seqin. \n"; #opens the input file or quits
    if (-e $seqout) {  #checks to see if the annotation output file already exists
        die "\nOutput file $seqout already appears to exist\.\n"; #quits if it does
    }
    if (-e $seqout."\.log") { #checks to see if the annotation log file already exists
        die "\nLog file $seqout\.log already appears to exist\.\n"; #quits if it does
    }
    if (-e $seqout."\.sql") { #checks to see if the SQL script file already exists
        die "\nLog file $seqout\.sql already appears to exist\.\n"; #quits if it does
    }
    
 

    open(OUTPUT, '>'.$seqout) or die "\nError opening output file $seqout\.\n";  #Creates the file for annotation output
    open(LOGFILE, '>'.$seqout."\.log") or die "\nError opening log file $seqout\.log\.\n"; #Creates the log file
    my $header = <SSI>;  #reads the first line of the file (which should contain header data)  #reads the first seattleseq output line as the header
    my @header = split(/\t/, $header);  #splits the header on each tab into a separate variable in an array
    foreach my $identifier(@header){ #iterates through each entry in the header array
        $identifier =~ s/\R//g;  #removes any linebreaks from the header (including the break at the end)
        $identifier =~ s/\%/\_percent/; #replaces the percent symbol (forbidden in field identifiers) with the word
        $identifier =~ s/\W/\_/g; #replaces any disallowed characters in field identifiers with underscores
    }
    print OUTPUT "varMD5\tuniLoc\t";  #prints the first two fields in the header
    foreach my $identifier(@header){   #iterates through the header values in the header array
        print OUTPUT $identifier."\t";  #prints each one
    }
    print OUTPUT "ESPfreqs\tMAFinESP\n";  #prints the last two header values
    unless ($header =~ /chromosome\tposition/i){  #checks to see if the header looks anything like a seattleseq header
        die "File header for $seqin does not indicate chromosome and position data." #quits if it doesn't
    }
    my $progress = 0; #starts the progress counter at zero
    
    LINE: while (my $line = <SSI>){  #loop that goes through file while there is file to go through line by line
        
        if ($line =~ /^#/) {  #checks to see if the line starts with a #, indicating a data or comment line
            print LOGFILE "Skipped line $progress\; it appears to be a comment line\.\n"; #if so, prints it directly to the output file unmodified
            $progress ++; #increments the progress line
            next LINE; #and moves on to the next line
        }
        
        $line =~ s/\R//g;  #removes any linebreaks from the data line being read (including the one at the end)
        $progress ++;  #Increments the progress (lines processed) count
        my @line = split(/\t/, $line);  #splits the data line just like the header
        if (scalar(@line) != scalar(@header)) {  #tests to see if there as many data entries as there are headers
            print "\nMismatched header and data on line $progress.\n";  #Prints the error to the console
            print LOGFILE "\nMismatched header and data on line $progress.\n"; #Prints the error to the logfile
            next LINE;  #skips the line if there are not
        }
        
        my $index = 0; #resets the index of which column within the row is being moved to the hash
        my %linehash = ();  #creates a hash with the header entry as key and data as value
        
        foreach my $headervalue(@header){  #goes through each header entry
            $linehash{$header[$index] } = $line[$index];  #writes the header entry as key and data entry as value in the hash
            $index++;  #increments the index so that the next iteration will look at the next column in each line
        }
        #Creating a hash value for the variant here
        my $hashuniloc = $linehash{chromosome}.":".$linehash{position};
        my $hashref = $linehash{referenceBase};
        my $hashalts = $linehash{sampleAlleles};
        my @hashalts = split(/\//, $hashalts);  #splitting reported alt alleles on the slash that separates them (if any)
        my %hashalts = map {$_ => "meaningless value"} @hashalts;
        my @hashaltsfiltered = keys %hashalts;
        my $hashvariantstring = $hashuniloc."/".$hashref;
        my $alreadyaddedavariant = 0;
        foreach my $potentialvariant(@hashaltsfiltered){
            if ($potentialvariant ne $hashref) {
                if ($alreadyaddedavariant) {
                    $hashvariantstring = $hashvariantstring.",".$potentialvariant;
                }
                else{
                    $hashvariantstring = $hashvariantstring."/".$potentialvariant;
                    $alreadyaddedavariant = "True"
                }
            }
            
        }
        my $varhash = md5_hex($hashvariantstring);
        #Done creating a hash for the variant
        print OUTPUT $varhash."\t".$linehash{chromosome}."\:".$linehash{position}."\t";  #prints the first and second fields of a data row, generates the line number trough the counter and the uniLoc by concatenating together the chromosome and locus fields
        foreach my $headervalue(@header){  #iterates through each value in the header array
            print OUTPUT $linehash{$headervalue}."\t";  #prints to the file the corresponding value from the hash of line values
        }
        
        my %ESP = ();
        $ESP{A} = 0;
        $ESP{G} = 0;
        $ESP{C} = 0;
        $ESP{T} = 0; #initializes the ESP hash to all zeros
        my %ESPfreq = ();
        $ESPfreq{A} = 0;
        $ESPfreq{G} = 0;
        $ESPfreq{C} = 0;
        $ESPfreq{T} = 0; #initializes the ESP hash to all zeros  #initializes the ESPfreq hash to all zeros

        
        if ($linehash{genomesESP} =~ /unknown/i) {  #checks if the ESP data exists
            print OUTPUT "unknown A0.000000\/T0.000000\/G0.000000\/C0.000000\t0\.0\n";  #if not, enters zero values and marks it as unknown
        }
        else{
            my @ESP = split(/\//, $linehash{genomesESP});  #creates an array of ESP values by spliting the different allele listings
            foreach my $entry(@ESP){  #goes through each of the ESP alleles
                my @entry = split (/\=/, $entry);  #splits each one on the equal sign to get a base and number of times observed
                $ESP{$entry[0]} = $entry[1];  #puts this value in the hash under the appropriate key
            }
            my $ESPtotal = $ESP{A} + $ESP{T} + $ESP{G} + $ESP{C};  #counts the total number of times this locus was studied in ESP
            $ESPfreq{'A'} = sprintf '<%.6f>', $ESP{A}/$ESPtotal;  #returns a frequency in the form of a 6 decimal place number (surrounded by angle brackets)
            $ESPfreq{'A'} =~ s/[\<\>]//g; #removes angle brackets from the value
            $ESPfreq{'T'} = sprintf '<%.6f>', $ESP{T}/$ESPtotal;
            $ESPfreq{'T'} =~ s/[\<\>]//g;
            $ESPfreq{'G'} = sprintf '<%.6f>', $ESP{G}/$ESPtotal;
            $ESPfreq{'G'} =~ s/[\<\>]//g;
            $ESPfreq{'C'} = sprintf '<%.6f>', $ESP{C}/$ESPtotal;
            $ESPfreq{'C'} =~ s/[\<\>]//g;
            $ESPfreq{'N'} = 1;
            print OUTPUT "A".$ESPfreq{A}."T".$ESPfreq{T}."G".$ESPfreq{G}."C".$ESPfreq{C}."\t";  #prints the ESP data as a single field
        
            my $rarestfreq = 1.0;  #initializes the rarestfreq value to 1
            my @samplealleles = split(/\//, $linehash{sampleAlleles});  #splits the field containing the sample alleles into individual bases
            foreach my $allele(@samplealleles){  #iterates through each allele
                if ($ESPfreq{$allele} < $rarestfreq) {  #if the frequency of the allele currently being looked at is less than the previous lowest value, this writes that value to the variable
                    $rarestfreq = $ESPfreq{$allele};  #actually changes the variable
                }    
            }    
            print OUTPUT $rarestfreq."\n";  #prints the frequency in ESP of the rarest allele we observed
        }
    print "\rProcessed $progress lines\.";  #updates the progress counter

    }
    close OUTPUT;  #closes the output file
    close LOGFILE;  #closes the log file
    
    
    print "\nDone generating the input file\, now generating the SQL script\!\n";  #prints done in the console
    
    
    #my $dbname;  #declares a variable for database name     SECTION COMMENTED OUT FOR REPLACEMENT WITH A SECTION THAT CAN USE COMMANDLINE ARGUMENTS FOR NAMES COPIED FROM VCF2SQL4.1
    #my $tablename;  #declares a variable for table name
    #while (!defined $dbname) {  #until the table name is defined
    #    print "In which database will we create this table\?\n";
    #    $dbname = readline STDIN;  #reads the user input line
    #    chomp $dbname;  #eliminates trailing characters and linebreaks from the input
    #    if ($dbname =~ /\W/) {  #checks the input for any non-word characters (anything other than alphaneumeric and underscore)
    #        print "Database name contains invalid characters\.  Valid characters are alphaneumeric or underscore\.\n";
    #        undef $dbname;  #undefines the database name before returning to the beginning of the loop
    #    }   
    #}
    #print "Please be sure that the database is already created before running the table creation script generated here\.\n";
    #while (!defined $tablename) {  #this block of code does the same thing with the table name as we just did with the database name
    #    print "What do you wish to call the table\?\n";
    #    $tablename = readline STDIN;
    #    chomp $tablename;
    #    if ($tablename =~ /\W/) {
    #        print "Table name contains invalid characters\.  Valid characters are alphaneumeric or underscore\.\n";
    #        undef $tablename;
    #    }
    #}
    
    my $dbname; #declares a variable for database name
    if ($opts{d}) {
        $dbname = $opts{d};
        if ($dbname =~ /\W/) {  #checks the input for any non-word characters (anything other than alphaneumeric and underscore)
            print "Database name contains invalid characters\.  Valid characters are alphaneumeric or underscore\.\n";
            undef $dbname;  #undefines the database name before returning to the beginning of the loop
        }
        else {print "Using database ".$dbname."\n"}
    }
    my $tablename;  #declares a variable for table name
    if ($opts{t}) {
        $tablename = $opts{t};
        if ($tablename =~ /\W/) {  #checks the input for any non-word characters (anything other than alphaneumeric and underscore)
            print "Database name contains invalid characters\.  Valid characters are alphaneumeric or underscore\.\n";
            undef $tablename;  #undefines the database name before returning to the beginning of the loop
        }
        else {print "Using table ".$tablename."\n"}
    }
    

    unless (-e "seattleSequel.prefs.txt"){  #checks to see if the preferences file exists
        print "No field\-type preference file found\.  Generating default preference file\.\n";  #prints a message if not
        open (PREFS, ">seattleSequel.prefs.txt") or die "Unable to create preferences file.\n";  #creates a new preferences file using the two lines below
        print PREFS "lineNumber\tuniLoc\t__inDBSNPOrNot\tchromosome\tposition\treferenceBase\tsampleGenotype\tsampleAlleles\tallelesDBSNP\taccession\tfunctionGVS\tfunctionDBSNP\trsID\taminoAcids\tproteinPosition\tcDNAPosition\tpolyPhen\tgranthamScore\tscorePhastCons\tconsScoreGERP\tscoreCADD\tchimpAllele\tgeneList\tAfricanHapMapFreq\tEuropeanHapMapFreq\tAsianHapMapFreq\thasGenotypes\tdbSNPValidation\trepeatMasker\ttandemRepeat\tclinicalAssociation\tdistanceToSplice\tmicroRNAs\tkeggPathway\tcpgIslands\ttfbs\tgenomesESP\tPPI\tproteinSequence\tdbSNP_ID\tdbSNP_5_percent_in_all\tdbSNP_5_percent_any\tdbSNP_Mutation\tESPfreqs\tMAFinESP\tdbSNP_5_percentany\tdbSNP_5_percent_in_any\tCNV\tA1_AF_AFR\tA1_AF_AMR\tA1_AF_EAS\tA1_AF_FIN\tA1_AF_NFE\tA1_AF_SAS\tA1_AF_OTH\tA1_AF_max\tA1_HOMOF_AFR\tA1_HOMOF_AMR\tA1_HOMOF_EAS\tA1_HOMOF_FIN\tA1_HOMOF_NFE\tA1_HOMOF_SAS\tA1_HOMOF_OTH\tA1_HOMOF_max\tA2_AF_AFR\tA2_AF_AMR\tA2_AF_EAS\tA2_AF_FIN\tA2_AF_NFE\tA2_AF_SAS\tA2_AF_OTH\tA2_AF_max\tA2_HOMOF_AFR\tA2_HOMOF_AMR\tA2_HOMOF_EAS\tA2_HOMOF_FIN\tA2_HOMOF_NFE\tA2_HOMOF_SAS\tA2_HOMOF_OTH\tA2_HOMOF_max\tF_rarest_allele\tCombo_max\tvarMD5\tgenomesExAC\n";
        print PREFS "INT NOT NULL\tVARCHAR(15) NOT NULL\tVARCHAR(12) NULL\tENUM('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT') NOT NULL\tINT NOT NULL\tVARCHAR(45) NULL\tVARCHAR(100) NULL\tVARCHAR(45) NULL\tVARCHAR(255) NULL\tVARCHAR(45) NULL\tVARCHAR(45) NULL\tVARCHAR(100) NULL\tINT NULL\tVARCHAR(15) NULL\tVARCHAR(20) NULL\tVARCHAR(8) NULL\tVARCHAR(25) NULL\tVARCHAR(10) NULL\tFLOAT(8,5) NULL\tVARCHAR(10) NULL\tVARCHAR(10) NULL\tVARCHAR(10) NULL\tVARCHAR(255) NOT NULL\tVARCHAR(10) NULL\tVARCHAR(10) NULL\tVARCHAR(10) NULL\tENUM('yes','no') NULL\tVARCHAR(180) NULL\tVARCHAR(30) NULL\tVARCHAR(200) NULL\tTEXT NULL\tINT NULL\tVARCHAR(30) NULL\tTEXT NULL\tVARCHAR(45) NULL\tVARCHAR(255) NULL\tVARCHAR(60) NULL\tTEXT NULL\tVARCHAR(45) NULL\tVARCHAR(30) NULL\tVARCHAR(25) NULL\tVARCHAR(25) NULL\tVARCHAR(25) NULL\tVARCHAR(50) NULL\tFLOAT NULL\tVARCHAR(25) NULL\tVARCHAR(25)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tVARCHAR(255)\tFLOAT\tVARCHAR(255)\tCHAR(32)\tVARCHAR(255)";
        print "Field\-type preference file generated\, may be edited as a tab-delimited text if needed\.\n";
        close PREFS;  #closes the preferences file
    }
    
    open (PREFS, "seattleSequel.prefs.txt") or die "Unable to open preference file\.\n"; #opens the preferences file
    my $prefkey = <PREFS>;  #reads the first line of the file into prefkey
    my $prefvalue = <PREFS>;  #reads the second line of the file into prefvalue
    close PREFS;  #closes the file
    chomp $prefkey;  #eliminates leading or trailing spaces or linebreaks from this value
    chomp $prefvalue;  #same as above
    my @prefkey = split(/\t/, $prefkey);  #splits the prefkey string at each tab to create an array
    my @prefvalue = split(/\t/, $prefvalue);  #splits the prefvalue string at each tab to create an array
    my %prefhash = ();  #initializes an empty hash of preference values
    my $index = 0;  #initializes index
    my $fieldtype = 0;  #initializes the user-chosen field type value to 0
    foreach my $key(@prefkey){  #iterates through each entry in the prefkey array
        $prefhash{$prefkey[$index]} = $prefvalue[$index];  #writes each entry to the preferences hash for the corresponding key and value pairs
        $index ++;  #increments index
    }
    open(OUTPUT, $seqout) or die "\nError opening output file $seqout\.\n";  #reopens the file used to write the SQL-ready sequence
    my $headers = <OUTPUT>;  #reads the header (first) line from the file
    close OUTPUT;  #closes the file
    chomp $headers;  #eliminates any leading or trailing space or line break characters from the line
    my @headers = split(/\t/, $headers);  #splits the string into an array
    foreach my $headervalue(@headers){  #goes through each header value in the headers
        if (!defined $prefhash{$headervalue}) {  #if it is undefined (becaue it was not in the preferences file)
            print "No default value found for $headervalue\.  Please select from one of the following\:\n";
            print "1  TEXT\n2  VARCHAR\(255\)\n3  INTEGER\n4  FLOAT\n";
            $fieldtype = readline STDIN;  #takes the user input
            chomp $fieldtype;  #eliminates leading and trailing linebreaks and other spacers
            while ($fieldtype !~ /^[1-4]$/) {  #keeps this loop going unless the fieldtype has been set to 1, 2, 3, or 4
                print "Invalid choice\, please select again\.\n";
                $fieldtype = readline STDIN;  #takes in another entry for the fieldtype
                chomp $fieldtype;
            }
            if ($fieldtype == 1) {  #if field type was one
                $prefhash{$headervalue} = "TEXT NULL";   #sets the matching preference hash value accordingly (next 3 statements do similar things)
            }
            if ($fieldtype == 2) {
                $prefhash{$headervalue} = "VARCHAR(255) NULL";
            }
            if ($fieldtype == 3) {
                $prefhash{$headervalue} = "INT NULL";
            }
            if ($fieldtype == 4) {
                $prefhash{$headervalue} = "FLOAT NULL";
            }
            
        }
        
    }
    print "Generating mySQL script\.\n";
    open(SQL, '>'.$seqout."\.sql") or die "\nError opening output file $seqout\.sql\.\n";  #opens the file that will contain the SQL script
    print SQL "USE $dbname \;\n";  #writes the first line of the file that says which database schema to use
    print SQL "CREATE TABLE $dbname\.$tablename \(\n";  #writes the next line to tell it to generate the table
    my $hlindex = 0; #initializes an index for header position to look for the end
    foreach my $headervalue(@headers) {  #iterates through the array of headers
        print SQL "$headervalue $prefhash{$headervalue}";  #putting out a line of code for each one to create the appropriate field
        if ($headers[$hlindex+1]) {  #checks to see if the next header value is empty (indicating we have reached the end of the headers)
            print SQL "\,\n";
        }
        else {
            print SQL "\)\;\n";
        }
        $hlindex ++;
    }
    print SQL "LOAD DATA LOCAL INFILE \'$seqout\' INTO TABLE $tablename\nCOLUMNS TERMINATED BY \'\\t\'\nLINES TERMINATED BY \'\\n\'\nIGNORE 1 lines\;\n";  #writes the final lines of the SQL script that tell it where to get the data and how to read it
    print SQL "ALTER TABLE $tablename ADD INDEX \(varMD5\)\, ADD INDEX \(uniLoc\)\, ADD INDEX \(geneList\)\;\n";  #creates an index by locus, gene, and line number
    close SQL;  #closes the SQL file
    print "Done\!\n";  #prints that it is done
}

sub checkopts{
    my $opts = shift;  #dereferences the hash containing the options
    
    my $file = $opts->{"f"}; #puts the value in options under key F into a variable called file
    
    unless(defined($file) and (-e $file)){  #unless the file entered exists...
        print "Input file not found or not defined in commandline arguments.\n";
        return 0;  #this function will return a value of 0, which is false and signals bad options
    }
}

sub usage{  #This subroutine prints directions
    print "This program will prepare a SeattleSeq output for loading into mySQL by cleaning up the headers, adding a few convenient columns, and generating the script to load the charts\.\nPlease use full \(not relative\) filenames to ensure the SQL script executes properly\.\nSample commandline\:\nperl seattleSequel\.pl \-f file\.txt\.\n";
    die "";
}


main();
