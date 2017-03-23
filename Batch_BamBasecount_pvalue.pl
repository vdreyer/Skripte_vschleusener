#!/usr/bin/perl
#Skript muss im Ordner der BAM-Dateien aufgerufen werden
mkdir bamcount_pvalue; #eventuell neuer Ordner für die Ergebnisse
$num_args = $#ARGV + 1;
if ($num_args != 2) {
   print "\n Verwendung: Batch_BamBasecount_pvalue.pl /Pfad/zur/interval.list /Pfad/zur/RefuAlt.table \n";
    exit;
#interval.list Datei:
#M.tuberculosis_H37Rv\tstart\tstop1
#M.tuberculosis_H37Rv\tstart2\tstop2
#M.tuberculosis_H37Rv\tstart3\tstop3
#M.tuberculosis_H37Rv\tstart4\tstop4

#RefuAlt.table-Datei:
#Pos\tREF\tALT\Antibiotic

}
#Referenzdatei benennen; Pfad muss angepasst werden
$h37='/auto/Thecus_Analysis/SeqWork/DATA/M._tuberculosis_H37Rv_2013-02-15';
#Schleife über alle BAM-Dateien des Ordners
while($file1=<*.bam>){
$file1=~/^(.+).bam/ or die "Format unkorrekt: $file1\n";
$bamfile=$1;
print "\n$file1\n$bamfile\n\n";

while ($dir=<bamcount_pvalue/${bamfile}*.table>){  #Abfrage, ob die Datei bereits berechnet wurde
  print "$dir\n";
  $h{$1}=1;
  }
  unless(defined $h{$1}){  #wenn nicht, geht es los

#bam-readcount aufrufen; Pfad zum Algorithmus muss angepasst werden
$do="/home/vschleusener/Downloads/bam-readcount-master/bin/bam-readcount $file1 -b 15 -f $h37.fasta -l $ARGV[0] > bamcount_pvalue/$bamfile.txt ";
print "\n\n\n".$do."\n\n";
system($do); if ($?){die "$do did not work: $?\n";}

#bam-readcount Ergebnisse in Tabelle konvertieren , Pfad muss angepasst werden
$do="perl /home/vschleusener/Skripte/readcount2table_binomial.pl bamcount_pvalue/$bamfile.table bamcount_pvalue/$bamfile.txt $ARGV[1]";
print "\n\n\n".$do."\n\n";
system($do); if ($?){die "$do did not work: $?\n";}
}
}