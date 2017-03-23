#!/usr/bin/perl
$num_args = $#ARGV + 1;

if ($num_args != 1) {
   print "\n Verwendung: Bamcount_filter.pl critical_p-value\n";
    exit;
}
mkdir bamcount_filtered;

while ($file1=<*.table>){
$file1=~/^(.+).table/ or die "strange file format: $file1\n";
$bamfile=$1;
print "\n$file1\n$bamfile\n\n";



#Abfrage, ob DBSCAN bereits existiert, bzw Isolat schon berechnet wurde
#  while ($dir=<${out}/${bamfile}*.table>){
#  print "$dir\n";
#  $h{$1}=1;
#  }
#  unless(defined $h{$1}){ #wenn das nicht der Fall ist geht es los:

open(Fout,">bamcount_filtered/${bamfile}_${ARGV[0]}-filtered.table") or die "\n\n\"Offnen von $ARGV[0] nicht m\"oglich\n\n\n";
open(Fin,"<$file1") or die "\n\n\"Offnen von $file1 nicht m\"oglich\n\n\n";

print Fout "POS\tREF\tALT\tDP\t#ALT\tFreq(ALT)\tQual(ALT)\tp-value\tantibiotic\n";

$line=<Fin>;
    while($line=<Fin>){
    @array=split("\t", $line);
        if ($array[18]<$ARGV[0]){#Abfrage, ob der p-Wert kleiner dem zu filternden p-Wert ist
            if ($array[2] eq "A"){
            print Fout "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[16]\t$array[17]\t$array[18]\t$array[19]";}
            elsif($array[2] eq "C"){
            print Fout "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[7]\t$array[16]\t$array[17]\t$array[18]\t$array[19]";}
            elsif($array[2] eq "G"){
            print Fout "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[10]\t$array[16]\t$array[17]\t$array[18]\t$array[19]";}
            elsif($array[2] eq "T"){
            print Fout "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[13]\t$array[16]\t$array[17]\t$array[18]\t$array[19]";}
            else {print "$line\n";}
        }
    }   
}
