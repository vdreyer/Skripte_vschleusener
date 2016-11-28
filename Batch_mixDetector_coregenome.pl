#!/usr/bin/perl
use Statistics::R;
# Programm verlassen, solange nicht die korrekte Zahl an Uebergabeparameter vorhanden sind
$num_args = $#ARGV + 1;
if ($num_args != 2) {
   print "\n Verwendung: Batch_mixDetector.pl /path/to/output/filename.csv /output/path\n";
    exit;
}
open(Fout,">>$ARGV[0]") or die "\n\ncannot open $ARGV[0]\n\n\n"; #Ausgabedatei erstellen/Oeffnen falls schon vorhanden
while ($file1=<*.bam>){
 $file1=~/^(.+).bam/ or die "strange file format: $file1\n";
 $bamfile=$1;
 print "\n$file1\n$bamfile\n\n";
 $out=$ARGV[1];#Ausgabepfad aus Uebergabeparametern verwenden

#Abfrage, ob DBSCAN bereits existiert, bzw Isolat schon berechnet wurde
  while ($dir=<${out}/${bamfile}*.table>){
   print "$dir\n";
   $h{$1}=1;
   }
   unless(defined $h{$1}){ #wenn das nicht der Fall ist geht es los:

   $dir=$out;
   $h37='/auto/Thecus_Analysis/SeqWork/DATA/M._tuberculosis_H37Rv_2013-02-15'; #Pfad zur Referenzdatei mit Index-Dateien etc

   #Varianten fuer das Kerngenom detektieren aus der eingelesenen BAM-Datei
   $do="java -jar /home/vschleusener/Downloads/GenomeAnalysisTK_new/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${h37}.fasta -I $file1 -o ${dir}/${bamfile}_Varianten_GATK.vcf -dcov 10000 -L /home/vschleusener/Dokumente/Core-Genome.interval_list -glm SNP -mbq 13 -A BaseCounts -A Coverage -nct 4 -nt 4 -rf BadCigar";
   print "\n\n\n".$do."\n\n";
   system($do); if ($?){die "$do did not work: $?\n";}

   #--> Varianten in einen Tabelle konvertieren fuer R
   $do="java -jar /home/vschleusener/Downloads/GenomeAnalysisTK_new/GenomeAnalysisTK.jar -T VariantsToTable -R ${h37}.fasta -V ${dir}/${bamfile}_Varianten_GATK.vcf -F POS -F REF -F ALT -F BaseCounts -dcov 10000 -F DP --allowMissingData --splitMultiAllelic -o ${dir}/${bamfile}_Varianten.table";
   print "\n\n\n".$do."\n\n";
   system($do); if ($?){die "$do did not work: $?\n";}

   #BaseCounts ebenfalls als eigene Spalten machen und " " als seperator verwenden --> mehr Spalten als Ueberschriften
   $do="awk 'BEGIN{FS=\"[\\t]|[,]|[ []|[] ]\";OFS=\" \"};{\$1=\$1};{print \$0}' < ${dir}/${bamfile}_Varianten.table  > ${dir}/${bamfile}_Varianten2.table
   awk 'BEGIN{FS=\"  \";OFS=\" \"};{\$1=\$1};{print \$0}' < ${dir}/${bamfile}_Varianten2.table  > ${dir}/${bamfile}_Varianten.table
   rm ${dir}/${bamfile}_Varianten2.table";  #zweimal um eventuell splitMultiAllelic huebsch zu formatieren
   print "\n\n\n".$do."\n\n";
   system($do); if ($?){die "$do did not work: $?\n";}

   #R-Code
   my $R = Statistics::R->new();
   $R->startR ;
   $R->send('library("fpc", lib.loc="/usr/lib64/R/library")');#R-Library fuer DBSCAN
   print "Tabelle einlesen...\n";
   $R->run(qq'mix<-read.table(file="${dir}/${bamfile}_Varianten.table", sep=" ", header=F, skip=1);');#Tabelle mit Basecalls einlesen und Ueberschriften auslassen, da weniger Ueberschriften als Spalten
   print "Spalten umbenennen...\n";
   $R->run(q'colnames(mix)=c("POS", "REF", "ALT", "A", "C", "G", "T", "DP")'); #Spalten benennen
   $R->run(qq'write.table(mix, "${dir}/${bamfile}_Varianten.table",quote=FALSE, sep = " ", row.names=F)');
   print "Computing ratios...\n";
   # neue Spalten Anteil.Ref und Anteil.Alt berechnen
   $R->run('for (i in 1:dim(mix)[1]){ifelse(mix$REF[i] == "A", mix$Anteil.REF[i]<-mix$A[i]/mix$DP[i],ifelse(mix$REF[i] == "C",mix$Anteil.REF[i]<-mix$C[i]/mix$DP[i],ifelse(mix$REF[i] == "G",mix$Anteil.REF[i]<-mix$G[i]/mix$DP[i],ifelse(mix$REF[i] == "T"
   ,mix$Anteil.REF[i]<-mix$T[i]/mix$DP[i], mix$Anteil.REF[i]<-NA))))}
   for (i in 1:dim(mix)[1]){ifelse(mix$ALT[i] == "A", mix$Anteil.ALT[i]<-mix$A[i]/mix$DP[i],ifelse(mix$ALT[i] == "C",mix$Anteil.ALT[i]<-mix$C[i]/mix$DP[i],ifelse(mix$ALT[i] == "G",mix$Anteil.ALT[i]<-mix$G[i]/mix$DP[i],ifelse(mix$ALT[i] == "T",mix$Anteil.ALT[i]<-mix$T[i]/mix$DP[i], mix$Anteil.ALT[i]<-NA))))}');
   print "Computing DBSCAN with ratios...\n";
   $R->run('m<-cbind(mix$Anteil.ALT,mix$Anteil.REF)
   colnames(m) <- c("alt", "ref")
   m<-m[complete.cases(m),]
   db<-dbscan(m,0.06,55)'); #Matrix aus den Anteilen bilden und DBSCAN aufrufen
   print "Computing Cluster center...\n";
   $R->run('center<-rbind(colMeans(m[db$cluster==1, ]),colMeans(m[db$cluster==2, ]),colMeans(m[db$cluster==3, ]),colMeans(m[db$cluster==4, ]))');
   print "Plot graph and color cluster...\nSave plot...\n";
   $R->run(qq'png("${dir}/${bamfile}_DBSCAN.png")');#Grafik speichern
   $R->run('plot(mix$Anteil.ALT,mix$Anteil.REF, xlab="alternative allele", ylab="reference allele", main="Scatterplot of the ratio of alt and ref alleles")
   points(m, col=db$cluster)
   points(rbind(colMeans(m[db$cluster==1,]),colMeans(m[db$cluster==2,]),colMeans(m[db$cluster==3,]),colMeans(m[db$cluster==4,])),col="yellow", pch=8, cex=2)
   dev.off()');#DBSCAN Clusterzentren berechnen und alles zeichnen
   print "\nMixture detected: ";
   $R->run('if (max(db$cluster)>1){ismix<-"YES"
   anteil<-center[,1]} else {ismix<-"NO"
   anteil<-center[,1]}');
   $ismix=$R->get('ismix');
   $anteil=$R->get('anteil');
   @anteil2= grep {$_ ne 'NaN'} @$anteil;
   print "$ismix\n\nCluster center: @anteil2\n";
   print Fout "$bamfile\t$ismix\t@anteil2\n";  #in Ausgabetabelle schreiben
   $R->stopR() ;
 }
}