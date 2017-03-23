# Skripte_vschleusener
Perl-Skripte der Doktorarbeit

In diesem Repository findet sich eine Auswahl der wichtigsten Perl-Skripte, die während meiner Dissertation am Forschungszentrum Borstel entstanden sind.

Vorraussetzungen zum Ausführen sind R (https://cran.r-project.org/), Perl (https://www.perl.org/get.html), das Perlmodul 'Statistics::R', das GenomeAnalysisToolkit (https://software.broadinstitute.org/gatk/) und der Bam-readcount-Algorithmus (https://github.com/genome/bam-readcount).

Sind die Vorraussetzungen erfüllt können die Programme mit der Kommandozeile im Ordner der BAM-Dateien aufgerufen werden. Dabei werden jeweils alle BAM-Dateien prozessiert, die sich im Ordner befinden, sofern die Berechnung nicht bereits durchgeführt wurde.

Es wird empfohlen vorprozessierte BAM-Dateien zu verwenden, das heißt Alignments, in denen Duplikate entfernt, Artefakte um Indels repariert und die Basenqulität rekalibriert wurden. Zusätzlich ist die Interpretation der Ergebnisse bei einer mittleren Coverage < 50x unsicher.

Vor der ersten Anwednung müssen die Pfade zu den Programmen wie GATK oder zur Referenzsequenz in den Skripte angepasst werden.

 Verwendung:
 Mischdetektor --> Batch_mixDetector.pl /Pfad/zur/Ausgabe/Dateiname.csv /Ausgabe/Pfad
 
 In der CSV-Ausgabedatei befindet sich die Entscheidung für jede BAM-Datei (Mischung Ja/Nein), sowie die detektierten Cluster. Zusätzlich wird für jede Datei ein Scatterplot des DBSCAN-Algorithmus erstellt und im Ausgabepfad gespeichert.
 
 Binomialtest -->  Batch_BamBasecount_pvalue.pl /Pfad/zur/interval.list /Pfad/zur/RefuAlt.table
 
 In diesem Repository sind zwei mögliche Intervalllisten beziehungsweise RefuAlt-Tabellen enthalten. Einmal zur Bestimmung der Linienspezifischen SNPs (Phylosnps_Master.v28) und einmal zur Bestimmung von Varianten in resistenzvermittelnden Genen (Resisnps_Master.v28).
 Die Bestimmung der phylogenetischen Marker bietet sich im Anschluss zum Mischungsdetektor an, um die enthaltenen Genotypen feststellen zu können.
 Für den Mischungsdetektor ist eine Intervallliste für das Kerngenom inkludiert, sodass nur Varianten innerhalb des Kerngenoms in die Analyse einbezogen werden.
