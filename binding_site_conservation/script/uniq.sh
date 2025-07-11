cut -f 6 HUMAN_ARATH.tsv | sort | uniq > HUMAN_ARATH.uniq.txt
cut -f 6 HUMAN_MOUSE.tsv | sort | uniq > HUMAN_MOUSE.uniq.txt
cut -f 6 HUMAN_RAT.tsv | sort | uniq > HUMAN_RAT.uniq.txt
cut -f 6 HUMAN_XENLA.tsv | sort | uniq > HUMAN_XENLA.uniq.txt
cut -f 6 HUMAN_PONAB.tsv | sort | uniq > HUMAN_PONAB.uniq.txt
cut -f 6 HUMAN_DANRE.tsv | sort | uniq > HUMAN_DANRE.uniq.txt
cut -f 6 HUMAN_DROME.tsv | sort | uniq > HUMAN_DROME.uniq.txt
cut -f 6 HUMAN_MACFA.tsv | sort | uniq > HUMAN_MACFA.uniq.txt
cut -f 6 HUMAN_YEAST.tsv | sort | uniq > HUMAN_YEAST.uniq.txt
cut -f 6 HUMAN_CHICK.tsv | sort | uniq > HUMAN_CHICK.uniq.txt
#cat *txt | sort | uniq -c |  cut -d "E" -f 1 | sort | uniq -c | less
#cat *txt | sort | uniq -c |  grep "10 E" > core.txt
cut -f 6 HUMAN_ARATH.ncrna.tsv | sort | uniq > HUMAN_ARATH.ncrna.uniq.txt
cut -f 6 HUMAN_MOUSE.ncrna.tsv | sort | uniq > HUMAN_MOUSE.ncrna.uniq.txt
cut -f 6 HUMAN_RAT.ncrna.tsv | sort | uniq > HUMAN_RAT.ncrna.uniq.txt
cut -f 6 HUMAN_XENLA.ncrna.tsv | sort | uniq > HUMAN_XENLA.ncrna.uniq.txt
cut -f 6 HUMAN_PONAB.ncrna.tsv | sort | uniq > HUMAN_PONAB.ncrna.uniq.txt
cut -f 6 HUMAN_DANRE.ncrna.tsv | sort | uniq > HUMAN_DANRE.ncrna.uniq.txt
cut -f 6 HUMAN_DROME.ncrna.tsv | sort | uniq > HUMAN_DROME.ncrna.uniq.txt
cut -f 6 HUMAN_MACFA.ncrna.tsv | sort | uniq > HUMAN_MACFA.ncrna.uniq.txt
cut -f 6 HUMAN_YEAST.ncrna.tsv | sort | uniq > HUMAN_YEAST.ncrna.uniq.txt
cut -f 6 HUMAN_CHICK.ncrna.tsv | sort | uniq > HUMAN_CHICK.ncrna.uniq.txt
#cat *.ncrna.uniq.txt | sort | uniq -c |  cut -d "E" -f 1 | sort | uniq -c | less
#cat *.ncrna.uniq.txt | sort | uniq -c |  grep "10 E" >> core.txt