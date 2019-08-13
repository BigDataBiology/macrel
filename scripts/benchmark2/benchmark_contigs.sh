#time ~/Desktop/Lib/FACS16_conda_envs.sh -m c --fasta representative_genomes/repcont.fna.gz --outfolder analysis/ --outtag HMRP_contigs -t 3 --block 100M --log HMRP_contigs.log --mem 0.75

time ~/Desktop/Lib/FACS16_conda_envs.sh -m r --fwd metagenome/SRR172902_1.fastq.gz --outfolder analysis/ --outtag SRR172902 -t 3 --block 100M --log SRR172902.log --mem 0.75

time ~/Desktop/Lib/FACS16_conda_envs.sh -m r --fwd metagenome/SRR172903_1.fastq.gz --outfolder analysis/ --outtag SRR172903 -t 3 --block 100M --log SRR172903.log --mem 0.75

