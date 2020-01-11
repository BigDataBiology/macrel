#!/usr/bin/env bash


workDir=$( cd $(dirname $0); pwd )

outfolder_mode_p="out_peptides"
outtag_mode_p="example"

outfolder_mode_c="out_contigs"
outtag_mode_c="out_contigs"

outfolder_mode_r="out_metag"
outtag_mode_r="example_metag"

outfolder_mode_a_fasta="out_abundance"
outtag_mode_a_fasta="example_abundance"

outfolder_mode_a_ref="out_abundance_pred"
outtag_mode_a_ref="example_abundance_pred"


cd $HOME/FACS/
echo ''


echo "=================== test p mode with example_seqs/expep.faa.gz ==================="
echo "***** run FACS.sh in p mode *****"
echo ''

bash FACS.sh --mode p --fasta example_seqs/expep.faa.gz --outfolder $outfolder_mode_p --outtag $outtag_mode_p -t 4
echo ''

echo "***** finish running *****"
echo ''

echo "***** validating results *****"
mode_p_result1=$outfolder_mode_p/$outtag_mode_p.tsv.gz
mode_p_result1_unzip=$outfolder_mode_p/$outtag_mode_p.tsv
mode_p_result2=$outfolder_mode_p/$outtag_mode_p.ids.tsv.gz
mode_p_result2_unzip=$outfolder_mode_p/$outtag_mode_p.ids.tsv
mode_p_expect1=$workDir/expect.mode_p.results/expect.tsv

if [ -f $mode_p_result1 ] && [ -f $mode_p_result2 ];
then
	zcat $mode_p_result1 > $mode_p_result1_unzip
	diff $mode_p_result1_unzip $mode_p_expect1
	if [ $? -eq 0 ]
	then
		echo "the results are same with what we expect."
		echo "test passed."
		rm -rf $mode_p_result1_unzip
	else
		rm -rf $mode_p_result1_unzip
		echo "The results are different from what we expect.please check them."
		echo "test failed."
		echo "exit!"
		echo ''
		exit 1
	fi
else
	echo "sorry,$mode_p_result1 or $mode_p_result2 does not exist."
	echo "the test failed to get right results,"
	echo "exit!"
	echo ''
	exit 1
fi
echo "***** finish validating *****"

echo "=================== finish to test p mode with example_seqs/expep.faa.gz ==================="
echo ''


echo "=================== test c mode with example_seqs/excontigs.fna.gz ==================="
echo "***** run FACS.sh in c mode *****"
echo ''

bash FACS.sh --mode c --fasta example_seqs/excontigs.fna.gz --outfolder $outfolder_mode_c --outtag $outtag_mode_c -t 4
echo ''

echo "***** finish running *****"
echo ''

echo "***** validating results *****"
mode_c_result1=$outfolder_mode_c/$outtag_mode_c.tsv.gz
mode_c_result1_unzip=$outfolder_mode_c/$outtag_mode_c.tsv
mode_c_result2=$outfolder_mode_c/$outtag_mode_c.ids.tsv.gz
mode_c_result2_unzip=$outfolder_mode_c/$outtag_mode_c.ids.tsv
mode_c_expect1=$workDir/expect.mode_c.results/expect.tsv

if [ -f $mode_c_result1 ] && [ -f $mode_c_result2 ];
then
	zcat $mode_c_result1 > $mode_c_result1_unzip
	diff $mode_c_result1_unzip $mode_c_expect1
	if [ $? -eq 0 ]
	then
		echo "the results are same with what we expect."
		echo "test passed."
		rm -rf $mode_c_result1_unzip
	else
		rm -rf $mode_c_result1_unzip
		echo "The results are different from what we expect.please check them."
		echo "test failed."
		echo "exit!"
		echo ''
		exit 1
	fi
else
	echo "sorry,$mode_c_result1 or $mode_c_result2 does not exist."
	echo "the test failed to get right results,"
	echo "exit!"
	echo ''
	exit 1
fi
echo "***** finish validating *****"

echo "=================== finish to test c mode with example_seqs/excontigs.fna.gz ==================="
echo ''


echo "=================== test r mode with example_seqs/R1.fq.gz and example_seqs/R2.fq.gz ==================="
echo "***** run FACS.sh in r mode *****"
echo ''

bash FACS.sh --mode r --fwd example_seqs/R1.fq.gz --rev example_seqs/R2.fq.gz --outfolder $outfolder_mode_r --outtag $outtag_mode_r -t 4
echo ''

echo "***** finish running *****"
echo ''

echo "***** validating results *****"
mode_r_result1=$outfolder_mode_r/$outtag_mode_r.tsv.gz
mode_r_result1_unzip=$outfolder_mode_r/$outtag_mode_r.tsv
mode_r_result2=$outfolder_mode_r/$outtag_mode_r.ids.tsv.gz
mode_r_result2_unzip=$outfolder_mode_r/$outtag_mode_r.ids.tsv
mode_r_expect1=$workDir/expect.mode_r.results/expect.tsv

if [ -f $mode_r_result1 ] && [ -f $mode_r_result2 ];
then
	zcat $mode_r_result1 > $mode_r_result1_unzip
	diff $mode_r_result1_unzip $mode_r_expect1
	if [ $? -eq 0 ]
	then
		echo "the results are same with what we expect."
		echo "test passed."
		rm -rf $mode_r_result1_unzip
	else
		rm -rf $mode_r_result1_unzip
		echo "The results are different from what we expect.please check them."
		echo "test failed."
		echo "exit!"
		echo ''
		exit 1
	fi
else
	echo "sorry,$mode_r_result1 or $mode_r_result2 does not exist."
	echo "the test failed to get right results,"
	echo "exit!"
	echo ''
	exit 1
fi
echo "***** finish validating *****"

echo "=================== finish to test r mode with example_seqs/R1.fq.gz and example_seqs/R2.fq.gz ==================="
echo ''


echo "=================== test a mode with example_seqs/R1.fq.gz and example_seqs/ref.faa.gz ==================="
echo "***** run FACS.sh in a mode *****"
echo ''

bash FACS.sh --mode a --fwd example_seqs/R1.fq.gz --fasta example_seqs/ref.faa.gz --outfolder $outfolder_mode_a_fasta --outtag $outtag_mode_a_fasta -t 4
echo ''

echo "***** finish running *****"
echo ''

echo "***** validating results *****"
mode_a_fasta_result1=$outfolder_mode_a_fasta/$outtag_mode_a_fasta.xprs
mode_a_fasta_result2=$outfolder_mode_a_fasta/$outtag_mode_a_fasta.params.xprs

if [ -e $mode_a_fasta_result1 ] && [ -e $mode_a_fasta_result2 ];
then
	echo "test passed."
else
	echo "sorry,$mode_a_fasta_result1 or $mode_a_fasta_result2 does not exist."
	echo "the test failed to get results,"
	echo "exit!"
	echo ''
	exit 1
fi
echo "***** finish validating *****"

echo "=================== finish to test a mode with example_seqs/R1.fq.gz and example_seqs/ref.faa.gz ==================="
echo ''


echo "=================== test mode a with example_seqs/R1.fq.gz and result of mode r ==================="
echo "***** run FACS.sh in a mode *****"
echo ''

bash FACS.sh --mode a --fwd example_seqs/R1.fq.gz --ref $outfolder_mode_r/$outtag_mode_r.tsv.gz --outfolder $outfolder_mode_a_ref --outtag $outtag_mode_a_ref -t 4
echo ''

echo "***** finish running *****"
echo ''

echo "***** validating results *****"
mode_a_ref_result1=$outfolder_mode_a_ref/$outtag_mode_a_ref.xprs
mode_a_ref_result2=$outfolder_mode_a_ref/$outtag_mode_a_ref.params.xprs

if [ -e $mode_a_ref_result1 ] && [ -e $mode_a_ref_result2 ];
then
	echo "test passed."
else
	echo "sorry,$mode_a_ref_result1 or $mode_a_ref_result2 does not exist."
	echo "the test failed to get results,"
	echo "exit!"
	echo ''
	exit 1
fi
echo "***** finish validating *****"

echo "=================== finish to test mode a with example_seqs/R1.fq.gz and result of mode r ==================="
echo ''

echo "========================================== done =========================================="
cd $workDir
echo ''
