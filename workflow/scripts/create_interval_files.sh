#!/bin/bash

# Create interval files from target/bait bed files

ml picard
ml java

RAW_DIR=/data/khanlab3/gb_omics/projects/khan/ccrgb_106/raw/NCI_RhadoSeqFoc_TE-95039494_hg38/
PANEL_DIR=/data/khanlab3/gb_omics/projects/khan/ccrgb_106/analysis/interval_files/

BAIT_BED=Probes_merged_ok_shareable_NCI_RhadoSeqFoc_TE-95039494_hg38_251209142544.bed
TARGET_BED=All_target_regions_NCI_RhadoSeqFoc_TE-95039494_hg38_251209142544.bed

DICT=/fdb/genomebrowser/fasta/hg38.dict

java -Xmx4g -jar $PICARDJARPATH/picard.jar BedToIntervalList \
  I=$RAW_DIR$BAIT_BED \
  O="$PANEL_DIR/NCI_RhadoSeqFoc_TE.baits.interval_list" \
  SD="$DICT"

java -Xmx4g -jar $PICARDJARPATH/picard.jar BedToIntervalList \
  I=$RAW_DIR$TARGET_BED \
  O="$PANEL_DIR/NCI_RhadoSeqFoc_TE.targets.interval_list" \
  SD="$DICT"
