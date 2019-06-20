#!/usr/bin/env python
__author__ = 'Shlomo Shenzis'

import gzip
import os
import sys
import argparse
from argparse import RawDescriptionHelpFormatter
from operator import add
import bisect
import copy
import subprocess
import re

version = " ______________________________\n" \
          "|         SRCP v1.2.0          |\n" \
          "|______________________________|\n" \
          "| Short Reads circRNA Pipeline |\n" \
          "|______________________________|\n" \
          "|  Written by: Shlomo Shenzis  |\n" \
          "|      @ Kadener Lab 2015      |\n" \
          "|______________________________|"
# Magic numbers:
CHR = 0
START = 0
INDEX = 0
SIGNATURE_IN_SAM = 0
LOCA = 1
GENOME = 1
STRAND_IN_SAM = 1
LOCB = 2
EXACT_SCORE = 2
NAME_IN_SAM = 2
TRANSCRIPTOME = 2
GENE_NAME = 3
SERIAL_NUMBER = 3
ALIGN_START_IN_SAM = 3
NUM_READS = 4
STRAND = 5
CIGAR_IN_SAM = 5
NUM_EXONS = 6
BED6_NUM_COL = 6
RUN_LINEAR_SCORE = 6
TRANSCRIPT_NAME = 7
CIRCLE_EXONS_STARTS = 7
FIELDS_PER_SAMPLE = 7
RUN_LINEAR_ALTERNATIVE_TRANSCRIPTS = 7
CIRCLE_EXONS_ENDS = 8
SCORE = 9
SEQ_IN_SAM = 9
QUAL_IN_SAM = 10
CHOSEN_TRANSCRIPT = 10
RUN_LINEAR_TRANSCRIPT = 11
ALTERNATIVE_TRANSCRIPTS = 11
SCORE_IN_EXONS = 12
EXONS_STARTS = 13
EXONS_ENDS = 14
REVERSE_IN_SAM = 16
EDIT_DISTANCE_IN_SAM = 16
RUN_LINEAR_EXON_STARTS = 17
RUN_LINEAR_EXON_ENDS = 18
# For sequence comparison
NO_MATCH = -1
MATCH_SENSE = 0
MATCH_ANTI_SENSE = 16

COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
CIGAR_FLAGS = ['N', 'H', 'P', '=', 'X'] #avigayel: removed I and S to allow soft clipping and insertion
EDIT_DISTANCE_REGEX = re.compile(ur'NM:i:([0-9]+)', re.IGNORECASE)

#Shahar_V7
#Input as bed 12 annotatin
def bed_12_to_11():
    BED_START = 1
    BED_EXON_SIZES = 10
    BED_EXON_STARTS = 11
    BED_NUM_OF_EXONS = 9
    bed12annot = open(args.annotation, 'r').read().strip('\n').split('\n')
    bed11annot = open('tmp/BED11_annot.bed', 'w')
    for row in bed12annot:
        row = row.strip('\t').split('\t')

        # Calculating Exons start and end by the relative position and block-size received as input
        new_start = list()
        new_end = list()
        old_start = row[BED_EXON_STARTS][:-1].split(",")
        old_sizes = row[BED_EXON_SIZES][:-1].split(",")
        for i in xrange(int(row[BED_NUM_OF_EXONS])):
            new_start.append(str(int(row[BED_START]) + int(old_start[i])))
            new_end.append(str(int(new_start[i]) + int(old_sizes[i])))

        new_start_str = ",".join(new_start) + ","
        new_end_str = ",".join(new_end) + ","
        row = row[:8] + [row[9]] + [new_start_str] + [new_end_str]
        bed11annot.write("\t".join(map(str, row)) + '\n')
    bed11annot.close()
    args.annotation = 'tmp/BED11_annot.bed'

def sort_filtered_circs(filtered_circs_lst):

    annotation_file = open(args.annotation,'r')
    annotation = annotation_file.read().strip('\n').split('\n')
    chr_lst=[] #list of chromosomes annotation
                     # each element is a dictionary,
                     # each dictionary has all the annotated genes
                     # the key is the start position and the value is the end position

    #gene_index = chr_dict[chr_index].keys().index(key)
    chr_dict=[]
    chr=''
    for line in annotation:
       line = line.split('\t')
       chr=line[CHR]
       if chr not in chr_lst:
           #add the chromosome to the list of chromosomes
           chr_lst.append(chr)
           #initial the dictionary
           chr_dict.append({})

       #add the annotation to the correct dictionary
       chr_index = chr_lst.index(chr)
       key = line[1]
       value = line[2]
       chr_dict[chr_index][key] = value
    annotation_file.close()

    sorted_dicts=[]
    sorted_lst=[]
    for i in range(0,len(chr_dict)):
        sorted_lst = sorted(chr_dict[i].keys(),key=int)
        sorted_lst = [int(x) for x in sorted_lst]
        sorted_dicts.append(sorted_lst)

     #check three cases:
     ##1) one end of the circ is in one gene and the other end is in a different gene
     ##2) one end of the circ is in the gene and the other is not in any gene
     ##3) both ends of the circ are not in any gene

    partialy_not_annotated_circ=[]
    trans_splicing_circs=[]
    fully_not_annotated_circs=[]

    for circ in filtered_circs_lst:
       #circ = circ.split('\t')
       chr=circ[CHR]
       start_circ=circ[LOCA]
       end_circ= circ[LOCB]
       #look for the location of the circ in the genes in the chromosome annotations
       if chr not in chr_lst:
           fully_not_annotated_circs.append(circ)
           continue
       chr_index = chr_lst.index(chr)
       circ_location_in_dict= bisect.bisect_right(sorted_dicts[chr_index],int(start_circ))


       previous_gene_start=0
       if(circ_location_in_dict>0):
           previous_gene_start = sorted_dicts[chr_index][circ_location_in_dict-1]
           previous_gene_end = int(chr_dict[chr_index][str(previous_gene_start)])
       else:
           previous_gene_end = 'none'

       next_gene_start = len(sorted_dicts[chr_index])
       if circ_location_in_dict < len(sorted_dicts[chr_index]):
           next_gene_start = sorted_dicts[chr_index][circ_location_in_dict]
           next_gene_end = int(chr_dict[chr_index][str(next_gene_start)])
       else:
           next_gene_end = 'none'

        #check which case are we having
       start_circ = int(start_circ)
       end_circ = int(end_circ)
       if previous_gene_end == 'none' or next_gene_end== 'none':
           if previous_gene_end == 'none':
               if end_circ < next_gene_start:
                   fully_not_annotated_circs.append(circ)
               elif end_circ >= next_gene_start:
                   partialy_not_annotated_circ.append(circ)
           elif next_gene_end == 'none':
               if start_circ > previous_gene_end:
                   fully_not_annotated_circs.append(circ)
               elif start_circ<= previous_gene_end:
                   partialy_not_annotated_circ.append(circ)
       else:
           if start_circ > previous_gene_end and end_circ < next_gene_start:
               fully_not_annotated_circs.append(circ)
           elif start_circ<= previous_gene_end and end_circ >= next_gene_start:
               trans_splicing_circs.append(circ)
           elif start_circ <= previous_gene_end and end_circ <= next_gene_start:
               partialy_not_annotated_circ.append(circ)
           elif start_circ >= previous_gene_end and end_circ >= next_gene_start:
               partialy_not_annotated_circ.append(circ)


    output = open('tmp/exons.bed', 'a')
    #fully_out=open('tmp/fully_not_annotated','w')
    #part_out=open('tmp/partially_not_annotated','w')
    #trans_out=open('tmp/trans_splicing','w')

    for item in fully_not_annotated_circs:
        item[SERIAL_NUMBER] = item[SERIAL_NUMBER]+'_fully_not_annotated'
        output.write(item[CHR] +'\t' +str(item[LOCA])+'\t'+str(item[LOCB])+'\t'+item[SERIAL_NUMBER] +'\t'+'1000\t'+item[STRAND]+  '\t'+ str(item[LOCA])+'\t'+str(item[LOCB])+'\t0,0,255\t1\t'+str(int(item[LOCB])-int(item[LOCA]))+'\t0\t0\t0\n')
    for item in partialy_not_annotated_circ:
        item[SERIAL_NUMBER] = item[SERIAL_NUMBER]+'_partialy_not_annotated'
        output.write(item[CHR] +'\t' +str(item[LOCA])+'\t'+str(item[LOCB])+'\t'+item[SERIAL_NUMBER] +'\t'+'1000\t'+item[STRAND]+  '\t'+ str(item[LOCA])+'\t'+str(item[LOCB])+'\t0,0,255\t1\t'+str(int(item[LOCB])-int(item[LOCA]))+'\t0\t0\t0\n')
    for item in trans_splicing_circs:
        item[SERIAL_NUMBER] = item[SERIAL_NUMBER]+'_trans_splicing'
        output.write(item[CHR] +'\t' +str(item[LOCA])+'\t'+str(item[LOCB])+'\t'+item[SERIAL_NUMBER] +'\t'+'1000\t'+item[STRAND]+  '\t'+ str(item[LOCA])+'\t'+str(item[LOCB])+'\t0,0,255\t1\t'+str(int(item[LOCB])-int(item[LOCA]))+'\t0\t0\t0\n')

    output.close()


def annotate_circles_sense_antisense():
    os.system('awk \'{print $1"\t"$2"\t"$3"\t*\t*\t"$6}\' '+args.circ+' > tmp/_circs_for_intersect.txt')
    #get the intersection of all circles with annotation
    os.system('bedtools intersect -a  tmp/_circs_for_intersect.txt -b '+args.annotation +
              ' -wa -wb -f 1.0 > tmp/_intersect_circles_with_annotation.bed')
    #get the intersection of circles with the same strand as in the annotation
    os.system('bedtools intersect -a  tmp/_circs_for_intersect.txt -b '+args.annotation +
              ' -wa -wb -f 1.0 -s > tmp/_intersect_circles_with_annotation_stranded.bed')

    #find all circles in the annotated circles file which do not appear in than the stranded circles file
    os.system('grep -Fxv -f tmp/_intersect_circles_with_annotation_stranded.bed tmp/_intersect_circles_with_annotation.bed > tmp/_circles_not_stranded.bed')

    circles_not_stranded_file = open('tmp/_circles_not_stranded.bed','r')
    circles_not_stranded = circles_not_stranded_file.read().strip('\n').split('\n')
    stranded_circles_file = open('tmp/_intersect_circles_with_annotation_stranded.bed','r')
    stranded_circles = stranded_circles_file.read().strip('\n').split('\n')
    circles_not_stranded_file.close()
    stranded_circles_file.close()

    antisense_circs = []
    to_remove = []
    
    #print('####################################')
    if len(circles_not_stranded)>1:
      for i in circles_not_stranded: #each i is the circle+annotation(from BED12)
          circle = i.split('\t')[:6]
          #print circle
          if not i in antisense_circs: #we didn't add this circ to the list yet
              if i.split('\t')[5] != i.split('\t')[11]: #the strand in circle and in annotation are different
                  if not i in stranded_circles: #the circle does not have another annotation in the stranded intersection
                    antisense_circs.append(circle)
          else:
              to_remove.append(circle)
    #print('####################################')
    for i in to_remove:
        if i in antisense_circs:
           antisense_circs.remove(i)
    print "writing to file"
    antisense_circs_file =  open('tmp/antisense_circles.bed','w')
    for as_c in antisense_circs:
        #f.write("\n".join(map(lambda x: str(x), mylist)))
        antisense_circs_file.write("\t".join(map(lambda  x:str(x),as_c)))
        antisense_circs_file.write("\n")
    antisense_circs_file.close()
    os.system('uniq tmp/antisense_circles.bed > tmp/antisense_circles_uniq.bed')


    os.system('awk \'NR==FNR{c[$1$2$3$6]++;next};c[$1$2$3$6] > 0 \' tmp/antisense_circles_uniq.bed '+args.circ+
              ' |awk \'{print $1"\t"$2"\t"$3"\t"$4"_antisense\t"$5"\t"$6}\' > tmp/antisense_circles_uniq_with_original_names.bed')
    os.system('awk \'NR==FNR{c[$1$2$3$6]++;next};c[$1$2$3$6] > 0 \' tmp/antisense_circles_uniq.bed '+args.circ+
              ' |awk \'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}\' > tmp/antisense_circles_uniq_with_names.bed')
    os.system(' grep -Fxv -f tmp/antisense_circles_uniq_with_names.bed '+args.circ+' > tmp/sense_circles.bed')
    os.system('cat tmp/sense_circles.bed tmp/antisense_circles_uniq_with_original_names.bed > tmp/final_annotated_circles.bed')

    args.circ = 'tmp/final_annotated_circles.bed'
    os.remove('tmp/antisense_circles.bed')
    os.remove('tmp/antisense_circles_uniq.bed')
    os.remove('tmp/antisense_circles_uniq_with_names.bed')
    os.remove('tmp/antisense_circles_uniq_with_original_names.bed')
    os.remove('tmp/_circles_not_stranded.bed')
    os.remove('tmp/_circs_for_intersect.txt')
    os.remove('tmp/_intersect_circles_with_annotation.bed')
    os.remove('tmp/_intersect_circles_with_annotation_stranded.bed')
    os.remove('tmp/sense_circles.bed')


def find_exons():
    output = open('tmp/exons.bed', 'w')
    circsH = open(args.circ, 'r')
    out_circ = open('tmp/_CIRC_TEMP_.bed', 'w')
    # filtered_out = open('tmp/filtered_out_circles.bed', 'w')
    # filtered_out_ex = open('tmp/filtered_out_circles_EXACT.bed', 'w')
    # counter for filtered out circles
    numFiltered = 0
    # counter for filtered out circles in the EXACT case.
    numFiltered_ex = 0

    # Create temporary modified circles .bed file:
    circles = list()
    # We want each row in the new file have the <chr>\t<locA>\t<locB>\t<ID>
    # so we create a counter.
    counter = 0
    for circle in circsH:
        circle = circle.strip('\n').strip('\t').split('\t')
        # Filter out circles with zero reads.
        circles.append(circle[START:BED6_NUM_COL])
        # write the rows!
        out_circ.write(circle[CHR]+'\t'+circle[LOCA]+'\t'+circle[LOCB]+'\t'+str(counter)+'\n')
        counter += 1
    out_circ.close()

    # Intersect circ and annotation.
    # In the intersection file each circle has the annotation lines containing it (-f 1.0).

    os.system('bedtools intersect -a tmp/_CIRC_TEMP_.bed -b '+args.annotation +
              ' -wa -wb -f 1.0 > tmp/_INTERSECTION_.bed')
    # The following row of code is garbage deletion but leaving it is good for debugging.
    ###########os.remove('tmp/_CIRC_TEMP_.bed')

    # Open the intersection file we created.
    inerxH = open('tmp/_INTERSECTION_.bed', 'r')
    inerx = inerxH.read().strip('\n').split('\n')
    inerxH.close()
    ###########os.remove('tmp/_INTERSECTION_.bed')
    # No circles scanned yet (@ = 'none'):
    current_circle = '@'
    winner_exons = list()
    winner_introns = list()
    other_trans_names = list()
    winner_score = -1
    winner_name = ''
    #before looking for the best transcript,
    # we have to remove all the annotated anti-sense circles
    for row in inerx:
        row = row.strip('\n').strip('\t').split('\t')
        # We scan the rows, each time when the circle changes (ID)
        # We make the new circle the current circle and update
        # the winning transcript as the transcript of the prev circle exons.
        if current_circle != row[SERIAL_NUMBER]:
            # Just a small check that this is not the first circle
            # when we dont yet have any previous circles.
            if current_circle != '@':
                # add data to circles:
                # calculate blockStarts (offsets, one of the columns in the BED specification)
                blockStarts = list()
                blockStarts.append(0)
                for i in range(1, len(winner_exons)):
                    blockStarts.append(blockStarts[i-1]+int(winner_exons[i-1])+int(winner_introns[i-1]))
                # Save the number of exons in the winner transcript.
                ex_len = len(winner_exons)
                # Convert all data to comma separated strings.
                winner_exons = ",".join(winner_exons)
                blockStarts = ",".join(map(str, blockStarts))
                # Add all data to the circles(circle per row) array
                # It is - [0,           1,                 2,           3,...] -> circles list of lists.
                #          |             |                  |           |,...
                #          v             v                  v           v,...
                #  [circ. columns]  [circ. columns]  [circ. columns] ...
                circles[int(current_circle)].append(str(ex_len))
                circles[int(current_circle)].append(winner_exons)
                circles[int(current_circle)].append(blockStarts)
                circles[int(current_circle)].append(str(winner_score))
                circles[int(current_circle)].append(winner_name)
                circles[int(current_circle)].append(','.join(other_trans_names))
                #Set up next circle:
                winner_exons = list()
                winner_introns = list()
                winner_score = -1
                winner_name = ''
            current_circle = row[SERIAL_NUMBER]
            locA = row[LOCA]
            locB = row[LOCB]

        #scan exons and introns
        exons = list()
        introns = list()
        score = 0
        name = row[TRANSCRIPT_NAME]
        starts = row[EXONS_STARTS].strip(',').split(',')
        ends = row[EXONS_ENDS].strip(',').split(',')
        # Find the exon in the current transcript that starts
        # as close as possible from the "left"(lower) to the circle
        # start location.
        placeStart = bisect.bisect_left(map(int, starts), int(locA))
        # StartLoop -  keeps track if the exons scanning loop
        # should start at all, as sometimes we finish the whole
        # circle exons scan in the special case below.
        startLoop = True
        # The case where the circle starts after all the exons in this transcript
        # may be either IN the last exon (Internal) or not at all.
        if placeStart == len(starts):
            placeStart -= 1
        # The awesome bull's eye case (starts exactly at exon)
        # In this case the score is at least 1.
        if starts[placeStart] == locA:
            score += 1
            low = placeStart
        else:  # does not start at exon
            low = placeStart-1
            if int(ends[low]) < int(locA):  # we fell in intron-before
                # print "fell in intron before..."
                if placeStart+1 < len(starts):
                    # Next exon is out of our scope (finished here)
                    if int(starts[placeStart+1]) >= int(locB):
                        # maybe it ENDS at end of exon though...
                        if int(ends[placeStart]) == int(locB):
                            score += 1
                        exons.append(str(int(locB) - int(locA)))
                        introns.append('0')
                        # no need to loop over other exons, already done.
                        startLoop = False
                if startLoop:
                    # the case where the end is at the intron after the first exon.
                    if int(ends[placeStart]) >= int(locB):
                        # or might be just a bull's eye at the end.
                        if int(ends[placeStart]) == int(locB):
                            score += 1
                        exons.append(str(int(locB) - int(locA)))
                        introns.append('0')
                        startLoop = False
                    else:
                        # intron and exon as first circ exon.
                        exons.append(str(int(ends[placeStart]) - int(locA)))
                        introns.append(str(int(starts[placeStart+1]) - int(ends[placeStart])))
                        low = placeStart + 1
            else:  # circle start in middle of exon
                if low+1 < len(starts):
                    if int(starts[low+1]) >= int(locB):
                        if int(ends[low]) == int(locB):
                            score += 1
                        exons.append(str(int(locB) - int(locA)))
                        introns.append('0')
                        startLoop = False
                if startLoop:
                    if int(ends[low]) >= int(locB):
                        if int(ends[low]) == int(locB):
                            score += 1
                        #and the end is in the same exon!
                        exons.append(str(int(locB) - int(locA)))
                        introns.append('0')
                        startLoop = False
                    else:
                        # Sanity check:
                        if (int(ends[low]) - int(locA)) < 0:
                            print "Negative!!! (start mid ex) at " + locA +":"+locB
                        exons.append(str(int(ends[low]) - int(locA)))
                        introns.append(str(int(starts[placeStart]) - int(ends[low])))
                        low = placeStart
        # Now the actual loop. After all the start special cases,
        # we scan exon by exon until we finish the circle.
        # However, we also have to check here the END special cases.
        while startLoop:
            # No more exons left
            if low >= len(starts):
                break
            # Just some array out-of-bound problems fix.
            if low + 1 < len(starts):
                nextStart = starts[low+1]
            else:
                nextStart = ends[low]
            # ended in middle of exon
            if int(locB) < int(ends[low]):
                # Sanity check:
                if (int(locB) - int(starts[low])) < 0:
                    print "Negative!!! (end mid ex) at " + locA +":"+locB
                exons.append(str(int(locB) - int(starts[low])))
                introns.append('0')
                break
            # last exon
            if int(locB) < int(nextStart) or int(locB) == int(ends[low]):
                # Maybe bull's eye at end?
                if int(locB) == int(ends[low]):
                    score += 1
                # Another sanity check:
                if (int(locB) - int(starts[low])) < 0:
                    print "Negative!!! (last ex) at " + locA + ":" + locB
                # exon and intron as last exon
                exons.append(str(int(locB) - int(starts[low])))
                introns.append('0')
                break
            #Just another exon:
            # sanity check:
            if (int(ends[low]) - int(starts[low])) < 0:
                    print "Negative!!! (reg ex) at " + locA + ":" + locB
            exons.append(str(int(ends[low]) - int(starts[low])))
            introns.append(str(int(nextStart) - int(ends[low])))
            low += 1
        # Keep track of the winner transcript, if this one is better
        # make it the winner.
        #print name, score, locA, locB
        if score > winner_score:
            other_trans_names = [name]
            winner_score = score
            winner_name = name
            winner_exons = copy.deepcopy(exons)
            winner_introns = copy.deepcopy(introns)
        elif score == winner_score:  # Maybe >=  ???
            other_trans_names.append(name)


    # Now when all the calculations are done, we can write output
    # and find statistics:
    go_to_next_circ = False
    numOK = 0
    filtered_circs=[]
    for circle in circles:
        go_to_next_circ = False
        # If a circle has less than 7 columns, it means it should be filtered
        # As no transcripts where found for it:
        if len(circle) > 7:
            for i in range(len((circle[CIRCLE_EXONS_STARTS]).split(','))):
                if circle[CIRCLE_EXONS_STARTS].split(',')[i] == "0":
                    #do not print the circle
                    go_to_next_circ = True
                    break
            # Make the EXACT filtration.
            if go_to_next_circ == True:
                continue
            if circle[SCORE] != str(EXACT_SCORE):
                # filtered_out_ex.write(circle[CHR]+'\t'+circle[LOCA]+'\t'+circle[LOCB]+'\t'+circle[NUM_READS]+'\t'+
                #                       circle[GENE_NAME]+'\n')
                numFiltered_ex += 1
            # reorganize all circle data into BED format columns:
            # (4 and 10 are merged - those are the CG name and gene-name)
            circle_line = \
                circle[CHR]+'\t'+circle[LOCA]+'\t'+circle[LOCB]+'\t'+circle[GENE_NAME]+':' + \
                circle[CHOSEN_TRANSCRIPT]+'\t'+circle[NUM_READS]+'\t'+\
                circle[STRAND]+'\t'+circle[LOCA]+'\t'+circle[LOCB]+'\t0,0,255\t'+circle[NUM_EXONS]+'\t' + \
                circle[CIRCLE_EXONS_STARTS]+'\t' + \
                circle[CIRCLE_EXONS_ENDS]+'\t'+circle[SCORE]+'\t'+circle[ALTERNATIVE_TRANSCRIPTS]
            output.write(circle_line + '\n')
            numOK += 1
        else: # filter it
            filtered_circs.append(circle)
            numFiltered += 1
            numFiltered_ex += 1
    # filtered_out.close()
    # filtered_out_ex.close()
    output.close()
    # Write all possible outputs one may need:
    print "Creating EXACT (only score 2) output..."
    os.system("cat tmp/exons.bed | awk '$13 == 2 {print $0}' > tmp/exons_EXACT.bed")
    sort_filtered_circs(filtered_circs)

    # print "Creating BED12 file..."
    # os.system("cat tmp/exons.bed | cut -f 1-12 > tmp/exons_bed12.bed")
    # print "Creating BED6 file..."
    # os.system("cat tmp/exons.bed | cut -f 1-12 | bed12ToBed6 > tmp/exons_bed6.bed")
    # print "Creating modified circles file..."
    # os.system("grep -F -v -f tmp/filtered_out_circles.bed "+args.circ+" > tmp/filtered_"+args.circ)
    # print "Creating modified EXACT circles file..."
    # os.system("grep -F -v -f tmp/filtered_out_circles_EXACT.bed "+args.circ+" > tmp/filtered_EXACT_"+args.circ)
    print "Finished."
    print str(numOK) + " circles where annotated and splitted to exons."
    print "Of them, "+str(numFiltered_ex-numFiltered) + " where INTERNAL (non-EXACT)."
    print "Output can be found at -> exons.bed"
    print str(numFiltered) + " circles where filtered in the process"  # , they can be found at -> filtered_out_circles.bed"
    print str(numFiltered_ex) + " circles where filtered in the process for EXACT"  # , they can be found at " \
    #                           "-> filtered_out_circles_EXACT.bed"


def remove_temp():
    if not args.keep_temporary:
        os.system('rm -rf tmp')


def reverse_complement(seq):
    bases = list(seq.upper())
    bases = reversed([COMPLEMENTS[base] for base in bases])
    return ''.join(bases)


def get_sequence(chrm, start, end):
    return fasta[index.index(chrm)][start:end]


# Checking the quality of the read in .sam file
# @param:
#       row_split - the sam line split by '\t' delim.
# @return:
#       quality in given quality system as int
def check_quality(row_split):
    return reduce(lambda x, y: x+(ord(y)-int(args.quality)),
                  row_split[QUAL_IN_SAM], 0)/len(row_split[QUAL_IN_SAM]) > int(args.min_quality)


# Comparing two sequences, one form the res_filtered file and one from the index file
# @param:
#       res_seq - the entire sequence from the res_filtered row
#       res_start_pt - where the alignment starts int the res_filtered row
#       index_seq - the seed surrounding the junction
# @return:
#       -1 of no match, 0 if matches as SENSE and 16 if matches as ANTI-SENSE.
def compare_seq(res_seq, res_start_pt, index_seq, half):
    junction_area = index_seq[half-args.junction:half+args.junction].upper()
    res_start_pt = half-args.junction-res_start_pt
    if res_start_pt+2*args.junction+1 > len(res_seq) or res_start_pt+1 < 0:
        return NO_MATCH
    res_seq = res_seq[res_start_pt+1:res_start_pt+2*args.junction+1].upper()
    if junction_area == res_seq:
        return MATCH_SENSE
    if reverse_complement(junction_area) == res_seq:
        return MATCH_ANTI_SENSE
    return NO_MATCH

#def compare_seq(res_seq, res_start_pt, index_seq):
#    res_end_pt = res_start_pt + 2 * args.junction
#    if res_end_pt > len(index_seq):
#        return NO_MATCH
#    # Cutting the res_seq to be of the same length of index_seq from the alignment starting point.
#    aligned_seq = index_seq[res_start_pt-1:res_end_pt-1].upper()
#    res_seq = res_seq[:len(aligned_seq)].upper()
#    if aligned_seq == res_seq:
#        return MATCH_SENSE
#    if reverse_complement(aligned_seq) == res_seq:
#        return MATCH_ANTI_SENSE
#    return NO_MATCH

# Adding to the left and right linear index files the junctions from pcp pipeline.
# @param:
#       output_left - the linear index left output file.
#       output_right - the linear index right output file
#       half - the length of an index seq.
def index_based_on_pcp(output_left, output_right, half):
    if args.pcp is None:
        return
    print
    print "PCP data present. Adding to linear index..."
    dict_lin_starts = dict()
    dict_lin_ends = dict()
    circles = open(args.exons, 'r').read().strip('\n').split('\n')
    pcp_annot = open(args.pcp, 'r').read().strip('\n').split('\n')
    print
    print "Reading PCP linear annotation..."
    for line in pcp_annot:
        line = line.strip('\t').split('\t')
        chrm = line[CHR]
        tx_start = line[LOCA]
        tx_end = line[LOCB]
        strand = line[STRAND]
        stamp_start = (chrm, tx_start, strand)
        stamp_end = (chrm, tx_end, strand)
        # BY START
        if stamp_start not in dict_lin_starts:
            dict_lin_starts[stamp_start] = [line]
        else:
            dict_lin_starts[stamp_start].append(line)
        # BY END
        if stamp_end not in dict_lin_ends:
            dict_lin_ends[stamp_end] = [line]
        else:
            dict_lin_ends[stamp_end].append(line)
    # Scan circles:
    print "Scanning circles and adding to index..."
    left_added = right_added = 0
    for circ in circles:
        circ = circ.strip('\t').split('\t')
        chrm = circ[CHR]
        tx_start = circ[LOCA]
        tx_end = circ[LOCB]
        strand = circ[STRAND]
        #AVIGAYEL
        circle_header = circ[GENE_NAME]+';'+chrm+":"+tx_start+'-'+tx_end+";"+strand
        stamp_left = (chrm, tx_start, strand)
        stamp_right = (chrm, tx_end, strand)
        if stamp_right in dict_lin_starts:
            for lin in dict_lin_starts[stamp_right]:
                output_right.write('>'+circle_header+'\n')
                output_right.write(get_sequence(chrm, int(lin[LOCA])-half, int(lin[LOCA])) +
                                   get_sequence(chrm, int(lin[LOCB]), int(lin[LOCB])+half)+'\n')
                right_added += 1
                
        if stamp_left in dict_lin_ends:
            for lin in dict_lin_ends[stamp_left]:
                output_left.write('>'+circle_header+'\n')
                output_left.write(get_sequence(chrm, int(lin[LOCA])-half, int(lin[LOCA])) +
                                  get_sequence(chrm, int(lin[LOCB]), int(lin[LOCB])+half)+'\n')
                left_added += 1
                print get_sequence(chrm, int(lin[LOCA])-half, int(lin[LOCA])) + get_sequence(chrm, int(lin[LOCB]), int(lin[LOCB])+half)+'\n'
                
    #AVIGAYEL close files
    output_left.close()
    output_right.close()
    print "Added to left:", left_added
    print "Added to right:", right_added
    print "Finished populating linear index (left and right)."


def run_circular(samp):
    # Read all input data and open output file
    output = open(args.output.replace('.sam', '_circ.sam'), 'w')
    index = dict()
    print "Reading input data..."
    # remove last two columns so -split will work
    os.system('cut -f1-12 '+args.exons+' > '+args.exons+'12')

    os.system('bedtools getfasta -fi '+args.system+'.fa -bed '+args.exons+'12 -s -split -name -fo tmp/getFastaOutput.fa')
    #os.system('bedtools getfasta -fi '+args.system+'.fa -bed '+args.exons+'6 -s -split -name -fo tmp/getFastaOutput.fa')
    if args.data.endswith(".gz"):
        data = gzip.open(args.data, 'r')
    else:
        data = open(args.data, 'r')
    data.readline()
    read_size = len(data.readline())-1
    print "Read size is: ", read_size
    if args.read_len is None:
        half = read_size - args.junction
        data.close()
        print "Index size was determined to be:", half*2
    else:
        # No empty reads
        if args.read_len < 2:
            print "WARNING: Read length cannot be < 2 !!! Auto recalculation..."
            half = read_size - args.junction
            data.close()
            print "Index size was determined to be:", half*2
        else:
            half = args.read_len-args.junction
            print "Index size was FORCED-SET to be:", half*2
    # Create both file and in-memory dict index:
    output_index = open('tmp/__index__.fa', 'w')
    exons = open(args.exons, 'r').read().strip('\n').split('\n')
    seqs = open('tmp/getFastaOutput.fa', 'r').read().strip('\n').split('\n')
    print
    print "Building index in fasta format to: tmp/__index__.fa"
    exons = [exon.strip('\t').split('\t') for exon in exons]
    seqs = [seqs[i] for i in xrange(1, len(seqs), 2)]
    junction = [seqs[i][-half:]+seqs[i][:half] for i in range(len(seqs))]
    for i in range(len(exons)):
        #avigayel add ';'+exons[i][STRAND]
        index[exons[i][GENE_NAME]+';'+exons[i][CHR]+':'+exons[i][LOCA]+'-'+exons[i][LOCB]+';'+exons[i][STRAND]] = junction[i]
        #avigayel add  ';'+exons[i][STRAND]
        output_index.write('>'+exons[i][GENE_NAME]+';'+exons[i][CHR]+':'+exons[i][LOCA]+'-'+exons[i][LOCB]+';'+exons[i][STRAND]+'\n')
        output_index.write(junction[i]+'\n')
    output_index.close()

    # All the alignments:
    print
    print "Building BT2 index..."
    os.system("bowtie2-build -fq tmp/__index__.fa tmp/junctions_"+str(half))
    print
    print "Aligning to index..."
    params = "" if args.N == '0' else ("-N "+args.N)

    # align to index:--very-sensitive removed
    print '''bowtie2 '''+params+''' -x tmp/junctions_''' + str(half) + ''' -U '''+args.data+''' -p '''+str(args.numThread)+''' | samtools view -hbuS - | samtools view -F4 - > tmp/__res__'''+samp+'''__.sam\n'''
                            
    os.system('''bowtie2 '''+params+''' -x tmp/junctions_''' + str(half) +
              ''' -U '''+args.data+''' -p '''+str(args.numThread)+''' |
              samtools view -hbuS - | samtools view -F4 - > tmp/__res__'''+samp+'''__.sam''')

    # convert to fastq:
    os.system('''cat tmp/__res__'''+samp+'''__.sam | grep -v ^@ |
              awk '{print "@"$1"\\n"$10"\\n+\\n"$11}' > tmp/__res__.fq''')
    print
    print "Aligning to genome..."
    print "bowtie2 -x "+args.system +  " " + params + " -U tmp/__res__.fq |" + "samtools view -S - > tmp/__aligned_genome__"+samp+".sam\n"
    # Align to genome
    os.system("bowtie2 -x "+args.system +
              " " + params + " -U tmp/__res__.fq |"
              "samtools view -S - > tmp/__aligned_genome__"+samp+".sam")

    print
    print "Aligning to transcriptome..."
    print "bowtie2 -x " + args.transcriptome + " " + params + " -U tmp/__res__.fq |" + "samtools view -S - > tmp/__aligned_transcriptome__"+samp+".sam\n"          
    # Align to transcriptome:
    os.system("bowtie2 -x " + args.transcriptome +
              " " + params + " -U tmp/__res__.fq |"
              "samtools view -S - > tmp/__aligned_transcriptome__"+samp+".sam")

    triplets = zip(
        open('tmp/__res__'+samp+'__.sam', 'r').read().strip('\n').split('\n'),
        open('tmp/__aligned_genome__'+samp+'.sam', 'r').read().strip('\n').split('\n'),
        open('tmp/__aligned_transcriptome__'+samp+'.sam', 'r').read().strip('\n').split('\n')
    )
    allowed_edit_distance = int(args.mis)+2 #avigayel: to allow 2 indels
    res_filtered = open('tmp/__res_filtered__'+samp+'__.sam', 'w')  #the reads that we want to save!
    safety_file = open(args.output.replace('.sam', '_safetyDistances.txt'), 'w')
    for triplet in triplets:
        edit_distance_index = re.search(EDIT_DISTANCE_REGEX, triplet[INDEX])
        if edit_distance_index is None:
            continue
        else:
            edit_distance_index = int(edit_distance_index.group(1))


        edit_distance_genome = re.search(EDIT_DISTANCE_REGEX, triplet[GENOME])
        if edit_distance_genome is None:
            edit_distance_genome = sys.maxint
        else:
            edit_distance_genome = int(edit_distance_genome.group(1))

        edit_distance_trans = re.search(EDIT_DISTANCE_REGEX, triplet[TRANSCRIPTOME])
        if edit_distance_trans is None:
            edit_distance_trans = sys.maxint
        else:
            edit_distance_trans = int(edit_distance_trans.group(1))

        if edit_distance_index <= allowed_edit_distance:   #the read aligned to the index.
            safety_distance = min(edit_distance_genome - edit_distance_index,
                                  edit_distance_trans - edit_distance_index)
            if safety_distance >= args.safety:
                if edit_distance_trans == edit_distance_genome == sys.maxint: # did not align to genome or transcriptome
                    safety_distance = 'INF'                                                  #avigayel: if the read aligns only to the index then continue to next read
                res_filtered.write(triplet[INDEX]+'\n')
                head = triplet[INDEX].strip('\t').split('\t')
                safety_file.write(head[NAME_IN_SAM]+'\t' +
                                  head[SIGNATURE_IN_SAM]+'\t' +
                                  str(safety_distance)+'\n')
    res_filtered.close()
    safety_file.close()

    # Filter by junction seed containment:
    res = open('tmp/__res_filtered__'+samp+'__.sam', 'r').read().strip('\n').split('\n')
    log = open('tmp/__filteredReadsLog__'+samp+'.txt', 'w')
    num_not = 0
    out_place = open(args.output.replace('.sam', '_junctionLocation.txt'), 'w')
    output_counts = open(args.output.replace('.sam', '_counts.txt'), 'w')
    output_uniq_counts = open(args.output.replace('.sam', '_uniq_counts.txt'), 'w')
    count_dict = dict()
    uniq_dict = dict()
    for key in index:
        count_dict[key] = 0
        uniq_dict[key] = set()
    for row in res:
        row_split = row.split('\t')
        if len(row_split) < 2:
            #print "Filtered out: lenght"
            #print row_split
            log.write("Filtered out: lenght" + '\n')
            log.write('\t'.join(row_split))
            
            continue
        name = row_split[NAME_IN_SAM]
        seq = row_split[SEQ_IN_SAM]
        start_pt = int(row_split[ALIGN_START_IN_SAM])
        if not check_quality(row_split):
            #print "Filtered out: quality"
            #print row_split
            log.write("Filtered out: quality" + '\n')
            log.write('\t'.join(row_split) + "\n")
           
            continue
        if half-args.junction < 0 or half+args.junction > len(index[name]):
            #print "Filtered out: Junction fell outside"
            #print row_split
            log.write("Filtered out: Junction fell outside" + '\n')
            log.write('\t'.join(row_split) + "\n")
            
            continue
        junc = index[name][half-args.junction:half+args.junction]  # fixed junction size.
        # Check if .SAM cigar has any flag except for alignment mismatch ('M').
        only_mis = not any([cig in row_split[CIGAR_IN_SAM].upper() for cig in CIGAR_FLAGS])
        # Checking seed containment in res_filtered's sequences
        strand_type = compare_seq(seq, start_pt, index[name], half)
        if not only_mis:
            #print "Filtered out: CIGAR includes deletions and insertions"
            #print row_split
            log.write("Filtered out: CIGAR includes more than two indels" + '\n')
            log.write('\t'.join(row_split) + "\n")

        if strand_type == NO_MATCH:
            #print "Filtered out: NO SEED MATCH"
            #print row_split
            log.write("Filtered out: strand_type == NO_MATCH" + '\n')
            log.write('\t'.join(row_split) + "\n")


        #if only_mis and strand_type != NO_MATCH:  # check if matches, and if so if it's a SENSE or ANTI-SENSE strand
	if  strand_type != NO_MATCH:  # check if matches, and if so if it's a SENSE or ANTI-SENSE strand
            row_split[STRAND_IN_SAM] = str(strand_type)
            row = "\t".join(row_split)
            count_dict[name] += 1
            uniq_dict[name].add(seq.upper())
            output.write(row + '\n')
            outputAll.write(row + '\n')
            out_place.write(name+'\t'+str(float(seq.find(junc)+len(junc)/2)/float(len(seq)))+'\n')
        else:
            num_not += 1
    print
    print "Not containing junction: ", num_not
    for key in index:
        output_counts.write(key+'\t'+str(count_dict[key])+'\n')
        output_uniq_counts.write(key+'\t'+str(len(uniq_dict[key]))+'\n')
    output.close()
    out_place.close()
    output_counts.close()
    output_uniq_counts.close()


def run_linear(samp):
    output_left = open('tmp/__output_left__.fa', 'w')
    output_right = open('tmp/__output_right__.fa', 'w')

    # Read all input data and open output file
    index = dict()
    print "Reading input data..."
    if args.data.endswith(".gz"):
        data = gzip.open(args.data, 'r')
    else:
        data = open(args.data, 'r')
    data.readline()
    read_size = len(data.readline()) - 1
    print "Read size is: ", read_size
    if args.read_len is None:
        half = read_size - args.junction
        data.close()
        print "Index size was determined to be:", half*2
    else:
        # No empty reads.
        if args.read_len < 2:
            print "WARNING: Read length cannot be < 2 !!! Auto recalculation..."
            half = read_size - args.junction
            data.close()
            print "Index size was determined to be:", half*2
        else:
            half = args.read_len-args.junction
            print "Index size was FORCED-SET to be:", half*2
    # create intersection:
    intx = subprocess.check_output('bedtools intersect -wao -s -f 1.0 -a '+args.exons+' -b '+args.annotation +
                                   ' | cut -f 1-6,13-25', shell=True)
    intx = intx.strip('\n').split('\n')
    intx = [row.strip('\t').split('\t') for row in intx]
    for row in intx:
        strand = row[STRAND]
        chrm = row[CHR]
        circle = row[GENE_NAME]+';'+chrm+":"+row[LOCA]+'-'+row[LOCB]+';'+strand
        locA = int(row[LOCA])
        locB = int(row[LOCB])
        if '.' in row[RUN_LINEAR_EXON_STARTS]:
            continue
        exons_starts = map(int, row[RUN_LINEAR_EXON_STARTS].strip(',').split(','))
        exons_ends = map(int, row[RUN_LINEAR_EXON_ENDS].strip(',').split(','))
        # if args.only_annotated and int(row[RUN_LINEAR_SCORE]) < EXACT_SCORE:
        #     continue
        plausible_transcripts = row[RUN_LINEAR_ALTERNATIVE_TRANSCRIPTS].strip(',').split(',')
        if row[RUN_LINEAR_TRANSCRIPT] not in plausible_transcripts:
            continue
        # get exonic sequence from the end and from the start
        exonic_tail = ''
        take_serial = False
        for ind in range(len(exons_starts)):
            if not take_serial:
                if exons_starts[ind] <= locB < exons_ends[ind]:
                    exonic_tail += get_sequence(chrm, locB, exons_ends[ind])
                    take_serial = True
                if ind+1 < len(exons_starts):
                    if exons_ends[ind] <= locB < exons_starts[ind+1]:
                        take_serial = True
            else:
                exonic_tail += get_sequence(chrm, exons_starts[ind], exons_ends[ind])

        exonic_head = ''
        for ind in range(len(exons_starts)):
            if locA < exons_ends[ind]:
                if locA < exons_starts[ind]:
                    break
                else:
                    exonic_head += get_sequence(chrm, exons_starts[ind], locA)
            else:
                exonic_head += get_sequence(chrm, exons_starts[ind], exons_ends[ind])
        # Now we have exonic seq before and after the circle.
        #avigayel-add strand to circ name
        if strand == '+':
            if len(exonic_head) >= half:
                output_left.write('>'+circle+'\n')
                output_left.write(exonic_head[-half:] +
                                  get_sequence(chrm, locA, locA+half)+'\n')
            else:
                continue
            if len(exonic_tail) >= half:
                output_right.write('>'+circle+'\n')
                output_right.write(get_sequence(chrm, locB-half, locB) +
                                   exonic_tail[:half]+'\n')
            else:
                continue
        else:
            if len(exonic_tail) >= half:
                output_left.write('>'+circle+'\n')
                output_left.write(get_sequence(chrm, locB-half, locB) +
                                  exonic_tail[:half]+'\n')
            else:
                continue
            if len(exonic_head) >= half:
                output_right.write('>'+circle+'\n')
                output_right.write(exonic_head[-half:] +
                                   get_sequence(chrm, locA, locA+half)+'\n')
            else:
                continue
    index_based_on_pcp(output_left, output_right, half)
    output_right.close()
    output_left.close()
    ## linaer sam files go to tmp/
    output_left = open(args.output.replace('.sam', '_LinearLeft.sam'), 'w')
    output_right = open(args.output.replace('.sam', '_LinearRight.sam'), 'w')
    os.system('sh tmp/remove_duplicates_fa.sh output')
    left_dict = dict()
    right_dict = dict()
    left_file = open('tmp/__output_left_no_dupl_num__.fa', 'r').read().strip('\n').split('\n')
    right_file = open('tmp/__output_right_no_dupl_num__.fa', 'r').read().strip('\n').split('\n')
    for lefty_head, lefty_seq in zip(left_file[::2], left_file[1::2]):
        key = lefty_head[1:]
        if key in left_dict:
            left_dict[key].append(lefty_seq)
        else:
            left_dict[key] = [lefty_seq]
    for righty_head, righty_seq in zip(right_file[::2], right_file[1::2]):
        key = righty_head[1:]
        if key in right_dict:
            right_dict[key].append(righty_seq)
        else:
            right_dict[key] = [righty_seq]
    print
    print "Building index in fa format to: tmp/__output_left_no_dupl_num__.fa, tmp/__output_right_no_dupl_num__.fa"
    # All the alignments:
    print
    print "Building BT2 index..."
    os.system("bowtie2-build -fq tmp/__output_left_no_dupl_num__.fa tmp/junctions_lin_left_"+str(half))
    os.system("bowtie2-build -fq tmp/__output_right_no_dupl_num__.fa tmp/junctions_lin_right_"+str(half))
    print
    print "Aligning to index..."
    params = "" if args.N == '0' else ("-N "+args.N)
    # Left
    os.system('''bowtie2 --very-sensitive ''' + params + ''' -x tmp/junctions_lin_left_''' + str(half) +
              ''' -U '''+args.data+''' -p '''+str(args.numThread)+''' |
              samtools view -hbuS - | samtools view -F4 - > tmp/__res_lin_left__'''+samp+'''__.sam''')
    os.system('''cat tmp/__res_lin_left__'''+samp+'''__.sam | grep -v ^@ |
              awk '{print "@"$1"\\n"$10"\\n+\\n"$11}' > tmp/__res_lin_left__.fq''')
    # Now right
    os.system('''bowtie2 --very-sensitive ''' + params + ''' -x tmp/junctions_lin_right_''' + str(half) +
              ''' -U '''+args.data+''' -p '''+str(args.numThread)+''' |
              samtools view -hbuS - | samtools view -F4 - > tmp/__res_lin_right__'''+samp+'''__.sam''')
    os.system('''cat tmp/__res_lin_right__'''+samp+'''__.sam | grep -v ^@ |
              awk '{print "@"$1"\\n"$10"\\n+\\n"$11}' > tmp/__res_lin_right__.fq''')
    print
    print "Aligning to genome left..."
    # Left
    os.system("bowtie2 -x "+args.system +
              " " + params + " -U tmp/__res_lin_left__.fq |"
              "samtools view -hbuS - | samtools view -F4 -> tmp/__aligned_genome_left__.sam")
    print
    print "Filtering aligned to genome from index aligned..."
    os.system('''awk '{split($17,a,":"); if(a[3]<='''+args.mis+'''){print $0}}' tmp/__aligned_genome_left__.sam '''
              '''| grep -v ^@ | cut -f 1 | sed "s/$/\t/" > tmp/__aligned_genome_left_filtered__.sam''')
    # Now right
    print
    print "Aligning to genome right..."
    os.system("bowtie2 -x "+args.system +
              " " + params + " -U tmp/__res_lin_right__.fq |"
              "samtools view -hbuS - | samtools view -F4 -> tmp/__aligned_genome_right__.sam")
    print
    print "Filtering aligned to genome from index aligned..."
    os.system('''awk '{split($17,a,":"); if(a[3]<='''+args.mis+'''){print $0}}' tmp/__aligned_genome_right__.sam '''
              '''| grep -v ^@ | cut -f 1 | sed "s/$/\t/" > tmp/__aligned_genome_right_filtered__.sam''')
    print
    print "Removing those from tmp/__res_lin_<side>__"+samp+"__.sam ... "
    # Left
    os.system('''cat tmp/__res_lin_left__'''+samp+'''__.sam |'''
              ''' awk '{split($17,a,":"); if(a[3]<='''+args.mis+'''){print $0}}' |'''
              ''' grep -vFf tmp/__aligned_genome_left_filtered__.sam'''
              ''' > tmp/__res_lin_left_filtered__'''+samp+'''__.sam''')
    # Now right
    os.system('''cat tmp/__res_lin_right__'''+samp+'''__.sam | '''
              '''awk '{split($17,a,":"); if(a[3]<='''+args.mis+'''){print $0}}' |'''
              ''' grep -vFf tmp/__aligned_genome_right_filtered__.sam'''
              ''' > tmp/__res_lin_right_filtered__'''+samp+'''__.sam''')
    # Filter by junction seed containment:
    res_left = open('tmp/__res_lin_left_filtered__'+samp+'__.sam', 'r').read().strip('\n').split('\n')
    res_right = open('tmp/__res_lin_right_filtered__'+samp+'__.sam', 'r').read().strip('\n').split('\n')
    num_not = 0
    output_counts_left = open(args.output.replace('.sam', '_left_counts.txt'), 'w')
    output_counts_right = open(args.output.replace('.sam', '_right_counts.txt'), 'w')
    uniq_counts_left = open(args.output.replace('.sam', '_left_uniq_counts.txt'), 'w')
    uniq_counts_right = open(args.output.replace('.sam', '_right_uniq_counts.txt'), 'w')
    output_counts = open(args.output.replace('.sam', '_allCounts.txt'), 'w')
    count_dict_left = dict()
    count_dict_right = dict()
    uniq_dict_left = dict()
    uniq_dict_right = dict()
    for key in left_dict:
        #avigayel -  change 2 to 3
        count_dict_left[';'.join(key.split(';')[:3])] = 0
        uniq_dict_left[';'.join(key.split(';')[:3])] = set()
    for key in right_dict:
        count_dict_right[';'.join(key.split(';')[:3])] = 0
        uniq_dict_right[';'.join(key.split(';')[:3])] = set()

    # for left
    for row in res_left:
        row_split = row.split('\t')
        if len(row_split) < 2:
            continue
        name = row_split[NAME_IN_SAM]
        seq = row_split[SEQ_IN_SAM]
        start_pt = int(row_split[ALIGN_START_IN_SAM])
        if not check_quality(row_split):
            continue
        junc = [left_dict[name][ind] for ind in range(len(left_dict[name]))]
        # check if .SAM cigar has any flag except for alignment mismatch ('M').
        only_mis = not any([cig in row_split[CIGAR_IN_SAM].upper() for cig in CIGAR_FLAGS])
        #if only_mis:
	if only_mis:
            strand_type_list = map(lambda j: compare_seq(seq, start_pt, j, half), junc)
            try:
                strand_type_list.index(MATCH_SENSE)
                strand_type = MATCH_SENSE
            except ValueError:
                try:
                    strand_type_list.index(MATCH_ANTI_SENSE)
                    strand_type = MATCH_ANTI_SENSE
                except ValueError:
                    num_not += 1
                    continue
        else:
            num_not += 1
            continue
            #avigayel change 2 to 3
        name = ';'.join(name.split(';')[:3])
        count_dict_left[name] += 1
        uniq_dict_left[name].add(seq.upper())
        row_split[STRAND_IN_SAM] = str(strand_type)
        row = "\t".join(row_split)
        output_left.write(row + '\n')
    print
    print "Not containing junction in left: ", num_not
    num_not = 0
    #Right
    for row in res_right:
        row_split = row.split('\t')
        if len(row_split) < 2:
            continue
        name = row_split[NAME_IN_SAM]
        seq = row_split[SEQ_IN_SAM]
        start_pt = int(row_split[ALIGN_START_IN_SAM])
        if not check_quality(row_split):
            continue
        junc = [right_dict[name][ind] for ind in range(len(right_dict[name]))]
        # check if .SAM cigar has any flag except for alignment mismatch ('M').
        only_mis = not any([cig in row_split[CIGAR_IN_SAM].upper() for cig in CIGAR_FLAGS])
        if only_mis:
            strand_type_list = map(lambda j: compare_seq(seq, start_pt, j, half), junc)
            try:
                strand_type_list.index(MATCH_SENSE)
                strand_type = MATCH_SENSE
            except ValueError:
                try:
                    strand_type_list.index(MATCH_ANTI_SENSE)
                    strand_type = MATCH_ANTI_SENSE
                except ValueError:
                    num_not += 1
                    continue
        else:
            num_not += 1
            continue
            #avigayel: change 2 to 3
        name = ';'.join(name.split(';')[:3])
        count_dict_right[name] += 1
        uniq_dict_right[name].add(seq.upper())
        row_split[STRAND_IN_SAM] = str(strand_type)
        row = "\t".join(row_split)
        output_right.write(row + '\n')
    print "Not containing junction in right: ", num_not
    for key in count_dict_left:
        key = ';'.join(key.split(';')[:3])
        output_counts_left.write(key+'\t'+str(count_dict_left[key])+'\n')
        uniq_counts_left.write(key+'\t'+str(len(uniq_dict_left[key]))+'\n')
    for key in count_dict_right:
        key = ';'.join(key.split(';')[:3])
        output_counts_right.write(key+'\t'+str(count_dict_right[key])+'\n')
        uniq_counts_right.write(key+'\t'+str(len(uniq_dict_right[key]))+'\n')
    for key in count_dict_left:
        key = ';'.join(key.split(';')[:3])
        if key not in count_dict_right:
            output_counts.write(key+'\t'+str(count_dict_left[key])+'\n')
        else:
            output_counts.write(key+'\t'+str(max(count_dict_left[key], count_dict_right[key]))+'\n')
    for key in count_dict_right:
        key = ';'.join(key.split(';')[:3])
        if key not in count_dict_left:
            output_counts.write(key+'\t'+str(count_dict_right[key])+'\n')
    output_left.close()
    output_right.close()
    output_counts_left.close()
    output_counts_right.close()
    uniq_counts_left.close()
    uniq_counts_right.close()
    output_counts.close()


def find_ratios(samp):
    print "Running analysis on circular..."
    run_circular(samp)
    print "-----------------------------"
    print
    print "Running analysis on linear..."
    run_linear(samp)
    output_counts_lin_left = {row.split('\t')[0]: row.split('\t')[1] for row in open(
                                  args.output.replace('.sam', '_left_counts.txt'),
                                  'r'
                              ).read().strip('\n').split('\n')}
    output_uniq_counts_lin_left = {row.split('\t')[0]: row.split('\t')[1] for row in open(
                                       args.output.replace('.sam', '_left_uniq_counts.txt'),
                                       'r'
                                   ).read().strip('\n').split('\n')}
    output_uniq_counts_lin_right = {row.split('\t')[0]: row.split('\t')[1] for row in open(
                                       args.output.replace('.sam', '_right_uniq_counts.txt'),
                                       'r'
                                   ).read().strip('\n').split('\n')}
    output_counts_lin_right = {row.split('\t')[0]: row.split('\t')[1] for row in open(
                                   args.output.replace('.sam', '_right_counts.txt'),
                                   'r'
                               ).read().strip('\n').split('\n')}
    output_counts_lin = {row.split('\t')[0]: row.split('\t')[1] for row in open(
                             args.output.replace('.sam', '_allCounts.txt'),
                             'r'
                         ).read().strip('\n').split('\n')}
    output_counts_cir = {row.split('\t')[0]: row.split('\t')[1] for row in open(
                             args.output.replace('.sam', '_counts.txt'),
                             'r'
                         ).read().strip('\n').split('\n')}
    output_uniq_counts_cir = {row.split('\t')[0]: row.split('\t')[1] for row in open(
                             args.output.replace('.sam', '_uniq_counts.txt'),
                             'r'
                         ).read().strip('\n').split('\n')}

    # Accumulate all reads numbers and ratios for all samples.
    # Writing of the structure to file is done in analyze_samples.
    # Structure of samples_dict: |dictionary where:                                                          |
    #                            |             key = circID, value = [list of tuples where each tuple is:    |
    #                            |                                       a 7 dimensional vector holding      |
    #                            |                                       all data about the i'th sample  ]   |
    safe_division = lambda u, v: u/v if v != 0.0 else 'INF'
    intersection = reduce(lambda x, y: dict(x, **y), (
        {x: safe_division(float(output_counts_cir[x]), float(output_counts_lin[x])) for x in output_counts_cir
         if x in output_counts_lin},
        {x: 'INF' for x in output_counts_cir if x not in output_counts_lin},
        {x: 0.0 for x in output_counts_lin if x not in output_counts_cir}))
    for key in intersection:
        res_tuple = (output_uniq_counts_cir[key] if key in output_uniq_counts_cir else 0.0,
                     output_counts_cir[key] if key in output_counts_cir else 0.0,
                     output_counts_lin_left[key] if key in output_counts_lin_left else 0.0,
                     output_uniq_counts_lin_left[key] if key in output_uniq_counts_lin_left else 0.0,
                     output_counts_lin_right[key] if key in output_counts_lin_right else 0.0,
                     output_uniq_counts_lin_right[key] if key in output_uniq_counts_lin_right else 0.0,
                     intersection[key])
        if key not in samples_dict:
            samples_dict[key] = [res_tuple]
        else:
            samples_dict[key].append(res_tuple)

    # Accumulate all safety distances and extract for each circle the best one.
    # This is for the currently running sample. Writing a file for all samples is done, like for the reads,
    # in analyze_samples().
    # Structure of safety_dict: |dictionary where:                                                                |
    #                           |           key = sample_name, value =                                            |
    #                           |                                  |dictionary where:                           | |
    #                           |                                  |           key=circID, value=maximal_safety | |
    safety_dict[samp] = dict()
    for row in open(
            args.output.replace('.sam', '_safetyDistances.txt'),
            'r'
    ).read().strip('\n').split('\n'):
        if(os.stat(args.output.replace('.sam', '_safetyDistances.txt')).st_size == 0):
            break
        name = row.split('\t')[0]
        val = row.split('\t')[2]
        val = int(val) if val != 'INF' else sys.maxint
        if name in safety_dict[samp]:
            safety_dict[samp][name] = max(safety_dict[samp][name], val)
        else:
            safety_dict[samp][name] = val
    print "Done."


def analyze_samples():
    print "Starting analysis..."
    print
    print "Filtering by read quality. Using Phred+"+str(args.quality)
    print "Checking dependencies..."
    samples = args.data.strip().strip(',').split(',')
    missing_file = False
    for sfile in samples:
        if not os.path.isfile(sfile):
            missing_file = True
            sys.stderr.write("ERROR: data file "+sfile+" does not exists!\n")
    if not os.path.isfile("tmp/remove_duplicates_fa.sh"):
        sys.stderr.write("ERROR: The file remove_duplicates_fa.sh wasn't found.\n"
                         "This file is necessary for the script to run.\n")
        missing_file = True
    if not os.path.isdir(os.path.dirname(args.system)):
        missing_file = True
        sys.stderr.write("ERROR: system directory "+args.system+" does not exists!\n")
    if not os.path.isdir(os.path.dirname(args.transcriptome)):
        missing_file = True
        sys.stderr.write("ERROR: transcriptome directory "+args.transcriptome+" does not exists!\n")
    if not os.path.isfile(args.annotation):
        missing_file = True
        sys.stderr.write("ERROR: annotation file "+args.annotation+" does not exists!\n")
    if not os.path.isfile(args.exons):
        missing_file = True
        sys.stderr.write("ERROR: circles exons file "+args.exons+" does not exists!\n")
    if missing_file:
        sys.stderr.write("ABORTING: not all needed files found, please fix.\n")
        exit(-1)
    # the case where all OK:
    print "All dependencies seem to be present. Continuing..."
    print "-----------------------------"
    output_name_legacy = args.output
    for sample in samples:
        args.data = sample
        sample = os.path.split(sample)[1]
        args.output = 'tmp/'+output_name_legacy.replace('.sam', '_'+os.path.splitext(os.path.split(sample)[1])[0] +
                                                        '.sam')
        print "Analyzing", os.path.splitext(os.path.split(sample)[1])[0]+"..."
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        find_ratios(os.path.splitext(os.path.split(sample)[1])[0])
        print "-----------------------------"
        print
    print "Creating final statistics..."
    output_samples = open(output_name_legacy.replace('.sam', '_allSamplesResults.txt'), 'w')
    output_safeties = open(output_name_legacy.replace('.sam', '_allSafetyDistances.txt'), 'w')
    # OR Nice for R
    #AVIGAYEL-add strand to fourth column
    # output_samples.write("chr\tlocA\tlocB\tname\tstrand\t"+"_uniq_circ_reads\t".join(samples)+"_uniq_circ_reads\t" +
    #                      "_circ\t".join(samples)+"_circ\t"+"_linear_left\t".join(samples) + "_linear_left\t"+
    #                      "_linear_uniq_left_reads\t".join(samples)+"_linear_uniq_left_reads\t" +
    #                      "_linear_right\t".join(samples)+"_linear_right\t" +
    #                      "_linear_uniq_right_reads\t".join(samples)+"_linear_uniq_right_reads\t" +
    #                      "_ratio\t".join(samples) + "_ratio\n")
    ## output without uniq linear
    output_samples.write("chr\tlocA\tlocB\tname\tcirc_allSamples\tstrand\tLinear_allSamples\t"+"_circ\t".join(samples)+"_circ\t" +
                         "_circ_uniq\t".join(samples)+"_circ_uniq\t"+"_linear_left\t".join(samples) + "_linear_left\t"+
                         "_linear_right\t".join(samples)+"_linear_right\n")



    if(args.common_names!= None):
        #avigayel get the genes common names
        file_common_names = open(args.common_names,'r')
        gene_common_name={}
        for gene in file_common_names:
            gene = gene.strip('\n').strip('\t').split('\t')
            gene_common_name[gene[1]] = gene[0]

        file_common_names.close()

    sumCirc= 0
    sumLin= 0
    for key in samples_dict:
        # Parsing of the circle name, looking like: gene_name:transcript;chr:locA-locB
        name_loc = key.split(';')
        name = name_loc[0]
        if(args.common_names!=None):
            common_name = gene_common_name[name.split(':')[1]]
        chrm = name_loc[1].split(":")[0]
        locA = name_loc[1].split(":")[1].split("-")[0]
        locB = name_loc[1].split(":")[1].split("-")[1]
        #avigayel add the strand
        strand = name_loc[2]
        if(args.common_names!=None):
            output_samples.write(chrm+"\t"+locA+"\t"+locB+"\t"+name+"|"+common_name+"\t")
        else:
            output_samples.write(chrm+"\t"+locA+"\t"+locB+"\t"+name+"\t")
        # Pre-allocate a list that will contain all stats about a circle.
        list_temp = ['']*len(samples)*FIELDS_PER_SAMPLE
        for sample in range(len(samples)):
            # list_temp[sample] = str(samples_dict[key][sample][0])                 # uniq_circ
            # list_temp[sample+1*len(samples)] = str(samples_dict[key][sample][1])  # circ
            # list_temp[sample+2*len(samples)] = str(samples_dict[key][sample][2])  # linear_left
            # list_temp[sample+3*len(samples)] = str(samples_dict[key][sample][3])  # uniq_linear_left
            # list_temp[sample+4*len(samples)] = str(samples_dict[key][sample][4])  # linear_right
            # list_temp[sample+5*len(samples)] = str(samples_dict[key][sample][5])  # uniq_linear_right
            # list_temp[sample+6*len(samples)] = str(samples_dict[key][sample][6])  # ratio

            list_temp[sample] = str(samples_dict[key][sample][1])                # circ
            list_temp[sample+1*len(samples)] = str(samples_dict[key][sample][0])   # uniq_circ
            list_temp[sample+2*len(samples)] = str(samples_dict[key][sample][2])  # linear_left
            #list_temp[sample+3*len(samples)] = str(samples_dict[key][sample][3])  # uniq_linear_left
            #list_temp[sample+4*len(samples)] = str(samples_dict[key][sample][4])  # linear_right
            list_temp[sample+3*len(samples)] = str(samples_dict[key][sample][4])  # linear_right
            #list_temp[sample+5*len(samples)] = str(samples_dict[key][sample][5])  # uniq_linear_right
            #list_temp[sample+6*len(samples)] = str(samples_dict[key][sample][6])  # ratio
            sumCirc = sumCirc + int(samples_dict[key][sample][1])
            sumLin = sumLin + int(samples_dict[key][sample][2]) + int(samples_dict[key][sample][4])
        output_samples.write(str(sumCirc)+"\t"+str(name_loc[2])+"\t"+str(sumLin)+"\t"+"\t".join(list_temp)+"\n")
        sumCirc= 0
        sumLin= 0
    output_samples.close()

    # Write all sample safeties:
    output_safeties.write('chr\tlocA\tlocB\tname\tstrand\t'+'_maximal_safety'.join(samples)+'_maximal_safety\n')
    all_keys_safeties = reduce(add, [dic.keys() for dic in safety_dict.values()], [])
    for key in all_keys_safeties:
        name_loc = key.split(';')
        name = name_loc[0]
        if(args.common_names!=None):
            common_name = gene_common_name[name.split(':')[1]]
        chrm = name_loc[1].split(":")[0]
        locA = name_loc[1].split(":")[1].split("-")[0]
        locB = name_loc[1].split(":")[1].split("-")[1]

        #avigayel add the strand information
        strand= name_loc[2]
        if(args.common_names!=None):
            output_safeties.write(chrm+"\t"+locA+"\t"+locB+"\t"+name+"|"+common_name+"\t"+strand+"\t")
        else:
            output_safeties.write(chrm+"\t"+locA+"\t"+locB+"\t"+name+"\t"+strand+"\t")
        # Pre-allocate a list that will contain safeties for all samples for a circle
        list_temp = ['']*len(samples)
        for sample in range(len(samples)):
            sample_name = os.path.splitext(os.path.split(samples[sample])[1])[0]
            if key in safety_dict[sample_name]:
                list_temp[sample] = \
                    str(safety_dict[sample_name][key]) if safety_dict[sample_name][key] < sys.maxint else 'INF'
            else:
                list_temp[sample] = 'None'
        output_safeties.write("\t".join(list_temp)+"\n")
    output_safeties.close()

    # Remove tmp folder.
    remove_temp()
    print "-----------------------------"
    print "        Finished"
    print "-----------------------------"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="  _______________________________________\n"
                    " |  SRCP - Short Reads circRNA Pipeline  |\n"
                    " |_______________________________________|",
        formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-circ", "-ci", "-c", help="The circles locations file in BED6 format.", required=True)
    parser.add_argument("-annotation", "-an", "-a", help="The annotation for the species."
                                                         " BED12 format excluding the color column.",
                        required=True)
    parser.add_argument("-read_len", "-len", "-l", help="Minimum length of the read. To use only if dataset contains "
                                                        "reads of different length. DEFAULT=None", required=False, default=None,
                                                        type=int)
    parser.add_argument("-junction", "-junc", "-j", help="Length of single side of junction flanking region."
                                                         " DEFAULT=10", required=False, default=10, type=int)
    parser.add_argument("-system", "-sys", "-s", help="System to use. Path to directory including fasta file of the system's genome and bowtie2 reference. Path with prefix", required=False)
    parser.add_argument("-transcriptome", "-trans", "-tra",
                        help="Bowtie2 reference of multi-fasta files including transcriptome. Path with prefix.", required=True)
    parser.add_argument("-numThread", "-p", help="Number of threads to use. DEFAULT=1", required=False, default=1, type=int)
    parser.add_argument("-data", "-d", "-dat", help="The data files, comma separated.", required=True)
    parser.add_argument("-pcp", help="If given, the pcp junctions output will be used to create linear index."
                                     " DEFAULT=None", required=False, default=None)
    parser.add_argument("-mis", "-m", "-M", help="Number of mismatches allowed in genome and transcriptome alignments."
                                                 " DEAFULT=2", required=False, choices=['0', '1', '2','3'], default='2')
    parser.add_argument("-N", help="Allow N mismatch in seed of the bowtie alignments."
                                   "Should ONLY be used if reads length is < 32bp. DEFAULT = 0",
                        required=False, choices=['0', '1', '2'], default='0')
    parser.add_argument("-safety", "-safe", "-sa", help="The minimal allowed safety. DEFAULT = args.mis+1",
                        required=False, default=None, type=int)
    parser.add_argument("--only_annotated", help="Use only annotated circles (EXACT). DEFAULT=false",
                        action='store_true')
    parser.add_argument("--keep_temporary", "--keep", help="Keep temporary files when done. DEFAULT=false",
                        action='store_true')
    parser.add_argument("-quality", "-qual", "-q", help="The quality standard, DEFAULTS=Phred+33",required=False,
                        choices=[33, 64], default=33)
    parser.add_argument("-min_quality", "-min_qual", "-mq",
                        help="The minimum quality below which reads from dataset are filtered. (DEFAULT=20)",
                        required=False, default=20)
    parser.add_argument("-output", "-o", "-O", help="Sets the prefix of the output files. DEFAULT=outputIndex",
                        required=False, default="outputIndex")
    parser.add_argument("--version", "--v", help="Show the version and about", action='version', version=version)
    #avigayel- add option to give a file with common names
    parser.add_argument("-common_names","-cn",help="optional genes common name file",required=False)
    args = parser.parse_args()

    print "-----------------------------"
    print "       Started SRCP"
    print "-----------------------------"
    print "-----------------------------"
    print "Writing needed shell files..."
    # Here we write a shel lscript file. it is used in the linear index creation
    # to remove duplicate rows and add numbering to the rows.
    # In older versions this script was an external dependency.
    print "-----------------------------"
    print "Trying to create tmp directory..."
    if not os.path.isdir("tmp"):
        print "Directory was created."
        os.mkdir('tmp')
    else:
        print "Directory was present. Rewriting tmp directory."
	os.system('rm -rf tmp')
	os.mkdir('tmp')    
    rdf = open('tmp/remove_duplicates_fa.sh', 'w')
    rdf.write('''
    fname="$1"\n
    echo "Started duplicates removal on "$fname"..."\n
    echo "left processing..."\n
    sed -e '/^>/s/$/@/' -e 's/^>/#/' "tmp/__"$fname"_left__.fa"  |\\\n
    tr -d '\\n' | tr "#" "\\n" | tr "@" "\\t" |\\\n
    sort -u |\\\n
    sed -e 's/^/>/' -e 's/\\t/\\n/' > "tmp/__"$fname"_left_no_dupl__.fa"\n
    echo "Numbering lines..."\n
    awk '{if(NR%2==0)printf "%s;%d\\n", $0, NR/2;else print $0}' < "tmp/__"$fname"_left_no_dupl__.fa" > "tmp/__"$fname"_left_no_dupl_num__.fa"\n
    echo "$(tail -n +2 "tmp/__"$fname"_left_no_dupl_num__.fa")" > "tmp/__"$fname"_left_no_dupl_num__.fa"\n
    echo "right processing..."\n
    sed -e '/^>/s/$/@/' -e 's/^>/#/' "tmp/__"$fname"_right__.fa"  |\\\n
    tr -d '\\n' | tr "#" "\\n" | tr "@" "\\t" |\\\n
    sort -u |\\\n
    sed -e 's/^/>/' -e 's/\\t/\\n/' > "tmp/__"$fname"_right_no_dupl__.fa"\n
    echo "Numbering lines..."\n
    awk '{if(NR%2==0)printf "%s;%d\\n", $0, NR/2;else print $0}' < "tmp/__"$fname"_right_no_dupl__.fa" > "tmp/__"$fname"_right_no_dupl_num__.fa"\n
    echo "$(tail -n +2 "tmp/__"$fname"_right_no_dupl_num__.fa")" > "tmp/__"$fname"_right_no_dupl_num__.fa"\n
    echo "Done!"
    ''')
    rdf.close()
    #shahar added this function
    bed_12_to_11()
    #avigayel - antisenese annotation
    annotate_circles_sense_antisense()

    outputAll = open(args.output+'_circReadsAll.sam', 'w')
    # fix output filename:
    print "Checking output filename..."
    splited = args.output.split('.')
    if len(splited) == 1:
        args.output += ".sam"
    if splited[-1] != 'sam':
        args.output.replace("."+splited[-1], ".sam")
    samples_dict = dict()
    safety_dict = dict()
    print "-----------------------------"
    print "Reading genome..."
    ####################### I removed the default option #####################################
    #if args.system in ['dm3', 'mm9', 'hg19', 'rn6']:
    #    args.system = '/home/kadenerlab/Documents/reference/'+args.system+'/'+args.system
    #else:
    #    pass
    #def_transc = '/home/kadenerlab/Documents/reference/riboIP_con_transctriptom/myGTFdm3_multiFasta_all'
    #if args.system != '/home/kadenerlab/Documents/reference/dm3/dm3' and args.transcriptome == def_transc:
    #    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    #    print "!!!!!!!!There might be an error as the default trancriptome is for DM3!!!!!!!!!!"
    #    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    ############################################################################################
    fasta = open(args.system+'.fa', 'r').read().split('>')
    index = [ chrm.strip('\n').split('\n')[0] for chrm in fasta]
    index = index[1:]
    fasta = fasta[1:]
    fasta = [fasta[i][len(index[i])+1:].replace('\n', '') for i in range(len(fasta))]
    print "-----------------------------"
    print "Extracting circRNA exons..."
    find_exons()
    print "-----------------------------"
    if args.only_annotated:
        args.exons = 'tmp/exons_EXACT.bed'
    else:
        args.exons = 'tmp/exons.bed'
    if args.safety is None:
        args.safety = int(args.mis) + 1
    analyze_samples()
