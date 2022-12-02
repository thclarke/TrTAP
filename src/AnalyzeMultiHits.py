#! /usr/bin/env python
# Susan Corbett - Spring 2012 and transformed by Thomas Clarke in 20222
# This script will read the given .csv file, looking for query results with multiple hits. 
# It reports to the .txt file any that do not overlap by more than the specified number of bases.
# Parameters:
#      -csvIn   InputCSVFilename       (required)
#      -txtOut  OutputTXTFilename      (default is MultHits.txt)
#      -overlap NumBasesOverlap        (default is 10)
#      -evalue	Evaluecutoff	       (default is none)


import re
import csv
import sys
import argparse

class Rpt:
  ALWAYS=0
  DEBUG_ONLY=1


def CalcDistBetweenHits(HitLocations, Report, SeqID):
  if len(HitLocations) < 2:
        return -1
    # Check for complete overlap - j sequence is contained within i sequence
  MinOverlap = sys.maxsize
  
    # Check for overlap at the left end
  for i in HitLocations:
      LeftEnd = i[0]
      RightEnd = i[1]
      PrevId = i[2]
      WriteToReport(Report,"LeftEnd=%d RightEnd=%d\n"%(LeftEnd,RightEnd),Rpt.DEBUG_ONLY)
      for j in HitLocations:
        # Don't need to compare to itself
        if i == j:
          continue

        Overlap = 0
        WriteToReport(Report,"j[0]=%d j[1]=%d\n"%(j[0],j[1]),Rpt.DEBUG_ONLY)

        # Check for no overlap at all
        if RightEnd < j[0]:
       #   WriteToReport(Report, "Gap found for %s and %s in %s that is size %d - %d\n"%(PrevId, j[2], SeqID, j[0], RightEnd), Rpt.ALWAYS)
          return (j[0] - RightEnd)
        if LeftEnd > j[1]:
      #    WriteToReport(Report,"Gap found for %s and %s in %s that is size %d - %d\n"%(j[2], PrevId, SeqID, LeftEnd,j[1]), Rpt.ALWAYS)
          return (LeftEnd - j[1])

                                                                      # Check for complete overlap - i sequence is contained within j sequence
        if LeftEnd >= j[0] and RightEnd <= j[1]:
          WriteToReport(Report,"Complete overlap 1\n",Rpt.DEBUG_ONLY)
          Overlap = -2
        elif LeftEnd < j[0] and RightEnd >= j[0]:
          WriteToReport(Report,"Overlap at left end\n",Rpt.DEBUG_ONLY)
          Overlap = -1

                                                    # Check for overlap at the right end
        elif LeftEnd <= j[1] and RightEnd > j[1]:
          WriteToReport(Report,"Overlap at right end\n",Rpt.DEBUG_ONLY)
          Overlap = -1

          WriteToReport(Report,"Overlap=%d\n"%(Overlap),Rpt.DEBUG_ONLY)
        WriteToReport(Report,"MinOverlap=%d\n"%(MinOverlap),Rpt.DEBUG_ONLY)

        if Overlap < MinOverlap:
          MinOverlap = Overlap

  return MinOverlap

# To get debug info, just comment out the if line
def WriteToReport(Report,String,Code):
  
  if (Code == Rpt.ALWAYS):
    Report.write(String)


# Make sure each hit overlaps all others
def CalcMinOverlap(HitLocations,Report):

  if len(HitLocations) < 2:
    return -1
  
  MinOverlap = sys.maxsize
  
  for i in HitLocations:
    LeftEnd = i[0]
    RightEnd = i[1]
    WriteToReport(Report,"LeftEnd=%d RightEnd=%d\n"%(LeftEnd,RightEnd),Rpt.DEBUG_ONLY)
    
    for j in HitLocations:
      # Don't need to compare to itself
      if i == j:
        continue

      Overlap = 0
        
      WriteToReport(Report,"j[0]=%d j[1]=%d\n"%(j[0],j[1]),Rpt.DEBUG_ONLY)

      # Check for no overlap at all
      if RightEnd < j[0] or LeftEnd > j[1]:
        WriteToReport(Report,"No overlap\n",Rpt.DEBUG_ONLY)
        return 0
        
      # Check for complete overlap - i sequence is contained within j sequence
      if LeftEnd >= j[0] and RightEnd <= j[1]:
        WriteToReport(Report,"Complete overlap 1\n",Rpt.DEBUG_ONLY)
        Overlap = RightEnd-LeftEnd+1
        
      # Check for complete overlap - j sequence is contained within i sequence
      elif LeftEnd >= j[0] and RightEnd <= j[1]:
        WriteToReport(Report,"Complete overlap 2\n",Rpt.DEBUG_ONLY)
        Overlap =  j[1]-j[0]+1
        
      # Check for overlap at the left end
      elif LeftEnd < j[0] and RightEnd >= j[0]:
        WriteToReport(Report,"Overlap at left end\n",Rpt.DEBUG_ONLY)
        Overlap = RightEnd-j[0]+1

      # Check for overlap at the right end
      elif LeftEnd <= j[1] and RightEnd > j[1]:
        WriteToReport(Report,"Overlap at right end\n",Rpt.DEBUG_ONLY)
        Overlap = j[1]-LeftEnd+1

      WriteToReport(Report,"Overlap=%d\n"%(Overlap),Rpt.DEBUG_ONLY)
      WriteToReport(Report,"MinOverlap=%d\n"%(MinOverlap),Rpt.DEBUG_ONLY)
      if Overlap < MinOverlap:
        MinOverlap = Overlap
        
  return MinOverlap
# End of CheckForOverlap()


def CoalesceSubjectSeqLengths(HitLocations,Report):
  
  Deleted = False
  WriteToReport(Report,"Before Coalescing:\n",Rpt.DEBUG_ONLY)
  WriteToReport(Report,str(HitLocations)+"\n",Rpt.DEBUG_ONLY)
  
  PrevSubjectSeqID = ""
  
  for i in range(len(HitLocations)):
    if i == 0:
      PrevSubjectSeqID = HitLocations[0][2]
      PrevIndex = 0
      continue
    
#   if i[2] == PrevSubjectSeqID:
#     i[0] = min(HitLocations[PrevIndex][0],i[0])
#     i[1] = max(HitLocations[PrevIndex][1],i[1])
    WriteToReport(Report,"HitLocations[PrevIndex+1][2] is %s\n"%(HitLocations[PrevIndex+1][2]),Rpt.DEBUG_ONLY)
    WriteToReport(Report,"PrevSubjectSeqID is %s\n"%(PrevSubjectSeqID),Rpt.DEBUG_ONLY)
    if HitLocations[PrevIndex+1][2] == PrevSubjectSeqID:
      HitLocations[PrevIndex+1][0] = min(HitLocations[PrevIndex][0],HitLocations[PrevIndex+1][0])
      HitLocations[PrevIndex+1][1] = max(HitLocations[PrevIndex][1],HitLocations[PrevIndex+1][1])
      del HitLocations[PrevIndex]
      WriteToReport(Report,"Deleted an entry\n",Rpt.DEBUG_ONLY)
      Deleted = True
    else:
      PrevSubjectSeqID = HitLocations[PrevIndex+1][2]
      PrevIndex += 1

  if Deleted == True:
    WriteToReport(Report,"*****DELETION WAS MADE*****\nAfter Coalescing:\n",Rpt.DEBUG_ONLY)
    WriteToReport(Report,str(HitLocations)+"\n",Rpt.DEBUG_ONLY)
# End of CoalesceSubjectSeqLengths()


def AnalyzeMultiHits(InCSVFile, OutTXTFile, NumBasesOverlap, eValcutoff):

  BlastFile = open(InCSVFile, 'rU')
  BlastResults = BlastFile.readlines()
  Report = open(OutTXTFile, 'w')

  NumLines = 0
  NumUniqueQueries = 0
  PrevQuerySeqID = ""
  StartInQuery = 0
  EndInQuery = 0
  NumQueryHits = 0
  HitLocations = []
  

  for BRL in BlastResults:
    BRLine = BRL.split("\t")
    NumLines += 1
	
    QuerySeqID = BRLine[0]
    SubjectSeqID = BRLine[1]
    if float(BRLine[10]) <= eValcutoff:
       WriteToReport(Report,"QuerySeqID is %s\n"%(QuerySeqID),Rpt.DEBUG_ONLY)
       WriteToReport(Report,"SubjectSeqID is %s\n"%(SubjectSeqID),Rpt.DEBUG_ONLY)
       if QuerySeqID != PrevQuerySeqID:
          if NumQueryHits > 1:
            CoalesceSubjectSeqLengths(HitLocations,Report)
	    if (len(HitLocations) > 1):
              WriteToReport(Report,"\nCalling CalcMinOverlap() - NumQueryHits=%d\n"%(NumQueryHits),Rpt.DEBUG_ONLY)
              MinOverlap = CalcMinOverlap(HitLocations,Report)
              WriteToReport(Report,"NumQueryHits=%d\n"%(NumQueryHits),Rpt.DEBUG_ONLY)
              WriteToReport(Report,"MinOverlap=%d\n"%(MinOverlap),Rpt.DEBUG_ONLY)
              WriteToReport(Report,"NumBasesOverlap=%d\n"%(NumBasesOverlap),Rpt.DEBUG_ONLY)
              if MinOverlap < NumBasesOverlap:
                Gap = CalcDistBetweenHits(HitLocations, Report,PrevQuerySeqID)
                WriteToReport(Report,"%s\n"%(PrevQuerySeqID),Rpt.ALWAYS)
          #else:
          #  WriteToReport(Report,"Sufficient Overlap: Query=%s Overlap=%d\n"%(PrevQuerySeqID, MinOverlap),Rpt.ALWAYS)
          NumUniqueQueries += 1
          PrevQuerySeqID = QuerySeqID
          NumQueryHits = 0
        # Clear out the list
          del HitLocations[0:len(HitLocations)]

       if (int(BRLine[6]) < int(BRLine[7])):
         StartInQuery = int(BRLine[6])
         EndInQuery = int(BRLine[7])
       else:
         StartInQuery = int(BRLine[7])
         EndInQuery = int(BRLine[6])
       WriteToReport(Report,"StartInQuery=%d EndInQuery=%d\n"%(StartInQuery,EndInQuery),Rpt.DEBUG_ONLY)
       HitLocations.append([StartInQuery,EndInQuery,SubjectSeqID])
       NumQueryHits += 1
  WriteToReport(Report,"%d Total Lines Processed\n"%(NumLines),Rpt.DEBUG_ONLY)
  WriteToReport(Report,"%d Unique Queries Found\n"%(NumUniqueQueries),Rpt.DEBUG_ONLY)
# End of AnalyzeMultiHits()



def main():

  parser = argparse.ArgumentParser(description='Reads the given .csv file, looking for query results with multiple hits, and reports to the .txt file any that do not overlap by more than the specified number of bases.')
  parser.add_argument('-csvIn', dest='InCSVFile', action='store', required=True,
                     help='Specify the .csv input filename')
  parser.add_argument('-txtOut', default="MultHits.txt", 
                      dest='OutTXTFile', action='store',
                     help='Specify the .txt output filename (default: MultHits.txt)')
  parser.add_argument('-overlap', default=50,
                      dest='NumBasesOverlap', action='store', type=int,
                     help='Specify the required overlap between hit sequences. Default = 50.')
  parser.add_argument('-evalue', default=1e-5,
		     dest='eValcutoff', action='store', type=float,
		     help='Specify the evalue cutoff for hit to be used. Default = 1e-5')

  args = parser.parse_args()

  AnalyzeMultiHits(args.InCSVFile, args.OutTXTFile, args.NumBasesOverlap, args.eValcutoff)
# End of main()



if __name__ == '__main__':
    main()

