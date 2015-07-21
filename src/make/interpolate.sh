#!/bin/bash

###
# Interpolate paths to binaries from bins.cfg into target scripts
#

config=$1
buildDir=$2

#read bin locations from cfg.bin
declare -A bins=(
  ["Rbin"]=`grep Rscript $config|cut -d, -f2` 
  # ["BLASTbin"]=`grep blast $config|cut -d, -f2` 
  ["SGEbin"]=`grep qstat $config|cut -d, -f2`
  ["BowtieBin"]=`grep bowtie2 $config|cut -d, -f2`
  # ["JavaBin"]=`grep java $config|cut -d, -f2`
  # ["SNVerBin"]=`grep SNVer $config|cut -d, -f2`
  ["SAMtoolsBin"]=`grep samtools $config|cut -d, -f2`
  ["BCFtoolsBin"]=`grep bcftools $config|cut -d, -f2`
)

function interpolate {
  
  targetfile=$1

  #Check if file exists
  if [ ! -r $targetfile ] 
  then
    echo "Cannot read file $targetfile!"
    exit 1
  fi

  #Load file contents into variable
  contents=`cat $targetfile`

  #iterate over replacements
  for key in ${!bins[@]}
  do
    #escape slash characters (don't even think about messing with this!)
    value=`echo ${bins[$key]}|sed -e 's/\//\\\\\//g'`
    #replace current variable with its value
    contents=`echo "$contents"|sed "s/\\\$$key/$value/g"`
    # echo "$key => ${bins[$key]}"
  done

  #return results
  echo "$contents">$targetfile
}

# interpolate "$buildDir/lib/analyze2.R"
interpolate "$buildDir/lib/master.R"
interpolate "$buildDir/lib/demuxer.R"
interpolate "$buildDir/lib/consolidateAndAlign.R"
interpolate "$buildDir/lib/libyogiseq.R"
